#include <iostream>
#include <cmath>
#include <array>
#include <type_traits>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cstdint>

// smallest power of 2 that fits n
template <typename T>
constexpr T smallest_pow2(T n, T m = 1)
{
   return (m < n)? smallest_pow2(n, m << 1) : m;
}

std::uint32_t count_bits(std::uint32_t i)
{
   // GCC only!!!
   return __builtin_popcount(i);
}

std::uint64_t count_bits(std::uint64_t i)
{
   // GCC only!!!
   return __builtin_popcountll(i);
}

template <typename T = std::uint32_t>
struct bitstream
{
   static_assert(std::is_unsigned<T>::value, "T must be unsigned");
   static constexpr auto nbits = 8 * sizeof(T);

   bitstream(std::size_t size_)
   {
      size = smallest_pow2(size_);
      array_size = size / nbits;
      bits.resize(array_size, 0);
   }

   void clear()
   {
      std::fill(bits.begin(), bits.end(), 0);
   }

   void set(std::uint32_t i, bool val)
   {
      auto mask = 1 << (i % nbits);
      auto& ref = bits[i / nbits];
      ref ^= (-T(val) ^ ref) & mask;
   }

   bool get(std::uint32_t i) const
   {
      auto mask = 1 << (i % nbits);
      return (bits[i / nbits] & mask) != 0;
   }

   template <typename F>
   void auto_correlate(F f)
   {
      // The first will always be zero:
      f(0, 0);

      auto mid_array = (array_size / 2) - 1;
      auto mid_pos = size / 2;
      auto index = 0;
      auto shift = 1;

      for (auto pos = 1; pos != mid_pos; ++pos)
      {
         auto* p1 = bits.data();
         auto* p2 = bits.data() + index;
         auto count = 0;

         if (shift == 0)
         {
            for (auto i = 0; i != mid_array; ++i)
               count += count_bits(*p1++ ^ *p2++);
         }
         else
         {
            auto shift2 = nbits - shift;
            for (auto i = 0; i != mid_array; ++i)
            {
               auto v = *p2++ >> shift;
               v |= *p2 << shift2;
               count += count_bits(*p1++ ^ v);
            }
         }
         ++shift;
         if (shift == nbits)
         {
            shift = 0;
            ++index;
         }

         f(pos, count);
      }
   }

   std::vector<T> bits;
   std::size_t size;
   std::size_t array_size;
};

struct zero_cross
{
   bool operator()(float s)
   {
      if (s < -0.1f)
         y = 0;
      else if (s > 0.0f)
         y = 1;
      return y;
   }

   bool y = 0;
};

struct noise
{
   float operator()() const
   {
      return (float(rand()) / (RAND_MAX / 2)) - 1.0;
   }
};

int main ()
{
   constexpr auto pi = M_PI;
   constexpr auto sps = 44100;                     //  Samples per second

   ////////////////////////////////////////////////////////////////////////////
   // Generate a test signal

   constexpr float freq = 440; // 82.41;
   constexpr float period  = float(sps) / freq;    // period
   constexpr float noise_level = 0.0;              // Noise level
   constexpr float _1st_level = 0.3;               // First harmonic level
   constexpr float _2nd_level = 0.4;               // Second harmonic level
   constexpr float _3rd_level = 0.3;               // Third harmonic level

   std::size_t buff_size = smallest_pow2<std::size_t>(std::ceil(sps / freq)) * 2;

   std::vector<float> signal(buff_size);
   noise ns; // noise

   for (int i = 0; i < buff_size; i++)
   {
      signal[i] = 0; // noise_level * ns();                       // Noise
      signal[i] += _1st_level *  sin(2 * pi * i / period);  // First harmonic
      signal[i] += _2nd_level *  sin(4 * pi * i / period);  // Second harmonic
      signal[i] += _3rd_level *  sin(6 * pi * i / period);  // Third harmonic
   }

   ////////////////////////////////////////////////////////////////////////////
   // The bitstream

   bitstream<> bin(buff_size);

   ////////////////////////////////////////////////////////////////////////////
   // Zero crossing
   zero_cross zc;

   for (auto i = 0; i != buff_size; ++i)
      bin.set(i, zc(signal[i]));

   ////////////////////////////////////////////////////////////////////////////
   // Binary Auto-correlation

   std::uint32_t max_count = 0;
   std::uint32_t min_count = UINT32_MAX;
   std::size_t est_index = 0;
   std::vector<std::uint32_t> corr(buff_size / 2);
   bin.auto_correlate(
      [&corr, &max_count, &min_count, &est_index](auto pos, auto count)
      {
         corr[pos] = count;
         max_count = std::max<std::uint32_t>(max_count, count);
         if (pos && count < min_count)
         {
            min_count = count;
            est_index = pos;
         }
      }
   );

   ////////////////////////////////////////////////////////////////////////////
   // Print the signal, zero crossings and correlations
   // auto mid_pos = buff_size / 2;
   // for (int i = 0; i < buff_size; i++)
   // {
   //    if (i < mid_pos)
   //       std::cout << signal[i] << ", " << bin.get(i) << ", " << float(corr[i])/max_count << std::endl;
   //    else
   //       std::cout << signal[i] << ", " << bin.get(i) << std::endl;
   // }

   ////////////////////////////////////////////////////////////////////////////
   // Estimate the pitch

   // float mid   = corr[est_index];
   // float left  = corr[est_index-1];
   // float right = corr[est_index+1];

   // //  assert( 2*mid - left - right > 0.0 );

   // float shift = 0.5 * (right-left) / (2 * mid - left - right);
   // float est_freq = est_index + shift;

/////////////////////////////////////
   float prev = 0;
   auto p = signal.begin();
   for (; *p <= 0.0f; ++p)
      prev = *p;
   auto dy = *p - prev;
   auto dx1 = -prev / dy;

   auto p2 = signal.begin() + est_index - 1;
   for (; *p2 <= 0.0f; ++p2)
      prev = *p2;
   dy = *p2 - prev;
   auto dx2 = -prev / dy;

   float samples_period = (p2-p) + (dx2 - dx1);
   float est_freq = sps / samples_period;

   std::cout << "Actual Frequency: " << freq << " Hz" << std::endl;
   std::cout << "Estimated Frequency: " << est_freq << " Hz" << std::endl;
   std::cout << "Error: " << 1200.0 * std::log2(est_freq / freq) << " cents" << std::endl;

   return 0;
}
