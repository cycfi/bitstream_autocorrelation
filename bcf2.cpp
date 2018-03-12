/*=============================================================================
   Copyright (c) 2018 Cycfi Research. All rights reserved.

   Distributed under the MIT License [ https://opensource.org/licenses/MIT ]
=============================================================================*/
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
   void auto_correlate(std::size_t start_pos, F f)
   {
      auto mid_array = (array_size / 2) - 1;
      auto mid_pos = size / 2;
      auto index = start_pos / nbits;
      auto shift = start_pos % nbits;

      for (auto pos = start_pos; pos != mid_pos; ++pos)
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
   constexpr auto sps = 44100;               // 20000;
   constexpr auto min_freq = 50.0;
   constexpr auto max_freq = 500.0;
   constexpr float freq = 261.626;           // 82.41;

   // These are in samples
   constexpr float period = float(sps) / freq;
   constexpr float min_period = float(sps) / max_freq;
   constexpr float max_period = float(sps) / min_freq;

   ////////////////////////////////////////////////////////////////////////////
   // Generate a test signal

   constexpr float noise_level = 0.0;     // Noise level (dB)
   constexpr float _1st_level = 0.3;      // First harmonic level
   constexpr float _2nd_level = 0.4;      // Second harmonic level
   constexpr float _3rd_level = 0.3;      // Third harmonic level

   constexpr float offset = 0;
   // constexpr float offset = period - 1.3; // Initial offset (some odd number)
   std::size_t buff_size = smallest_pow2<std::size_t>(std::ceil(max_period)) * 2;

   std::vector<float> signal(buff_size);
   noise ns; // noise

   for (int i = 0; i < buff_size; i++)
   {
      auto angle = (i + offset) / period;
      signal[i] = noise_level * ns();                       // Noise
      signal[i] += _1st_level *  std::sin(2 * pi * angle);  // First harmonic
      signal[i] += _2nd_level *  std::sin(4 * pi * angle);  // Second harmonic
      signal[i] += _3rd_level *  std::sin(6 * pi * angle);  // Third harmonic
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

#define PRINTING 0

   std::uint32_t max_count = 0;
   std::uint32_t min_count = UINT32_MAX;
   std::size_t est_index = 0;
   std::vector<std::uint32_t> corr(buff_size / 2);
   bin.auto_correlate(PRINTING? 0 : min_period,
      [&corr, &max_count, &min_count, &est_index](auto pos, auto count)
      {
         corr[pos] = count;
         max_count = std::max<std::uint32_t>(max_count, count);
         if (count < min_count)
         {
            min_count = count;
            est_index = pos;
         }
      }
   );

   ////////////////////////////////////////////////////////////////////////////
   // Print the signal, zero crossings and correlations (for graphing)
#if PRINTING
   auto mid_pos = buff_size / 2;
   for (int i = 0; i < buff_size; i++)
   {
      if (i < mid_pos)
         std::cout << signal[i] << ", " << bin.get(i) << ", " << float(corr[i])/max_count << std::endl;
      else
         std::cout << signal[i] << ", " << bin.get(i) << std::endl;
   }

   return 0; // return now
#endif

   ////////////////////////////////////////////////////////////////////////////
   // Handle harmonics
   auto sub_threshold = 0.1 * max_count;
   int max_div = est_index / min_period;
   for (int div = max_div; div != 0; div--)
   {
      bool all_strong = true;
      float mul = 1.0f / div;

      for (int k = 1; k != div; k++)
      {
         int sub_period = k * est_index * mul;
         if (corr[sub_period] > sub_threshold)
         {
            all_strong = false;
            break;
         }
      }

      if (all_strong)
      {
         est_index = est_index * mul;
         break;
      }
   }

   ////////////////////////////////////////////////////////////////////////////
   // Estimate the pitch

   float prev = 0;
   auto first = signal.begin();
   for (; *first <= 0.0f; ++first)
      prev = *first;
   auto dy = *first - prev;
   auto dx1 = -prev / dy;

   auto last = signal.begin() + est_index - 1;
   for (; *last <= 0.0f; ++last)
      prev = *last;
   dy = *last - prev;
   auto dx2 = -prev / dy;

   float samples_period = (last-first) + (dx2 - dx1);
   float est_freq = sps / samples_period;

   std::cout << "Actual Frequency: " << freq << " Hz" << std::endl;
   std::cout << "Estimated Frequency: " << est_freq << " Hz" << std::endl;
   std::cout << "Error: " << 1200.0 * std::log2(est_freq / freq) << " cents" << std::endl;

   return 0;
}
