#include <iostream>
#include <cmath>
#include <array>
#include <type_traits>
#include <cstdlib>

template <typename T>
constexpr bool is_pow2(T n)
{
   return (n & (n - 1)) == 0;
}

template <std::uint32_t N, typename T = std::uint32_t>
struct bitstream
{
   static_assert(is_pow2(N), "N must be a power of 2, except 0");
   static_assert(std::is_unsigned<T>::value, "T must be unsigned");

   static constexpr auto nbits = 8 * sizeof(T);
   static constexpr auto array_size = N / nbits;

   void clear()
   {
      bits.fill(0);
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

   static std::uint32_t count_bits(std::uint32_t i)
   {
      // GCC only!!!
      return __builtin_popcount(i);
   }

   static std::uint64_t count_bits(std::uint64_t i)
   {
      // GCC only!!!
      return __builtin_popcountll(i);
   }

   template <typename F>
   void auto_correlate(F f)
   {
      // The first will always be zero:
      f(0, 0);

      constexpr auto mid_array = (array_size / 2) - 1;
      constexpr auto mid_pos = N / 2;

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

   std::array<T, array_size> bits;
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
   constexpr auto sps = 20000;   //  Samples per second
   constexpr auto buff_size = 1024;

   ////////////////////////////////////////////////////////////////////////////
   // Generate a test signal
   std::array<float, buff_size> signal;
   noise ns; // noise

   const float f = 82.41;
   float p  = float(sps) / f;

   for (int i = 0; i < buff_size; i++)
   {
      signal[i] = 0.1 * ns();                   // Noise
      signal[i] += 0.3 *  sin(2 * pi * i / p);  // First harmonic
      signal[i] += 0.4 *  sin(4 * pi * i / p);  // Second harmonic
      signal[i] += 0.3 *  sin(6 * pi * i / p);  // Third harmonic
   }

   ////////////////////////////////////////////////////////////////////////////
   // Zero crossing
   zero_cross zc;

   bitstream<buff_size> bin;
   bin.clear();

   for (auto i = 0; i != buff_size; ++i)
      bin.set(i, zc(signal[i]));

   ////////////////////////////////////////////////////////////////////////////
   // Binary Auto-correlation
   std::array<std::uint32_t, buff_size / 2> corr;
   bin.auto_correlate(
      [&corr](auto pos, auto count)
      {
         corr[pos] = count;
      }
   );

   ////////////////////////////////////////////////////////////////////////////
   // Print the signal, zero crossings and correlations
   constexpr auto mid_pos = buff_size / 2;
   for (int i = 0; i < buff_size; i++)
   {
      if (i < mid_pos)
         std::cout << signal[i] << ", " << bin.get(i) << ", " << float(corr[i])/mid_pos << std::endl;
      else
         std::cout << signal[i] << ", " << bin.get(i) << std::endl;
   }

   return 0;
}
