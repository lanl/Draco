//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Endian.hh
 * \author Mike Buksas
 * \date   Tue Oct 23 14:15:55 2007
 * \brief  Function declarations for endian conversions
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef dsxx_Endian_hh
#define dsxx_Endian_hh

#include "ds++/config.h"
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>

//------------------------------------------------------------------------------------------------//
/*!
 * Endian conversion functions.
 *
 * The endian nature of a data representation describes the order in which the constitutent bytes of
 * a multi-byte data value are ordered. We are concerned with converting between big and little
 * endian orderings on platforms where the char data type is one byte in size.
 *
 * If there are other endians out there, I seriously do not want to know about them.
 *
 * To convert between big and little endian data we intrepret the data to be converted as a
 * character array by casting a pointer to the data to (char*). We then manipulate the order, but
 * not the contents, of the character data.
 *
 * Note that we are implicitly assuming that the size of char on each platform is one byte.
 *
 * In order for these functions to work on floating point data, we are assuming that the floating
 * point representations are identical on the two architectures _except_ for the difference in
 * endianness. Also, the sign and exponent information of the floating point representation must fit
 * within a single byte of data, so that it does not require extra steps at the bit-level for
 * conversion.
 */
//------------------------------------------------------------------------------------------------//

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Elemetary byte-swapping routine.
 *
 * \arg The data to byte-swap, represented as character data.
 * \arg The size of the data array.
 *
 * This is a core routine used by other functions to convert data between endian representations.
 *
 * It swaps the elements of a character array of length n. Element 0 is swapped with element n, 1
 * with n-1 etc... The contents of the individual elements are not changed, only their order.
 *
 * For example, consider the unsigned integer value: \c 0xDEADBEEF.  (\c 0x means this is a
 * hexidecimal value) Two hexidecimal digits is a single byte (16^2 = 2^8) so the layout of the
 * value in big endian style is:
 * \verbatim
 *       0        1        2        3
 *     D  E     A  D     B  E     E  F
 *  |--------|--------|--------|--------|
 *       ^        ^        ^        ^
 *       |        +--------+        |
 *       +--------------------------+
 *                 swapped
 * \endverbatim
 * The conversion to little endian involves the swap operations pictured in the diagram above. The
 * resulting value (if still interpreted as big-endian) is \c 0xEFBEADDE.
 *
 * We provide two versions for signed and unsigned character data. Internally, we use
 * unsigned. Certain applications use signed char data, and the second form is provided if they need
 * to manipulate the character data directly, instead of using one of the byte_swap functions.
 */
inline void char_byte_swap(unsigned char *data, int n) {
  unsigned char *end = data + n - 1;
  while (data < end)
    std::swap(*data++, *end--);
}

inline void char_byte_swap(char *data, int n) {
  char *end = data + n - 1;
  while (data < end)
    std::swap(*data++, *end--);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief General byte-swapping routine
 *
 * This function operates in place on its argument.
 */

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#if (DBS_GNUC_VERSION >= 70000)
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif
#endif

template <typename T> void byte_swap(T &value) {
  char_byte_swap((unsigned char *)(&value), sizeof(T));
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

//------------------------------------------------------------------------------------------------//
/*!
 * \brief General byte-swapping routine.
 *
 * This function returns a byte-swapped copy of the argument.
 */
template <typename T> T byte_swap_copy(T value) {
  byte_swap(value);
  return value;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Does this platform use big or little endianness
 * \return true if platform uses big endian format
 */
bool is_big_endian();

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Does this platform support IEEE float representation?
 *
 * Some older Cray machines did not support the IEEE float representation.  This simple test will
 * identify machines that are IEEE compliant.
 *
 * \return true if we support IEEE float representation.
 */
bool has_ieee_float_representation();

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Tim Kelley's specialized byte-swapping routines from ds++/swap.hh.
 *
 * \param[in] input Value to be considered for byte-swapping.
 * \return byte-swapped value.
 *
 * Do byte-swapping one of two ways: either use GNU extended asm for x86 (really 486+), or use the
 * "poor people's" method of digging out one byte at a time and moving it to the right place. The
 * poor method is really not that bad: GCC for example seems to optimize it down to about 5
 * instructions in the 32 bit case, including loads, versus 2 instructions (including load) for the
 * inline asm case.
 */
template <> inline uint32_t byte_swap_copy<uint32_t>(uint32_t const input) {
#ifdef __use_x86_gnu_asm
  uint32_t output = input;
  asm("bswap %0" : "+g"(output) :);
#else
  uint32_t byte, output;
  byte = input & 255U;
  output = (byte << 24U);

  byte = input & 65280U; // 255 << 8
  output = output | (byte << 8U);

  byte = input & 16711680U; // 255 << 16
  output = output | (byte >> 8U);

  byte = input & 4278190080U;      // 255 << 24
  output = output | (byte >> 24U); // look out--algebraic shift r.

#endif // __use_x86_asm
  return output;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Tim Kelley's specialized byte-swapping routines from ds++/swap.hh.
 *
 * \param[in] input Value to be considered for byte-swapping.
 * \return byte-swapped value.
 *
 * Do byte-swapping one of two ways: either use GNU extended asm for x86 (really 486+), or use the
 * "poor people's" method of digging out one byte at a time and moving it to the right place. The
 * poor method is really not that bad: GCC for example seems to optimize it down to about 5
 * instructions in the 32 bit case, including loads, versus 2 instructions (including load) for the
 * inline asm case.
 */
template <> inline double byte_swap_copy<double>(double const input) {
#ifdef __use_x86_gnu_asm
  double output = input;
  asm("bswap %0" : "+g"(output) :);
#else
  union {
    double d;
    uint64_t u;
  } b64 = {0.0};

  uint64_t byte{0};
  uint64_t tmp{0};
  uint64_t uinput{0};

  // change meaning of input bits to uint64_t:
  b64.d = input;
  uinput = b64.u;

  // 1
  byte = uinput & 255U;
  tmp = (byte << 56U);
  // 2
  byte = uinput & 65280U; // 255 << 8
  tmp = tmp | (byte << 40U);
  // 3
  byte = uinput & 16711680U; // 255 << 16
  tmp = tmp | (byte << 24U);
  // 4
  byte = uinput & 4278190080U; // 255 << 24
  tmp = tmp | (byte << 8U);
  // 5
  byte = uinput & 1095216660480ULL; // 255 << 32
  // byte = uinput & 1095216660480; // 255 << 32
  tmp = tmp | (byte >> 8U);
  // 6
  byte = uinput & 280375465082880ULL; // 255 << 40
  tmp = tmp | (byte >> 24U);
  // 7
  byte = uinput & 71776119061217280ULL; // 255 << 48
  tmp = tmp | (byte >> 40U);
  // 8
  byte = uinput & 18374686479671623680ULL; // 255 << 56
  // byte = uinput & 0xff00000000000000; // 255 << 56
  tmp = tmp | (byte >> 56U);

  // change meaning of bits in b64.
  b64.u = tmp;
  double output = b64.d;
#endif // __use_x86_gnu_asm

  return output;
}

} // end namespace rtt_dsxx

//! These versions can be called by Fortran.  They wrap the C++ implementation.
extern "C" {
int dsxx_is_big_endian();
void dsxx_byte_swap_int(int &value);
void dsxx_byte_swap_int64_t(int64_t &value);
void dsxx_byte_swap_double(double &value);
}

#endif // dsxx_Endian_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/Endian.hh
//------------------------------------------------------------------------------------------------//
