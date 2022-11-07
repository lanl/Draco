//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Endian.cc
 * \author Kelly Thompson
 * \date   Wed Nov 09 14:15:14 2011
 * \brief  Function declarations for endian conversions
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Endian.hh"
#include <array>

namespace rtt_dsxx {
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Does this platform use big or little endianness
 * \return true if platform uses big endian format
 */
bool is_big_endian() {
  union {
    uint32_t i;
    std::array<char, 4> c;
  } data = {0x01020304};

  return data.c[0] == 1;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Does this platform support IEEE float representation?
 *
 * Some older Cray machines did not support the IEEE float representation.  This simple test will
 * identify machines that are IEEE compliant.
 *
 * \return true if we support IEEE float representation.
 */
bool has_ieee_float_representation() {
  // start by assume IEEE platform (i.e.: not a Cray machine).
  bool i_am_ieee(true);

  // Create a double precision value that will be used to test bit representations.
  double d_two(2.0);
  size_t constexpr size(sizeof(double));
  // Generate a bit-by-bit view of the double precision value:
  std::array<char, size> char_two{'\0'};
  std::memcpy(&char_two, &d_two, size);

  // IEEE reference value:
  std::array<char, size> ieee64_two = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40};
  if (is_big_endian())
    char_byte_swap(&ieee64_two[0], size);

  // Cray reference value:
  //
  // Note The 5th value of the actual Cray representation causes overflow warnings on IEEE machines:
  // char cray64_two[size] = {0x00,0x00,0x00,0x00,0x00,0x80,0x02,0x40};

  for (size_t i = 0; i < size; ++i)
    if (char_two[i] != ieee64_two[i]) {
      i_am_ieee = false;
      break;
    }

  return i_am_ieee;
}

} // end namespace rtt_dsxx

//! These versions can be called by Fortran.  They wrap the C++ implementation.
extern "C" {

int dsxx_is_big_endian() {
  if (rtt_dsxx::is_big_endian())
    return 1;
  return 0;
}
void dsxx_byte_swap_int(int &value) {
  rtt_dsxx::byte_swap(value);
  return;
}
void dsxx_byte_swap_int64_t(int64_t &value) {
  rtt_dsxx::byte_swap(value);
  return;
}
void dsxx_byte_swap_double(double &value) {
  rtt_dsxx::byte_swap(value);
  return;
}

} // end extern "C"

//------------------------------------------------------------------------------------------------//
// end of ds++/Endian.cc
//------------------------------------------------------------------------------------------------//
