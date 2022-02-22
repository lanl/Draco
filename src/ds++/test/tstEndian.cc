//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstEndian.cc
 * \author Mike Buksas
 * \date   Tue Oct 23 16:20:59 2007
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/Endian.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <array>
#include <limits>
#include <sstream>

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void test_char_data(ScalarUnitTest &ut) {
  array<unsigned char, 3> data = {'a', 'b', 'c'};
  unsigned int length = sizeof(data) / sizeof(unsigned char);

  char_byte_swap(data.data(), length);

  if ((data[0] != 'c') || (data[1] != 'b') || (data[2] != 'a'))
    FAILMSG("unsigned char_byte_swap function failed");

  array<char, 3> pdata = {'a', 'b', 'c'};
  unsigned int plength = sizeof(pdata) / sizeof(char);

  char_byte_swap(pdata.data(), plength);

  if ((pdata[0] != 'c') || (pdata[1] != 'b') || (pdata[2] != 'a'))
    FAILMSG("plain char_byte_swap function failed");
}

//------------------------------------------------------------------------------------------------//
void test_integer(ScalarUnitTest &ut) {
  // Integer. This value overflows unsigned ints.
  int moo = 0xDEADBEEF;

  byte_swap(moo);

  if (static_cast<unsigned>(moo) != 0xEFBEADDE)
    FAILMSG("byte_swap failed for for integer type");

  // Unsigned integer
  unsigned int u_moo = 0xDEADBEEF;

  byte_swap(u_moo);

  if (u_moo != 0xEFBEADDE)
    FAILMSG("byte_swap failed for for unsigned integer type");

  // uint32_t, to test the specialized version of byte_swap_copy
  uint32_t uint32_moo = 0xDEADBEEF;

  uint32_t uint32_moo_swapped = byte_swap_copy(uint32_moo);

  if (uint32_moo_swapped != 0xEFBEADDE)
    FAILMSG("byte_swap_copy failed for for unsigned integer type");
}

//------------------------------------------------------------------------------------------------//
void test_int64(ScalarUnitTest &ut) {
  // Integer.
  int64_t moo = 0xFADEDDEADBEEFBAD;

  byte_swap(moo);

  if (static_cast<uint64_t>(moo) != 0xADFBEEDBEADDDEFA)
    FAILMSG("byte_swap failed for for int64 type");

  // Unsigned integer
  uint64_t u_moo = 0xFADEDDEADBEEFBAD;

  byte_swap(u_moo);

  if (u_moo != 0xADFBEEDBEADDDEFA)
    FAILMSG("byte_swap failed for for uint64 integer type");

  byte_swap(u_moo);
  if (u_moo != 0xFADEDDEADBEEFBAD)
    FAILMSG("2x byte_swap failed for for uint64 integer type");

  // Swap again, using byte_swap_copy
  uint64_t u_moo_swapped = byte_swap_copy(u_moo);

  if (u_moo_swapped != 0xADFBEEDBEADDDEFA)
    FAILMSG("byte_swap_copy failed for for uint64 integer type");
}

//------------------------------------------------------------------------------------------------//
/* This test demonstrates that two applications of byte-swap in succession
 * return the original value.
 *
 * To do this, we sweep over a lot of double values. We use a non-integral
 * multiplier for successive values to avoid small subsets of the available
 * patterns of bits. E.g. multiples of 2.
 */
void test_idempotence(ScalarUnitTest &ut) {

  for (double value = 1.0;                               // NOLINT
       value < std::numeric_limits<double>::max() / 4.0; // divide by 4 to
       value *= 3.4)                                     // prevent overflow
  {
    // Use the in-place version to test positive values.
    double local = value;
    byte_swap(local);
    byte_swap(local);

    // These numbers should be identical, so I'm testing for equality.
    if (std::abs(local - value) > std::numeric_limits<double>::epsilon())
      FAILMSG("byte_swap failed to reproduce original number");

    // Use the copy-generating version to test negative numbers.
    const double neg_local = byte_swap_copy(byte_swap_copy(-value));

    if (std::abs(neg_local + value) > std::numeric_limits<double>::epsilon())
      FAILMSG("byte_swap failed to reproduce original number");
  }

  return;
}

//------------------------------------------------------------------------------------------------//
void test_ieee_float(ScalarUnitTest &ut) {
  // These tests always pass, but they may print different messages.

  // Endianess
  if (is_big_endian())
    PASSMSG("This machine uses big endian byte ordering.");
  else
    PASSMSG("This machine uses little endian byte ordering.");

  // IEEE floating point?
  if (has_ieee_float_representation()) {
    std::ostringstream msg;
    msg << "Looks like we are on a platform that supports IEEE "
        << "floating point representation.";
    PASSMSG(msg.str());
  } else {
    std::ostringstream msg;
    msg << "This platform does not support IEEE floating point "
        << "representation.";
    PASSMSG(msg.str());
  }
}

//------------------------------------------------------------------------------------------------//
void test_externc(ScalarUnitTest &ut) {
  int result(42);
  result = dsxx_is_big_endian();
  if (result < 0 || result > 1)
    ITFAILS;

  result = 0xDEADBEEF;
  dsxx_byte_swap_int(result);
  if (result != static_cast<int>(0xEFBEADDE))
    ITFAILS;

  int64_t i64(0xFADEDDEADBEEFBAD);
  dsxx_byte_swap_int64_t(i64);
  if (i64 != static_cast<int64_t>(0xADFBEEDBEADDDEFA))
    ITFAILS;

  // Hexedecimal floating-point
  double d(42);
  dsxx_byte_swap_double(d);
  // should not be 42
  if (rtt_dsxx::soft_equiv(d, 42.0, 1.0e-15))
    ITFAILS;
  dsxx_byte_swap_double(d);
  // double swap should return 42
  if (!rtt_dsxx::soft_equiv(d, 42.0, 1.0e-15))
    ITFAILS;

  return;
}
//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, release);
  try {
    test_char_data(ut);
    test_integer(ut);
    test_int64(ut);
    test_idempotence(ut);
    test_ieee_float(ut);
    test_externc(ut);
    PASSMSG("Just Because.");
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstEndian.cc
//------------------------------------------------------------------------------------------------//
