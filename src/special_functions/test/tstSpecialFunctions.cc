//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/test/tstSpecialFunctions.cc
 * \author Kelly Thompson
 * \date   Tue Sep 27 12:49:39 2005
 * \brief  Unit tests for kronecker_delta and factorial.
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "special_functions/Factorial.hh"
#include "special_functions/KroneckerDelta.hh"

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstKdelta(rtt_dsxx::UnitTest &ut) {
  using rtt_sf::kronecker_delta;
  if (kronecker_delta(0, 0) == 1)
    PASSMSG("Found kronecker_delta(0,0) == 1, kronecker_delta is working.");
  else
    FAILMSG("Found kronecker_delta(0,0) != 1, kronecker_delta is not working.");
  if (kronecker_delta(0, 1) == 0)
    PASSMSG("Found kronecker_delta(0,1) == 0, kronecker_delta is working.");
  else
    FAILMSG("Found kronecker_delta(0,1) != 0, kronecker_delta is not working.");
  if (kronecker_delta(1, 1) == 1)
    PASSMSG("Found kronecker_delta(1,1) == 1, kronecker_delta is working.");
  else
    FAILMSG("Found kronecker_delta(1,1) != 1, kronecker_delta is not working.");
  if (kronecker_delta(1, 0) == 0)
    PASSMSG("Found kronecker_delta(1,0) == 0, kronecker_delta is working.");
  else
    FAILMSG("Found kronecker_delta(1,0) != 0, kronecker_delta is not working.");
  if (kronecker_delta(-1, 0) == 0)
    PASSMSG("Found kronecker_delta(-1,0) == 0, kronecker_delta is working.");
  else
    FAILMSG("Found kronecker_delta(-1,0) != 0, kronecker_delta is not working.");
  if (kronecker_delta(-1, -1) == 1)
    PASSMSG("Found kronecker_delta(-1,-1) == 1, kronecker_delta is working.");
  else
    FAILMSG("Found kronecker_delta(-1,-1) != 1, kronecker_delta is not working.");

  unsigned uZero(0);
  unsigned uOne(1);
  long lZero(0);
  long lOne(1);

  if (kronecker_delta(uOne, uZero) == uZero)
    PASSMSG("Found kronecker_delta<unsigned>(uOne,uZero) == uZero, "
            "kronecker_delta is working.");
  else
    FAILMSG("Found kronecker_delta<unsigned>(uOne,uZero) != uZero, "
            "kronecker_delta is not working.");

  if (kronecker_delta(lOne, lZero) == static_cast<unsigned>(lZero))
    PASSMSG("Found kronecker_delta<long>(uOne,uZero) == uZero, "
            "kronecker_delta is working.");
  else
    FAILMSG("Found kronecker_delta<long>(uOne,uZero) != uZero, "
            "kronecker_delta is not working.");

  float fZero(0.0);
  float fOne(1.0);
  double dZero(0.0);
  double dOne(1.0);

  if (kronecker_delta(fOne, fZero) == uZero)
    PASSMSG("Found kronecker_delta<float>(fOne,fZero) == uZero, "
            "kronecker_delta is working.");
  else
    FAILMSG("Found kronecker_delta<float>(fOne,fZero) != uZero, "
            "kronecker_delta is not working.");

  if (kronecker_delta(dOne, dZero) == uZero)
    PASSMSG("Found kronecker_delta<double>(dOne,dZero) == uZero, "
            "kronecker_delta is working.");
  else
    FAILMSG("Found kronecker_delta<double>(dOne,dZero) != uZero, "
            "kronecker_delta is not working.");

  return;
}

//------------------------------------------------------------------------------------------------//
void tstFactorial(rtt_dsxx::UnitTest &ut) {
  using rtt_sf::factorial;

  // Test factorial

  if (factorial(0) == 1)
    PASSMSG("Found factorial(0) == 1, factorial is working.");
  else
    FAILMSG("Found factorial(0) != 1, factorial is not working.");
  if (factorial(1) == 1)
    PASSMSG("Found factorial(1) == 1, factorial is working.");
  else
    FAILMSG("Found factorial(1) != 1, factorial is not working.");
  if (factorial(2) == 2)
    PASSMSG("Found factorial(2) == 2, factorial is working.");
  else
    FAILMSG("Found factorial(2) != 2, factorial is not working.");
  if (factorial(3) == 6)
    PASSMSG("Found factorial(3) == 6, factorial is working.");
  else
    FAILMSG("Found factorial(3) != 6, factorial is not working.");

  try {
    factorial(13);
    FAILMSG("factorial(13) failed to trigger out of range error");
  } catch (std::range_error &) {
    PASSMSG("factorial(13) correctly triggered out of range error");
  }

  if (factorial(-3) == 1)
    PASSMSG("Found factorial(-3) == 1, factorial is working.");
  else
    FAILMSG("Found factorial(-3) != 1, factorial is not working.");

  unsigned uOne(1);
  long lOne(1);

  if (factorial(uOne) == uOne)
    PASSMSG("Found factorial<unsigned>(1) == unsigned(1), factorial is working.");
  else
    FAILMSG("Found factorial<unsigned>(1) != unsigned(1), factorial is not "
            "working.");
  if (factorial(lOne) == lOne)
    PASSMSG("Found factorial<long>(1) == long(1), factorial is working.");
  else
    FAILMSG("Found factorial<long>(1) != long(1), factorial is not working.");

  return;
}

//------------------------------------------------------------------------------------------------//
void tstFF(rtt_dsxx::UnitTest &ut) {
  using rtt_sf::factorial_fraction;

  // Test factorial_fraction

  if (rtt_dsxx::soft_equiv(factorial_fraction(1, 1), 1.0))
    PASSMSG("factorial_fraction(1,1) == 1");
  else
    FAILMSG("factorial_fraction(1,1) != 1");

  if (rtt_dsxx::soft_equiv(factorial_fraction(6, 4), 30.0))
    PASSMSG("factorial_fraction(6,4) == 30");
  else
    FAILMSG("factorial_fraction(6,4) != 30");

  if (rtt_dsxx::soft_equiv(factorial_fraction(1, 2), 0.5))
    PASSMSG("factorial_fraction(1,2) == 0.5");
  else
    FAILMSG("factorial_fraction(1,2) != 0.5");

  if (rtt_dsxx::soft_equiv(factorial_fraction(100, 99), 100.0))
    PASSMSG("factorial_fraction(100,99) == 100");
  else
    FAILMSG("factorial_fraction(100,99) != 100");

  long lOne(1);
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(factorial_fraction(lOne, lOne), 1.0));

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  using namespace rtt_sf;
  using namespace std;
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tstFactorial(ut);
    tstKdelta(ut);
    tstFF(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of test_sf.cc
//------------------------------------------------------------------------------------------------//
