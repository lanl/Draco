//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/test/tstExpInt.cc
 * \author Paul Talbot
 * \date   Thu Jul 28 09:20:34 2011
 * \brief  Tests the ExpInt for correct solutions in each routine
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "special_functions/ExpInt.hh"
#include <sstream>

using rtt_dsxx::soft_equiv;
using namespace rtt_sf;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

//
// Ei(x) Tests
//

//test Ei(x) with x < -log(eps) (small x routine)
void tstEi_low(rtt_dsxx::UnitTest &ut) {
  double x = 0.1;
  double val = Ei(x);
  double expVal = -1.6228128298555;

  std::ostringstream msg;
  if (soft_equiv(val, expVal)) {
    msg << "Ei(" << x << ") returned the expected value, " << expVal << ".";
    ut.passes(msg.str());
  } else {
    msg << "Ei(" << x << ") did NOT return the expected value.\n"
        << "\tExpected " << expVal << ", but got " << val;
    ut.failure(msg.str());
  }
  return;
}

//test Ei(x) with x > -log(eps) (small x routine)
void tstEi_high(rtt_dsxx::UnitTest &ut) {
  double x = 40;
  double val = Ei(x);
  double const ten(10.0);
  double expVal = 6.0397182636112 * std::pow(ten, 15);

  std::ostringstream msg;
  if (soft_equiv(val, expVal)) {
    msg << "Ei(" << x << ") returned the expected value, " << expVal << ".";
    ut.passes(msg.str());
  } else {
    msg << "Ei(" << x << ") did NOT return the expected value.\n"
        << "\tExpected " << expVal << ", but got " << val;
    ut.failure(msg.str());
  }
  return;
}

//test Ei(x) with x < 0 (negative x routine, -E_1(x)
void tstEi_neg(rtt_dsxx::UnitTest &ut) {
  double x = -5.2;
  double val = Ei(x);
  double expVal = -0.00090862161244866;

  std::ostringstream msg;
  if (soft_equiv(val, expVal)) {
    msg << "Ei(" << x << ") returned the expected value, " << expVal << ".";
    ut.passes(msg.str());
  } else {
    msg << "Ei(" << x << ") did NOT return the expected value.\n"
        << "\tExpected " << expVal << ", but got " << val;
    ut.failure(msg.str());
  }
  return;
}

//test E_0(x) (special case)
void tstE0(rtt_dsxx::UnitTest &ut) {
  unsigned n = 0;
  double x = 3.14;
  double val = En(n, x);
  double expVal = 0.013784330542027;

  std::ostringstream msg;
  if (soft_equiv(val, expVal)) {
    msg << "E_" << n << "(" << x << ") returned the expected value, " << expVal
        << ".";
    ut.passes(msg.str());
  } else {
    msg << "E_" << n << "(" << x << ") did NOT return the expected value.\n"
        << "\tExpected " << expVal << ", but got " << val;
    ut.failure(msg.str());
  }
  return;
}

//test E_3(0) (special case)
void tstE3_0(rtt_dsxx::UnitTest &ut) {
  unsigned n = 3;
  double x = 0;
  double val = En(n, x);
  double expVal = 0.5;

  std::ostringstream msg;
  if (soft_equiv(val, expVal)) {
    msg << "E_" << n << "(" << x << ") returned the expected value, " << expVal
        << ".";
    ut.passes(msg.str());
  } else {
    msg << "E_" << n << "(" << x << ") did NOT return the expected value.\n"
        << "\tExpected " << expVal << ", but got " << val;
    ut.failure(msg.str());
  }
  return;
}

//test E_4(x), x < 1
void tstE4_low(rtt_dsxx::UnitTest &ut) {
  unsigned n = 4;
  double x = 0.3;
  double val = En(n, x);
  double expVal = 0.2169352242375;

  std::ostringstream msg;
  if (soft_equiv(val, expVal)) {
    msg << "E_" << n << "(" << x << ") returned the expected value, " << expVal
        << ".";
    ut.passes(msg.str());
  } else {
    msg << "E_" << n << "(" << x << ") did NOT return the expected value.\n"
        << "\tExpected " << expVal << ", but got " << val;
    ut.failure(msg.str());
  }
  return;
}

//test E_1(x), x > 1
void tstE1_high(rtt_dsxx::UnitTest &ut) {
  unsigned n = 1;
  double x = 5;
  double val = En(n, x);
  double expVal = 0.0011482955912753;

  std::ostringstream msg;
  if (soft_equiv(val, expVal)) {
    msg << "E_" << n << "(" << x << ") returned the expected value, " << expVal
        << ".";
    ut.passes(msg.str());
  } else {
    msg << "E_" << n << "(" << x << ") did NOT return the expected value.\n"
        << "\tExpected " << expVal << ", but got " << val;
    ut.failure(msg.str());
  }
  return;
}

//------------------------------------------------------------------------------------------------//
// RUN TESTS
//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tstEi_low(ut);
    tstEi_high(ut);
    tstEi_neg(ut);
    tstE0(ut);
    tstE3_0(ut);
    tstE4_low(ut);
    tstE1_high(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstExpInt.cc
//------------------------------------------------------------------------------------------------//
