//--------------------------------------------*-C++-*---------------------------------------------//
/*!
* \file   linear/test/tstfnorm.cc
* \author Kent Budge
* \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "linear/fnorm.hh"
#include <cmath> // cos, sin
#include <vector>

namespace {

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void func(const std::vector<double> &x, std::vector<double> &fvec) {
  fvec.resize(2);
  fvec[0] = x[1] * sin(x[0]);
  fvec[1] = x[1] * cos(x[0]);
}

//------------------------------------------------------------------------------------------------//
void tst_f_norm(rtt_dsxx::UnitTest &ut) {
  PASSMSG("Starting the test.");
  std::vector<double> const x = {0.235, 3.2};
  std::vector<double> fvec;
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(rtt_linear::fnorm(x, fvec, &func), 0.5 * 3.2 * 3.2));
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tst_f_norm(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstfnorm.cc
//------------------------------------------------------------------------------------------------//
