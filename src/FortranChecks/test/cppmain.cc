//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   FortranChecks/test/cppmain.cc
 * \author Kelly Thompson
 * \date   Tuesday, Jun 12, 2012, 16:03 pm
 * \brief  Test C++ main linking a Fortran library.
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include <iostream>

//------------------------------------------------------------------------------------------------//
// forward declaration of f90 functions
extern "C" void sub1(double alpha, size_t *numPass, size_t *numFail);

//------------------------------------------------------------------------------------------------//
void test_isocbinding_sub1(rtt_dsxx::UnitTest &ut) {
  double alpha = 1.0;
  size_t np(ut.numPasses);
  size_t nf(ut.numFails);
  // Call fortran subroutine
  sub1(alpha, &np, &nf);
  ut.numPasses = static_cast<unsigned>(np);
  ut.numFails = static_cast<unsigned>(nf);
  std::cout << ut.numPasses << " " << ut.numFails << std::endl;
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_isocbinding_sub1(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of cppmain.cc
//------------------------------------------------------------------------------------------------//
