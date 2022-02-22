//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstfc4_hw.cc
 * \author Kelly Thompson
 * \date   Thursday, Nov 12, 2015, 10:35 am
 * \brief  Test C++ main linking a Fortran library that uses MPI
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include <iostream>

//------------------------------------------------------------------------------------------------//
// forward declaration of f90 functions
extern "C" void tst_mpi_hw(size_t *numFail);

//------------------------------------------------------------------------------------------------//
void test_mpi_hw(rtt_dsxx::UnitTest &ut) {
  // size_t np(ut.numPasses);
  size_t nf(ut.numFails);
  // Call Fortran subroutine
  tst_mpi_hw(&nf);
  ut.numPasses = 1;
  ut.numFails = static_cast<unsigned>(nf);
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_mpi_hw(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of c4/test/tstfc4_hw.cc
//------------------------------------------------------------------------------------------------//
