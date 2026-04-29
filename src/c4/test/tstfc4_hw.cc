//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstfc4_hw.cc
 * \author Kelly Thompson
 * \brief  Test C++ main linking a Fortran library that uses MPI
 * \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"
#include <cstddef> // size_t

//------------------------------------------------------------------------------------------------//
// forward declaration of f90 functions
extern "C" void tst_mpi_hw(size_t *numFail);

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    size_t nf(ut.numFails);
    // Call Fortran subroutine
    tst_mpi_hw(&nf);
    ut.numPasses = 1;
    ut.numFails = static_cast<unsigned>(nf);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of c4/test/tstfc4_hw.cc
//------------------------------------------------------------------------------------------------//
