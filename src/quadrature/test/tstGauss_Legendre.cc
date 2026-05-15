//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/test/tstGauss_Legendre.cc
 * \author Kent G. Budge
 * \note   Copyright (C) 2012-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"
#include "quadrature_test.hh"
#include "quadrature/Gauss_Legendre.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_quadrature;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    Gauss_Legendre quadrature(8); // SN order = 8
    quadrature_test(ut, quadrature);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstGauss_Legendre.cc
//------------------------------------------------------------------------------------------------//
