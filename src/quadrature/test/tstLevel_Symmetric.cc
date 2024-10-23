//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/test/tstLevel_Symmetric.cc
 * \author Kent G. Budge
 * \date   Tue Nov  6 13:08:49 2012
 * \note   Copyright (C) 2012-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"
#include "quadrature_test.hh"
#include "quadrature/Level_Symmetric.hh"
#include <iostream>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_quadrature;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    for (unsigned N = 2; N < 25; N += 2)
    // Pathological in curvilinear geometry at N>=18
    {
      cout << "Order " << N << ":\n";
      Level_Symmetric quadrature(N);
      FAIL_IF_NOT(quadrature.sn_order() == N);
      bool const cartesian_tests_only(N >= 18);
      quadrature_test(ut, quadrature, cartesian_tests_only);
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstLevel_Symmetric.cc
//------------------------------------------------------------------------------------------------//
