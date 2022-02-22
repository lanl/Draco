//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/test/tstLevel_Symmetric.cc
 * \author Kent G. Budge
 * \date   Tue Nov  6 13:08:49 2012
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "quadrature_test.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "quadrature/Level_Symmetric.hh"

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
      cout << "Order " << N << ':' << endl;
      Level_Symmetric quadrature(N);
      if (quadrature.sn_order() != N) {
        ut.failure("NOT correct SN order");
      }
      bool cartesian_tests_only(N >= 18);
      quadrature_test(ut, quadrature, cartesian_tests_only);
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstLevel_Symmetric.cc
//------------------------------------------------------------------------------------------------//
