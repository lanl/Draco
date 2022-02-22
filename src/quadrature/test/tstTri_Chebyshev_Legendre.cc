//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/test/tstTri_Chebyshev_Legendre.cc
 * \author Kent G. Budge
 * \date   Tue Nov  6 13:08:49 2012
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "quadrature_test.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "quadrature/Tri_Chebyshev_Legendre.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_quadrature;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    {
      Tri_Chebyshev_Legendre quadrature(8); // SN order = 8
      if (quadrature.sn_order() != 8) {
        ut.failure("NOT correct SN order");
      }
      quadrature_test(ut, quadrature);
    }
    {
      Tri_Chebyshev_Legendre quadrature(8, 1, 2); // SN order = 8, mu=1, eta=2
      quadrature_test(ut, quadrature);
    }
    Tri_Chebyshev_Legendre quadrature4(4);
    quadrature_integration_test(ut, quadrature4);
    Tri_Chebyshev_Legendre quadrature8(8);
    quadrature_integration_test(ut, quadrature8);
    Tri_Chebyshev_Legendre quadrature12(12);
    quadrature_integration_test(ut, quadrature12);
    Tri_Chebyshev_Legendre quadrature16(16);
    quadrature_integration_test(ut, quadrature16);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstTri_Chebyshev_Legendre.cc
//------------------------------------------------------------------------------------------------//
