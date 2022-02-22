//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/test/tstTri_Chebyshev_Legendre.cc
 * \author Kent G. Budge
 * \date   Tue Nov  6 13:08:49 2012
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//
// $Id: template_test.cc 5830 2011-05-05 19:43:43Z kellyt $
//------------------------------------------------------------------------------------------------//

#include "quadrature_test.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "quadrature/Square_Chebyshev_Legendre.hh"

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
      Square_Chebyshev_Legendre quadrature(8); // SN order = 8
      if (quadrature.sn_order() != 8) {
        ut.failure("NOT correct SN order");
      }
      quadrature_test(ut, quadrature);
    }
    {
      Square_Chebyshev_Legendre quadrature(8, 1, 2);
      quadrature_test(ut, quadrature);
    }

    Square_Chebyshev_Legendre quadrature4(4);
    quadrature_integration_test(ut, quadrature4);
    Square_Chebyshev_Legendre quadrature8(8);
    quadrature_integration_test(ut, quadrature8);
    Square_Chebyshev_Legendre quadrature10(10);
    quadrature_integration_test(ut, quadrature10);
    Square_Chebyshev_Legendre quadrature12(12);
    quadrature_integration_test(ut, quadrature12);
    Square_Chebyshev_Legendre quadrature16(16);
    quadrature_integration_test(ut, quadrature16);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstSquare_Chebyshev_Legendre.cc
//------------------------------------------------------------------------------------------------//
