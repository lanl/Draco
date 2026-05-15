//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/test/tstProduct_Chebyshev_Legendre.cc
 * \author James S. Warsa
 * \note   Copyright (C) 2013-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"
#include "quadrature_test.hh"
#include "quadrature/Product_Chebyshev_Legendre.hh"

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
      // SN order = 8, Azimuthal order = 12

      Product_Chebyshev_Legendre quadrature(8, 12);
      if (quadrature.sn_order() != 8) {
        FAILMSG("NOT correct SN order");
      }
      quadrature_test(ut, quadrature);
    }
    {
      // SN order = 8, Azimuthal order = 12

      Product_Chebyshev_Legendre quadrature(8, 12, 1, 2);
      quadrature_test(ut, quadrature);
    }

    Product_Chebyshev_Legendre quadrature4(4, 8);
    quadrature_integration_test(ut, quadrature4);
    Product_Chebyshev_Legendre quadrature8(8, 16);
    quadrature_integration_test(ut, quadrature8);
    Product_Chebyshev_Legendre quadrature10(10, 20);
    quadrature_integration_test(ut, quadrature10);
    Product_Chebyshev_Legendre quadrature12(12, 24);
    quadrature_integration_test(ut, quadrature12);
    Product_Chebyshev_Legendre quadrature16(16, 32);
    quadrature_integration_test(ut, quadrature16);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstProduct_Chebyshev_Legendre.cc
//------------------------------------------------------------------------------------------------//
