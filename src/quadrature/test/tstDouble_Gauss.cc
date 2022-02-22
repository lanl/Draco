//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/test/tstDouble_Gauss.cc
 * \author Kent G. Budge
 * \date   Tue Nov  6 13:08:49 2012
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "quadrature_test.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "quadrature/Double_Gauss.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_quadrature;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    Double_Gauss quad8(8); // SN order = 8
    quadrature_test(ut, quad8);

    Double_Gauss quad2(2); // SN order = 2
    quadrature_test(ut, quad2);

    // Test
    {
      FAIL_IF_NOT(quadrature_interpolation_model_as_text(SN) == "SN");
      FAIL_IF_NOT(quadrature_interpolation_model_as_text(GQ1) == "GQ1");
      FAIL_IF_NOT(quadrature_interpolation_model_as_text(GQ2) == "GQ2");
      FAIL_IF_NOT(quadrature_interpolation_model_as_text(GQF) == "GQF");
      try {
        FAIL_IF_NOT(quadrature_interpolation_model_as_text(SVD) == "SVD");
      } catch (rtt_dsxx::assertion const & /*error*/) {
        PASSMSG("assertion caught for unlisted QIM == SVD.");
      }
      try {
        FAIL_IF_NOT(quadrature_interpolation_model_as_text(END_QIM) == "END_QIM");
      } catch (rtt_dsxx::assertion const & /*error*/) {
        PASSMSG("assertion caught for invalid QIM == END_QIM.");
      }
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstDouble_Gauss.cc
//------------------------------------------------------------------------------------------------//
