//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/test/tstGeneral_Octant_Quadrature.cc
 * \author Kent G. Budge
 * \date   Tue Nov  6 13:08:49 2012
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "quadrature_test.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "quadrature/General_Octant_Quadrature.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_quadrature;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    double const V = 1 / sqrt(3.0);
    vector<double> mu(1, V), eta(1, V), xi(1, V), wt(1, 1.0);
    General_Octant_Quadrature quadrature(2, mu, eta, xi, wt, 2, TRIANGLE_QUADRATURE);

    quadrature_test(ut, quadrature);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstGeneral_Octant_Quadrature.cc
//------------------------------------------------------------------------------------------------//
