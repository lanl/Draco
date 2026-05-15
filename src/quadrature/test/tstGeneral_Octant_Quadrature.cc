//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/test/tstGeneral_Octant_Quadrature.cc
 * \author Kent G. Budge
 * \note   Copyright (C) 2012-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"
#include "quadrature_test.hh"
#include "quadrature/General_Octant_Quadrature.hh"
#include "quadrature/Quadrature_Class.hh"
#include <numbers>
#include <vector>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_quadrature;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    double const V = std::numbers::inv_sqrt3;
    vector<double> const mu(1, V), eta(1, V), xi(1, V), wt(1, 1.0);
    General_Octant_Quadrature quadrature(2, mu, eta, xi, wt, 2, TRIANGLE_QUADRATURE);

    quadrature_test(ut, quadrature);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstGeneral_Octant_Quadrature.cc
//------------------------------------------------------------------------------------------------//
