//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/test/tstLobatto.cc
 * \author Kent G. Budge
 * \date   Tue Nov  6 13:08:49 2012
 * \note   Copyright (C) 2012-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"
#include "quadrature_test.hh"
#include "quadrature/Lobatto.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_quadrature;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    for (unsigned n = 4; n < 16; n += 2) {
      Lobatto quadrature(n);

      quadrature_test(ut, quadrature);
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstLobatto.cc
//------------------------------------------------------------------------------------------------//
