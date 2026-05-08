//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/test/tstbtridag.cc
 * \author Kent Budge
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "linear/btridag.hh"
#include <vector>

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    PASSMSG("Starting the test.");
    std::vector<double> a(2), b = {2., 5.}, c(2), r = {1., 2.}, u(2);
    c[0] = 3.;
    a[1] = 1.;
    rtt_linear::btridag(a, b, c, r, 2U, 1U, u);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(1.0, u[0] * 2 + u[1] * 3));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(2.0, u[0] * 1 + u[1] * 5));
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstbtridag.cc
//------------------------------------------------------------------------------------------------//
