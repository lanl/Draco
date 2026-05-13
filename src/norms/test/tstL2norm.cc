//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   norms/test/tstL2norm.cc
 * \author Kent Budge
 * \brief  Test the L2norm function template.
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "c4/C4_Functions.hh"
#include "c4/ParallelUnitTest.hh"
#include "norms/L2norm.hh"
#include <cmath>
#include <vector>

namespace {

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_norms;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstL2norm(UnitTest &ut) {
  unsigned const N = 2;
  unsigned const n = rtt_c4::nodes();

  vector<double> const x(N, rtt_c4::node() + 1);

  double norm = L2norm(x);
  ut.check(soft_equiv(norm, sqrt(1. / 6 + n * (0.5 + n / 3.))), "L2norm is correct");

  norm = L2norm_diff(x, x);
  ut.check(soft_equiv(norm, 0.0), "L2norm_diff of vector with self is zero");

  vector<double> const y(N, 2 * rtt_c4::node() + 2);
  norm = L2norm_diff(x, y);
  ut.check(soft_equiv(norm, sqrt(1. / 6 + n * (0.5 + n / 3.))), "L2norm_diff is correct");

  return;
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, release);
  try {
    tstL2norm(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstL2norm.cc
//------------------------------------------------------------------------------------------------//
