//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   norms/test/tstL2norm.cc
 * \author Kent Budge
 * \date   Tue Sep 18 09:06:26 2007
 * \brief  Test the L2norm function template.
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include "norms/L2norm.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_norms;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstL2norm(UnitTest &ut) {
  unsigned const N = 2;
  unsigned const n = rtt_c4::nodes();

  vector<double> x(N, rtt_c4::node() + 1);

  double norm = L2norm(x);
  ut.check(soft_equiv(norm, sqrt(1. / 6 + n * (0.5 + n / 3.))), "L2norm is correct");

  norm = L2norm_diff(x, x);
  ut.check(soft_equiv(norm, 0.0), "L2norm_diff of vector with self is zero");

  vector<double> y(N, 2 * rtt_c4::node() + 2);
  norm = L2norm_diff(x, y);
  ut.check(soft_equiv(norm, sqrt(1. / 6 + n * (0.5 + n / 3.))), "L2norm_diff is correct");

  return;
}

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
