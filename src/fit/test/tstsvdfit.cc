//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   fit/test/tstsvdfit.cc
 * \author Kent Budge
 * \date   Tue Aug 26 12:02:36 2008
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "fit/svdfit.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_fit;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void funcs(double const x, vector<double> &y) {
  y.resize(3);
  y[0] = 1;
  y[1] = x;
  y[2] = x * x;
}

//------------------------------------------------------------------------------------------------//
void tstsvdfit(UnitTest &ut) {
  unsigned const N = 10;
  vector<double> x(N), y(N), sig(N);
  vector<double> a(3), u, v, w;
  double chisq;

  for (unsigned i = 0; i < 10; i++) {
    x[i] = i;
    y[i] = 3.2 + i * (1.7 + 2.1 * i);
    sig[i] = 1;
  }

  svdfit(x, y, sig, a, u, v, w, chisq, funcs, 1.0e-12);

  if (chisq < 1.0e-25 && soft_equiv(a[0], 3.2) && soft_equiv(a[1], 1.7) && soft_equiv(a[2], 2.1)) {
    ut.passes("fit is good");
  } else {
    ut.failure("fit is NOT good");
  }

  for (unsigned i = 0; i < 10; ++i) {
    y[i] += (i % 2 ? -1.0e-5 : 1.0e-5);
  }

  svdfit(x, y, sig, a, u, v, w, chisq, funcs, 1.0e-12);

  if (chisq < 1.0e-8 && soft_equiv(a[0], 3.2, 1.0e-5) && soft_equiv(a[1], 1.7, 1.0e-5) &&
      soft_equiv(a[2], 2.1, 1.0e-5)) {
    ut.passes("fit is good");
  } else {
    ut.failure("fit is NOT good");
  }
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstsvdfit(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstsvdfit.cc
//------------------------------------------------------------------------------------------------//
