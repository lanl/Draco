//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstInterpolate.cc
 * \author Kent G. Budge, Kendra P. Long
 * \date   Tues Nov 15 2022
 * \brief  Test functions defined in ds++/Interpolate.hh.
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/Interpolate.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/dbc.hh"

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tst1dlin(rtt_dsxx::UnitTest &ut) {
  const auto initfails = ut.numFails;
  // function y = 2.5 * x - 1.0

  // define boundary points
  double x1{1.0};
  double y1{2.5 * x1 - 1.0};
  double x2{3.0};
  double y2{2.5 * x2 - 1.0};

  double x = 1.452;
  double y = rtt_dsxx::interpolate::linear_1d(x1, x2, y1, y2, x);
  double ref = 2.5 * x - 1.0;

  FAIL_IF_NOT(rtt_dsxx::soft_equiv(y, ref));

  // try another one
  x1 = 1.45;
  y1 = 2.5 * x1 - 1.0;
  x2 = 1.1;
  y2 = 2.5 * x2 - 1.0;

  x = 1.33;
  y = rtt_dsxx::interpolate::linear_1d(x1, x2, y1, y2, x);
  ref = 2.5 * x - 1.0;

  FAIL_IF_NOT(rtt_dsxx::soft_equiv(y, ref));

  if ((ut.numFails - initfails) == 0)
    PASSMSG("1D linear interpolation test passes.");
  else
    FAILMSG("1D linear interpolation test fails.");
  return;
}

void tst3dlin(rtt_dsxx::UnitTest &ut) {
  const auto initfails = ut.numFails;
  // define boundary points
  double x0{1.0};
  double x1{2.0};
  double y0{3.0};
  double y1{4.0};
  double z0{5.0};
  double z1{6.0};
  double f000{1.0};
  double f100{1.0};
  double f001{3.0};
  double f101{3.0};
  double f010{5.0};
  double f110{5.0};
  double f011{7.0};
  double f111{7.0};

  double x{1.75};
  double y{3.25};
  double z{5.5};

  double f = rtt_dsxx::interpolate::linear_3d(x0, x1, y0, y1, z0, z1, f000, f100, f001, f101, f010,
                                              f110, f011, f111, x, y, z);
  double ref = 3.0;
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(f, ref));

  // try another one
  x = 1.5;
  y = 3.99;
  z = 5.01;

  f = rtt_dsxx::interpolate::linear_3d(x0, x1, y0, y1, z0, z1, f000, f100, f001, f101, f010, f110,
                                       f011, f111, x, y, z);
  ref = 4.98;
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(f, ref));

  if ((ut.numFails - initfails) == 0)
    PASSMSG("3D linear interpolation test passes.");
  else
    FAILMSG("3D linear interpolation test fails.");
  return;
}

//------------------------------------------------------------------------------------------------//
void tst1dlag(rtt_dsxx::UnitTest &ut) {
  const auto initfails = ut.numFails;

  const size_t n_break = 2;
  const size_t n_local = 3;
  const std::vector<double> points{0.0, 0.5, 1.0};
  const std::vector<double> yvals{0.0, 2.5, 5.0};
  std::vector<double> multipliers =
      rtt_dsxx::interpolate::lagrange_multipliers(n_break, n_local, points);
  double x{0.75};
  double fx = rtt_dsxx::interpolate::lagrange_1d(yvals, points, multipliers, x);
  double ref{5.0 * x};
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(fx, ref));

  x = 0.33;
  fx = rtt_dsxx::interpolate::lagrange_1d(yvals, points, multipliers, x);
  ref = 5.0 * x;
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(fx, ref));

  if ((ut.numFails - initfails) == 0)
    PASSMSG("1D lagrange interpolation test passes.");
  else
    FAILMSG("1D lagrange interpolation test fails.");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tst1dlin(ut);
    tst3dlin(ut);
    tst1dlag(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstInterpolate.cc
//------------------------------------------------------------------------------------------------//
