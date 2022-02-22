//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/test/tstqrupdt.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:39:20 2004
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "linear/qr_unpack.hh"
#include "linear/qrdcmp.hh"
#include "linear/qrupdt.hh"
#include <array>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void tstqrupdt(UnitTest &ut) {
  vector<double> A(2 * 2);
  A[0 + 2 * 0] = 2.;
  A[0 + 2 * 1] = 3.;
  A[1 + 2 * 0] = 1.;
  A[1 + 2 * 1] = 5.;

  vector<double> C, D;

  qrdcmp(A, 2, C, D);

  // Unpack
  vector<double> QT;
  qr_unpack(A, 2, C, D, QT);

  // Now do a Householder update.
  vector<double> s(2), t(2);
  vector<double> u(2), v(2);

  s[0] = -0.1;
  s[1] = 0.1;
  t[0] = 0.1;
  t[1] = 0.2;

  u[0] = QT[0 + 2 * 0] * s[0] + QT[0 + 2 * 1] * s[1];
  u[1] = QT[1 + 2 * 0] * s[0] + QT[1 + 2 * 1] * s[1];
  v = t;

  qrupdt(A, QT, 2, u, v);

  // Check the update

  std::array<double, 2 * 2> QR;
  QR[0 + 2 * 0] = QT[0 + 2 * 0] * A[0 + 2 * 0] + QT[0 + 2 * 1] * A[1 + 2 * 0];
  QR[0 + 2 * 1] = QT[0 + 2 * 0] * A[0 + 2 * 1] + QT[0 + 2 * 1] * A[1 + 2 * 1];
  QR[1 + 2 * 0] = QT[1 + 2 * 0] * A[0 + 2 * 0] + QT[1 + 2 * 1] * A[1 + 2 * 0];
  QR[1 + 2 * 1] = QT[1 + 2 * 0] * A[0 + 2 * 1] + QT[1 + 2 * 1] * A[1 + 2 * 1];

  if (soft_equiv(QR[0 + 2 * 0], 2.0 + s[0] * t[0])) {
    ut.passes("0,0 is correct");
  } else {
    ut.failure("0,0 is NOT correct");
  }
  if (soft_equiv(QR[0 + 2 * 1], 3.0 + s[0] * t[1])) {
    ut.passes("0,1 is correct");
  } else {
    ut.failure("0,1 is NOT correct");
  }
  if (soft_equiv(QR[1 + 2 * 0], 1.0 + s[1] * t[0])) {
    ut.passes("1,0 is correct");
  } else {
    ut.failure("1,0 is NOT correct");
  }
  if (soft_equiv(QR[1 + 2 * 1], 5.0 + s[1] * t[1])) {
    ut.passes("1,1 is correct");
  } else {
    ut.failure("1,1 is NOT correct");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstqrupdt(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstqrupdt.cc
//------------------------------------------------------------------------------------------------//
