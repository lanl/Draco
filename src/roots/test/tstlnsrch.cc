//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   roots/test/tstlnsrch.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:39:20 2004
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "roots/lnsrch.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;
using namespace rtt_roots;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void func(const vector<double> &x, vector<double> &fvec) {
  fvec.resize(2);
  fvec[0] = 7.2 * x[0] + 3.5 * x[1] + 2.3;
  fvec[1] = -2.2 * x[0] + 2.7 * x[1] + 5.4;
}

//------------------------------------------------------------------------------------------------//
void tstlnsrch(UnitTest &ut) {
  vector<double> x = {0.0, 0.0};
  vector<double> xold = x;
  vector<double> fvec;
  double fold = fnorm(x, fvec, &func);

  vector<double> g(2);
  g[0] = 7.2 * fvec[0] - 2.2 * fvec[1];
  g[1] = 3.5 * fvec[0] + 2.7 * fvec[1];

  vector<double> p = {-0.4675753, 1.6190125};

  double f;
  bool check;
  lnsrch(xold, fold, g, p, x, f, check, fvec, &func, 0.01, 0.0);

  if (!check) {
    FAILMSG("lnsrch did NOT bomb gracefully");
  } else {
    PASSMSG("lnsrch bombed gracefully");
  }

  for (double &pvalue : p)
    pvalue = -pvalue;
  lnsrch(xold, fold, g, p, x, f, check, fvec, &func, 0.01, 0.0);

  if (check || f > 1.0e-10) {
    FAILMSG("lnsrch did NOT succeed");
  } else {
    PASSMSG("lnsrch successful");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  try {
    ScalarUnitTest ut(argc, argv, release);
    tstlnsrch(ut);
  } catch (exception &err) {
    cout << "ERROR: While testing tstlnsrch, " << err.what() << endl;
    return 1;
  } catch (...) {
    cout << "ERROR: While testing tstlnsrch, "
         << "An unknown exception was thrown." << endl;
    return 1;
  }
  return 0;
}

//------------------------------------------------------------------------------------------------//
// end of tstlnsrch.cc
//------------------------------------------------------------------------------------------------//
