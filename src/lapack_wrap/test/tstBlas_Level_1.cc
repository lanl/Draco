//----------------------------------*-C++-*----------------------------------//
/*!
* \file   lapack_wrap/test/tstBlas_Level_1.cc
* \brief  Test Blas level 1 wrap.
* \note   Copyright (C) 2016 Los Alamos National Security, LLC
*         All rights reserved.
* \version $Id$
*/
//---------------------------------------------------------------------------//

#include "lapack_wrap/Blas.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"

using namespace std;

using rtt_lapack_wrap::blas_copy;
using rtt_lapack_wrap::blas_scal;
using rtt_lapack_wrap::blas_dot;
using rtt_lapack_wrap::blas_axpy;
using rtt_lapack_wrap::blas_nrm2;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
template <typename T> void tst_copy(rtt_dsxx::UnitTest &ut) {
  vector<T> x(10, 0.0);
  vector<T> y(10, 0.0);

  for (int i = 0; i < 10; i++)
    x[i] = 1.2 + i;

  blas_copy(10, &x[0], 1, &y[0], 1);
  if (!soft_equiv(x.begin(), x.end(), y.begin(), y.end()))
    ITFAILS;

  std::fill(y.begin(), y.end(), static_cast<T>(0.0));
  if (soft_equiv(x.begin(), x.end(), y.begin(), y.end()))
    ITFAILS;

  blas_copy(x, 1, y, 1);
  if (!soft_equiv(x.begin(), x.end(), y.begin(), y.end()))
    ITFAILS;

  string ttype(typeid(T).name());
  if (ut.numFails == 0)
    PASSMSG(std::string("BLAS copy (") + ttype +
            std::string("float) tests ok."));
  else
    PASSMSG(std::string("BLAS copy (") + ttype + std::string("float) fails."));
  return;
}

//---------------------------------------------------------------------------//
template <typename T> void tst_scal(rtt_dsxx::UnitTest &ut) {
  vector<double> x(10, 0.0);
  vector<double> y(10, 0.0);
  vector<double> ref(10, 0.0);

  double alpha = 10.0;

  for (int i = 0; i < 10; i++) {
    y[i] = 1.2 + i;
    ref[i] = alpha * y[i];
  }

  x = y;
  blas_scal(10, alpha, &x[0], 1);
  if (!soft_equiv(x.begin(), x.end(), ref.begin(), ref.end()))
    ITFAILS;

  x = y;
  blas_scal(alpha, x, 1);
  if (!soft_equiv(x.begin(), x.end(), ref.begin(), ref.end()))
    ITFAILS;

  string ttype(typeid(T).name());
  if (ut.numPasses > 0 && ut.numFails == 0)
    PASSMSG(std::string("BLAS scal (") + ttype +
            std::string("float) tests ok."));
  else
    PASSMSG(std::string("BLAS scal (") + ttype + std::string("float) fails."));
  return;
}

//---------------------------------------------------------------------------//
template <typename T> void tst_dot(rtt_dsxx::UnitTest &ut) {
  vector<double> x(10, 0.0);
  vector<double> y(10, 0.0);

  double ref = 0.0;
  double dot = 0.0;

  for (int i = 0; i < 10; i++) {
    x[i] = i + 1.5;
    y[i] = (x[i] + 2.5) / 2.1;
    ref += x[i] * y[i];
  }

  dot = blas_dot(10, &x[0], 1, &y[0], 1);
  if (!soft_equiv(dot, ref))
    ITFAILS;
  dot = 0.0;

  dot = blas_dot(x, 1, y, 1);
  if (!soft_equiv(dot, ref))
    ITFAILS;
  dot = 0.0;

  string ttype(typeid(T).name());
  if (ut.numPasses > 0 && ut.numFails == 0)
    PASSMSG(std::string("BLAS dot (") + ttype +
            std::string("float) tests ok."));
  else
    PASSMSG(std::string("BLAS dot (") + ttype + std::string("float) fails."));
  return;
}

//---------------------------------------------------------------------------//
template <typename T> void tst_axpy(rtt_dsxx::UnitTest &ut) {
  vector<double> x(10, 0.0);
  vector<double> y(10, 1.0);
  vector<double> ref(10, 0.0);

  double alpha = 10.0;

  for (size_t i = 0; i < x.size(); i++) {
    x[i] = static_cast<double>(i) + 1.0;
    ref[i] = alpha * x[i] + y[i];
  }

  blas_axpy(10, alpha, &x[0], 1, &y[0], 1);
  if (!soft_equiv(y.begin(), y.end(), ref.begin(), ref.end()))
    ITFAILS;

  fill(y.begin(), y.end(), 1.0);
  blas_axpy(alpha, x, 1, y, 1);
  if (!soft_equiv(y.begin(), y.end(), ref.begin(), ref.end()))
    ITFAILS;

  string ttype(typeid(T).name());
  if (ut.numPasses > 0 && ut.numFails == 0)
    PASSMSG(std::string("BLAS axpy (") + ttype +
            std::string("float) tests ok."));
  else
    PASSMSG(std::string("BLAS axpy (") + ttype + std::string("float) fails."));
  return;
}

//---------------------------------------------------------------------------//
template <typename T> void tst_nrm2(rtt_dsxx::UnitTest &ut) {
  vector<double> x(10, 0.0);

  double ref = 0.0;
  double nrm = 0.0;

  for (size_t i = 0; i < x.size(); i++) {
    x[i] = 1.25 + (1.0 - i * 0.5);
    ref += x[i] * x[i];
  }
  ref = sqrt(ref);

  nrm = blas_nrm2(10, &x[0], 1);
  if (!soft_equiv(nrm, ref))
    ITFAILS;
  nrm = 0.0;

  nrm = blas_nrm2(x.begin(), x.end());
  if (!soft_equiv(nrm, ref))
    ITFAILS;
  nrm = 0.0;

  nrm = blas_nrm2(x, 1);
  if (!soft_equiv(nrm, ref))
    ITFAILS;

  string ttype(typeid(T).name());
  if (ut.numPasses > 0 && ut.numFails == 0)
    PASSMSG(std::string("BLAS nrm2 (") + ttype +
            std::string("float) tests ok."));
  else
    PASSMSG(std::string("BLAS nrm2 (") + ttype + std::string("float) fails."));
  return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tst_copy<float>(ut);
    tst_scal<float>(ut);
    tst_dot<float>(ut);
    tst_axpy<float>(ut);
    tst_nrm2<float>(ut);
    tst_copy<double>(ut);
    tst_scal<double>(ut);
    tst_dot<double>(ut);
    tst_axpy<double>(ut);
    tst_nrm2<double>(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstBlas_Level_1.cc
//---------------------------------------------------------------------------//
