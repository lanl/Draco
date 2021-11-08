//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   min/test/tstmrqmin.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:39:20 2004
 * \note   Copyright (C) 2010-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "min/mrqmin.hh"
#include <algorithm>
#include <fstream>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_min;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

using Model = void(vector<double> const &x, vector<double> const &a, double &y,
                   vector<double> &dyda);

void model(vector<double> const &x, vector<double> const &a, double &y, vector<double> &dyda) {
  Require(x.size() == 4);
  Require(a.size() == 9);
  dyda.resize(9);

  double const A = a[0];
  double const B = a[1];
  double const C = a[2];
  double const D = a[3];
  double const E = a[4];
  double const F = a[5];
  double const G = a[6];
  double const H = a[7];
  double const J = a[8];

  double const p = x[0];
  double const n = x[1];
  double const c = x[2];
  double const an = x[3];

  y = A + B * c + C * c * an + n * (D + E * c + F * c * an + p * (G + H * c + J * c * an));
  dyda[0] = 1;
  dyda[1] = c;
  dyda[2] = c * an;
  dyda[3] = n;
  dyda[4] = n * c;
  dyda[5] = n * c * an;
  dyda[6] = n * p;
  dyda[7] = n * p * c;
  dyda[8] = n * p * c * an;
}

//------------------------------------------------------------------------------------------------//
// tstmrqmin
//------------------------------------------------------------------------------------------------//

void tstmrqmin(UnitTest &ut) {
  std::cout << "\nTesting tstmrqmin()...\n" << std::endl;

  // Parse an input file
  std::string filename(ut.getTestSourcePath() + std::string("mrqmin_data.inp"));

  ifstream data(filename.c_str());
  vector<double> x, y, sig;
  for (;;) {
    double ymin = 1e100;
    for (unsigned i = 0; i < 3; ++i) {
      double s, n, c, a, t;
      data >> s >> n >> c >> a >> t;
      if (!data) {
        break;
      }
      if (i == 0) {
        x.push_back(s);
        x.push_back(n);
        x.push_back(c);
        x.push_back(a);
      }
      ymin = min(t, ymin);
    }
    if (!data) {
      break;
    }
    {
      y.push_back(ymin);
      sig.push_back(ymin);
    }
  }

  vector<double> a(9, 0.0);
  vector<bool> ia(9, false);
  double chisq{0.0};
  double alamda{0.0};
  vector<double> covar, alpha;
  unsigned iopt(0);
  double copt = 1e100;
  for (unsigned i = 0; i < 9; ++i) {
    ia[i] = true;
    for (unsigned j = 0; j < 3; ++j) {
      if (j == 0) {
        alamda = -1;
      }
      if (j == 2) {
        alamda = 0;
      }
      Check(y.size() < UINT_MAX);
      mrqmin(x, y, sig, static_cast<unsigned>(y.size()), 4U, a, ia, covar, alpha, 9U, chisq, model,
             alamda);
    }
    if (chisq < copt) {
      iopt = i;
      copt = chisq;
    }
    ia[i] = false;
    a[i] = 0.0;
  }
  cout << "Optimum 1-var is " << sqrt(copt / static_cast<double>(y.size())) << " for variable "
       << iopt << endl;

  copt = 1e100;
  unsigned i1(0), i2(0);
  for (unsigned i = 0; i < 9; ++i) {
    ia[i] = true;
    for (unsigned j = i + 1; j < 9; ++j) {
      ia[j] = true;
      for (unsigned ii = 0; ii < 3; ++ii) {
        if (ii == 0) {
          alamda = -1;
        }
        if (ii == 2) {
          alamda = 0;
        }
        Check(y.size() < UINT_MAX);
        mrqmin(x, y, sig, static_cast<unsigned>(y.size()), 4U, a, ia, covar, alpha, 9U, chisq,
               model, alamda);
      }
      if (chisq < copt) {
        i1 = i;
        i2 = j;
        copt = chisq;
      }
      ia[j] = false;
      a[j] = 0.0;
    }
    ia[i] = false;
    a[i] = 0.0;
  }
  cout << "Optimum 2-var is " << sqrt(copt / static_cast<double>(y.size())) << " for variables "
       << i1 << i2 << endl;

  copt = 1e100;
  unsigned i3(0);
  for (unsigned i = 0; i < 9; ++i) {
    ia[i] = true;
    for (unsigned j = i + 1; j < 9; ++j) {
      ia[j] = true;
      for (unsigned k = j + 1; k < 9; ++k) {
        ia[k] = true;
        for (unsigned ii = 0; ii < 3; ++ii) {
          if (ii == 0) {
            alamda = -1;
          }
          if (ii == 2) {
            alamda = 0;
          }
          Check(y.size() < UINT_MAX);
          mrqmin(x, y, sig, static_cast<unsigned>(y.size()), 4U, a, ia, covar, alpha, 9U, chisq,
                 model, alamda);
        }
        if (chisq < copt) {
          i1 = i;
          i2 = j;
          i3 = k;
          copt = chisq;
        }
        ia[k] = false;
        a[k] = 0.0;
      }
      ia[j] = false;
      a[j] = 0.0;
    }
    ia[i] = false;
    a[i] = 0.0;
  }
  cout << "Optimum 3-var is " << sqrt(copt / static_cast<double>(y.size())) << " for variables "
       << i1 << i2 << i3 << endl;

  copt = 1e100;
  unsigned i4(0);
  for (unsigned j1 = 0; j1 < 9; ++j1) {
    ia[j1] = true;
    for (unsigned j2 = j1 + 1; j2 < 9; ++j2) {
      ia[j2] = true;
      for (unsigned j3 = j2 + 1; j3 < 9; ++j3) {
        ia[j3] = true;
        for (unsigned j4 = j3 + 1; j4 < 9; ++j4) {
          ia[j4] = true;
          for (unsigned ii = 0; ii < 3; ++ii) {
            if (ii == 0) {
              alamda = -1;
            }
            if (ii == 2) {
              alamda = 0;
            }
            Check(y.size() < UINT_MAX);
            mrqmin(x, y, sig, static_cast<unsigned>(y.size()), 4U, a, ia, covar, alpha, 9U, chisq,
                   model, alamda);
          }
          if (chisq < copt) {
            i1 = j1;
            i2 = j2;
            i3 = j3;
            i4 = j4;
            copt = chisq;
          }
          ia[j4] = false;
          a[j4] = 0.0;
        }
        ia[j3] = false;
        a[j3] = 0.0;
      }
      ia[j2] = false;
      a[j2] = 0.0;
    }
    ia[j1] = false;
    a[j1] = 0.0;
  }
  cout << "Optimum 4-var is " << sqrt(copt / static_cast<double>(y.size())) << " for variables "
       << i1 << i2 << i3 << i4 << endl;

  copt = 1e100;
  unsigned i5(0);
  for (unsigned j1 = 0; j1 < 9; ++j1) {
    ia[j1] = true;
    for (unsigned j2 = j1 + 1; j2 < 9; ++j2) {
      ia[j2] = true;
      for (unsigned j3 = j2 + 1; j3 < 9; ++j3) {
        ia[j3] = true;
        for (unsigned j4 = j3 + 1; j4 < 9; ++j4) {
          ia[j4] = true;
          for (unsigned j5 = j4 + 1; j5 < 9; ++j5) {
            ia[j5] = true;
            for (unsigned ii = 0; ii < 3; ++ii) {
              if (ii == 0) {
                alamda = -1;
              }
              if (ii == 2) {
                alamda = 0;
              }
              Check(y.size() < UINT_MAX);
              mrqmin(x, y, sig, static_cast<unsigned>(y.size()), 4U, a, ia, covar, alpha, 9U, chisq,
                     model, alamda);
            }
            if (chisq < copt) {
              i1 = j1;
              i2 = j2;
              i3 = j3;
              i4 = j4;
              i5 = j5;
              copt = chisq;
            }
            ia[j5] = false;
            a[j5] = 0.0;
          }
          ia[j4] = false;
          a[j4] = 0.0;
        }
        ia[j3] = false;
        a[j3] = 0.0;
      }
      ia[j2] = false;
      a[j2] = 0.0;
    }
    ia[j1] = false;
    a[j1] = 0.0;
  }
  cout << "Optimum 5-var is " << sqrt(copt / static_cast<double>(y.size())) << " for variables "
       << i1 << i2 << i3 << i4 << i5 << endl;

  copt = 1e100;
  unsigned i6(0);
  for (unsigned j1 = 0; j1 < 9; ++j1) {
    ia[j1] = true;
    for (unsigned j2 = j1 + 1; j2 < 9; ++j2) {
      ia[j2] = true;
      for (unsigned j3 = j2 + 1; j3 < 9; ++j3) {
        ia[j3] = true;
        for (unsigned j4 = j3 + 1; j4 < 9; ++j4) {
          ia[j4] = true;
          for (unsigned j5 = j4 + 1; j5 < 9; ++j5) {
            ia[j5] = true;
            for (unsigned j6 = j5 + 1; j6 < 9; ++j6) {
              ia[j6] = true;
              for (unsigned ii = 0; ii < 3; ++ii) {
                if (ii == 0) {
                  alamda = -1;
                }
                if (ii == 2) {
                  alamda = 0;
                }
                Check(y.size() < UINT_MAX);
                mrqmin(x, y, sig, static_cast<unsigned>(y.size()), 4U, a, ia, covar, alpha, 9U,
                       chisq, model, alamda);
              }
              if (chisq < copt) {
                i1 = j1;
                i2 = j2;
                i3 = j3;
                i4 = j4;
                i5 = j5;
                i6 = j6;
                copt = chisq;
              }
              ia[j6] = false;
              a[j6] = 0.0;
            }
            ia[j5] = false;
            a[j5] = 0.0;
          }
          ia[j4] = false;
          a[j4] = 0.0;
        }
        ia[j3] = false;
        a[j3] = 0.0;
      }
      ia[j2] = false;
      a[j2] = 0.0;
    }
    ia[j1] = false;
    a[j1] = 0.0;
  }
  cout << "Optimum 6-var is " << sqrt(copt / static_cast<double>(y.size())) << " for variables "
       << i1 << i2 << i3 << i4 << i5 << i6 << endl;

  copt = 1e100;
  unsigned i7(0);
  for (unsigned j1 = 0; j1 < 9; ++j1) {
    ia[j1] = true;
    for (unsigned j2 = j1 + 1; j2 < 9; ++j2) {
      ia[j2] = true;
      for (unsigned j3 = j2 + 1; j3 < 9; ++j3) {
        ia[j3] = true;
        for (unsigned j4 = j3 + 1; j4 < 9; ++j4) {
          ia[j4] = true;
          for (unsigned j5 = j4 + 1; j5 < 9; ++j5) {
            ia[j5] = true;
            for (unsigned j6 = j5 + 1; j6 < 9; ++j6) {
              ia[j6] = true;
              for (unsigned j7 = j6 + 1; j7 < 9; ++j7) {
                ia[j7] = true;
                for (unsigned ii = 0; ii < 3; ++ii) {
                  if (ii == 0) {
                    alamda = -1;
                  }
                  if (ii == 2) {
                    alamda = 0;
                  }
                  Check(y.size() < UINT_MAX);
                  mrqmin(x, y, sig, static_cast<unsigned>(y.size()), 4U, a, ia, covar, alpha, 9U,
                         chisq, model, alamda);
                }
                if (chisq < copt) {
                  i1 = j1;
                  i2 = j2;
                  i3 = j3;
                  i4 = j4;
                  i5 = j5;
                  i6 = j6;
                  i7 = j7;
                  copt = chisq;
                }
                ia[j7] = false;
                a[j7] = 0.0;
              }
              ia[j6] = false;
              a[j6] = 0.0;
            }
            ia[j5] = false;
            a[j5] = 0.0;
          }
          ia[j4] = false;
          a[j4] = 0.0;
        }
        ia[j3] = false;
        a[j3] = 0.0;
      }
      ia[j2] = false;
      a[j2] = 0.0;
    }
    ia[j1] = false;
    a[j1] = 0.0;
  }
  cout << "Optimum 7-var is " << sqrt(copt / static_cast<double>(y.size())) << " for variables "
       << i1 << i2 << i3 << i4 << i5 << i6 << i7 << endl;

  ia.resize(0);
  ia.resize(9, false);
  //    ia[0] = true;
  //    ia[1] = true;
  //    ia[2] = true;
  ia[3] = true;
  //    ia[4] = true;
  ia[5] = true;
  //    ia[6] = true;
  //    ia[7] = true;
  //    ia[8] = true;
  for (unsigned ii = 0; ii < 3; ++ii) {
    if (ii == 0) {
      alamda = -1;
    }
    if (ii == 2) {
      alamda = 0;
    }
    Check(y.size() < UINT_MAX);
    mrqmin(x, y, sig, static_cast<unsigned>(y.size()), 4U, a, ia, covar, alpha, 9U, chisq, model,
           alamda);
  }

  cout << endl;
  cout << "A = " << a[0] << " +/- " << sqrt(covar[0 + 9 * 0]) << endl;
  cout << "B = " << a[1] << endl;
  cout << "C = " << a[2] << endl;
  cout << "D = " << a[3] << endl;
  cout << "E = " << a[4] << " +/- " << sqrt(covar[4 + 9 * 4]) << endl;
  cout << "F = " << a[5] << " +/- " << sqrt(covar[5 + 9 * 5]) << endl;
  cout << "G = " << a[6] << " +/- " << sqrt(covar[6 + 9 * 6]) << endl;
  cout << "H = " << a[7] << " +/- " << sqrt(covar[7 + 9 * 7]) << endl;
  cout << "J = " << a[8] << " +/- " << sqrt(covar[8 + 9 * 8]) << endl;
  cout << "rms deviation = " << sqrt(chisq / static_cast<double>(y.size())) << endl;

  Check(y.size() < UINT_MAX);
  auto const N = static_cast<unsigned>(y.size());
  vector<double> xx(4);
  vector<double> dyda(9);
  double maxerr = 0;
  unsigned ierr(0);
  for (unsigned i = 0; i < N; ++i) {
    for (unsigned j = 0; j < 4; ++j) {
      xx[j] = x[j + 4 * i];
    }
    double yfit;
    model(xx, a, yfit, dyda);
    double yerr = fabs(yfit - y[i]) / y[i];
    if (yerr > maxerr) {
      maxerr = fabs(yerr);
      ierr = i;
    }
  }
  cout << "Max relative deviation = " << maxerr << " at i = " << ierr << endl;

  if (sqrt(chisq / static_cast<double>(y.size())) < 1) {
    ut.passes("fit is correct");
  } else {
    ut.failure("fit is NOT correct");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstmrqmin(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstmrqmin.cc
//------------------------------------------------------------------------------------------------//
