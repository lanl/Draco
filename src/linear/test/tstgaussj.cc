//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/test/tstgaussj.cc
 * \author Kent Budge
 * \date   Wed Sep 13 11:46:05 2006
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "linear/gaussj.hh"
#include <complex>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstgaussj(UnitTest &ut) {
  {
    unsigned const N = 4, M = 2;
    vector<double> a, ai(N * N), b(N * M);
    for (unsigned i = 0; i < N; i++) {
      for (unsigned j = 0; j < N; j++) {
        ai[i + N * j] = (i + j + 1) * (i + j + 1);
      }
    }
    ai[N * N - 1]++;
    a = ai;
    for (unsigned i = 0, k = 0; i < N; i++) {
      for (unsigned j = 0; j < M; j++) {
        b[i + N * j] = k++;
      }
    }
    bool count = false;
    gaussj(ai, N, b, M);
    for (unsigned i = 0; i < N; i++) {
      for (unsigned j = 0; j < N; j++) {
        double sum = 0;
        for (unsigned k = 0; k < N; k++) {
          sum += a[i + N * k] * ai[k + N * j];
        }
        count =
            count || (((i == j) && !soft_equiv(sum, 1.0)) || ((i != j) && !soft_equiv(sum, 0.0)));
      }
    }
    for (unsigned i = 0, k = 0; i < N; i++) {
      for (unsigned j = 0; j < M; j++) {
        double sum = 0;
        for (unsigned kp = 0; kp < N; kp++) {
          sum += a[i + N * kp] * b[kp + N * j];
        }
        count = count || !soft_equiv(sum, static_cast<double>(k++));
      }
    }
    if (count)
      ut.failure("did NOT correctly solve equations");
    else
      ut.passes("correctly solved equations");
  }
  {
    unsigned const N = 4;
    vector<vector<double>> a, ai(N);
    vector<double> b(N);
    for (unsigned i = 0; i < N; i++) {
      ai[i].resize(N);
      for (unsigned j = 0; j < N; j++) {
        ai[i][j] = (i + j + 1) * (i + j + 1);
      }
    }
    ai[N - 1][N - 1]++;
    a = ai;
    for (unsigned i = 0, k = 0; i < N; i++) {
      b[i] = k++;
    }
    bool count = false;
    gaussj(ai, b);
    for (unsigned i = 0; i < N; i++) {
      for (unsigned j = 0; j < N; j++) {
        double sum = 0;
        for (unsigned k = 0; k < N; k++) {
          sum += a[i][k] * ai[k][j];
        }
        count =
            count || (((i == j) && !soft_equiv(sum, 1.0)) || ((i != j) && !soft_equiv(sum, 0.0)));
      }
    }
    for (unsigned i = 0, k = 0; i < N; i++) {
      double sum = 0;
      for (unsigned kp = 0; kp < N; kp++) {
        sum += a[i][kp] * b[kp];
      }
      count = count || !soft_equiv(sum, static_cast<double>(k++));
    }
    if (count)
      ut.failure("did NOT correctly solve equations");
    else
      ut.passes("correctly solved equations");
  }
  {
    unsigned const N = 4, M = 2;
    vector<double> a, ai(N * N), b(N * M);
    for (unsigned i = 0; i < N; i++) {
      for (unsigned j = 0; j < N; j++) {
        ai[i + N * j] = (i + j + 1) * (i + j + 1);
      }
    }
    ai[N * N - 1]++;
    a = ai;
    for (unsigned i = 0, k = 0; i < N; i++) {
      for (unsigned j = 0; j < M; j++) {
        b[i + N * j] = k++;
      }
    }
    bool count = false;
    gaussj(ai, N, b, M);
    for (unsigned i = 0; i < N; i++) {
      for (unsigned j = 0; j < N; j++) {
        double sum = 0;
        for (unsigned k = 0; k < N; k++) {
          sum += a[i + N * k] * ai[k + N * j];
        }
        count =
            count || (((i == j) && !soft_equiv(sum, 1.0)) || ((i != j) && !soft_equiv(sum, 0.0)));
      }
    }
    for (unsigned i = 0, k = 0; i < N; i++) {
      for (unsigned j = 0; j < M; j++) {
        double sum = 0;
        for (unsigned kp = 0; kp < N; kp++) {
          sum += a[i + N * kp] * b[kp + N * j];
        }
        count = count || !soft_equiv(sum, static_cast<double>(k++));
      }
    }
    if (count)
      ut.failure("did NOT correctly solve equations");
    else
      ut.passes("correctly solved equations");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstgaussj(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstgaussj.cc
//------------------------------------------------------------------------------------------------//
