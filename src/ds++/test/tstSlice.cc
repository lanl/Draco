//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstSlice.cc
 * \author Kent Budge
 * \date   Thu Jul  8 08:02:51 2004
 * \brief  Test the Slice subset container class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Slice.hh"
#include <array>

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstSlice(UnitTest &ut) {
  vector<unsigned> v = {0, 1, 2, 3, 4};
  Slice<vector<unsigned>::iterator> s = slice(v.begin() + 1, 2, 2);
  FAIL_IF_NOT(s.size() == 2);
  FAIL_IF_NOT(s[1] == 3);
  FAIL_IF_NOT(s.begin() < s.end());
  FAIL_IF_NOT(!(s.end() < s.begin()));

  Slice<vector<unsigned>::iterator>::iterator i = s.begin();
  FAIL_IF_NOT(*i == 1);
  FAIL_IF_NOT(i[0] == 1);
  FAIL_IF_NOT(*(i + 1) == 3);
  ++i;
  FAIL_IF_NOT(*i == 3);
  FAIL_IF_NOT(i - s.begin() == 1);
  i++;
  FAIL_IF(i != s.end());

  Slice<vector<unsigned>::iterator>::const_iterator ci = s.begin();
  FAIL_IF_NOT(ci.first() == v.begin() + 1);
  FAIL_IF_NOT(ci.offset() == 0);
  FAIL_IF_NOT(ci.stride() == 2);
  FAIL_IF_NOT(*ci == 1);
  FAIL_IF_NOT(ci[0] == 1);
  FAIL_IF_NOT(*(ci + 1) == 3);
  ++ci;
  FAIL_IF_NOT(*ci == 3);
  FAIL_IF_NOT(s.begin() < ci);
  FAIL_IF_NOT(ci - s.begin() == 1);
  ci++;
  FAIL_IF(ci != s.end());

  Slice<vector<unsigned>::iterator> s2(v.begin(), 3, 2);
  FAIL_IF_NOT(s2.size() == 3);
  FAIL_IF_NOT(s2[1] == 2);

  array<double, 6> da = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  Slice<double *> das(da.data(), 2, 3);
  FAIL_IF_NOT(das.size() == 2);

  vector<double> db_vector(6);
  double *const db = &db_vector[0];
  db[0] = 0;
  Slice<vector<double>::iterator> dbs(db_vector.begin(), 2, 3);
  FAIL_IF_NOT(dbs.size() == 2);

  Slice<vector<unsigned>::iterator> const cs = s;
  FAIL_IF_NOT(cs[1] == 3);
  FAIL_IF_NOT(cs.front() == 1);
  FAIL_IF_NOT(cs.back() == 3);
  FAIL_IF_NOT(cs.begin() < cs.end());
  FAIL_IF_NOT(!(cs.end() < cs.begin()));
  FAIL_IF_NOT(!cs.empty());
  PASSMSG("Done with tests.");
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstSlice(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstSlice.cc
//------------------------------------------------------------------------------------------------//
