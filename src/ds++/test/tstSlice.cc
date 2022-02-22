//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstSlice.cc
 * \author Kent Budge
 * \date   Thu Jul  8 08:02:51 2004
 * \brief  Test the Slice subset container class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
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
  if (s.size() == 2)
    PASSMSG("size of vector Slice is correct");
  else
    FAILMSG("size of vector Slice is NOT correct");
  if (s[1] == 3)
    PASSMSG("indexing of Slice is correct");
  else
    FAILMSG("size of Slice is NOT correct");
  if (s.begin() < s.end())
    PASSMSG("ordering of Slice is correct");
  else
    FAILMSG("ordering of Slice is NOT correct");
  if (!(s.end() < s.begin()))
    PASSMSG("ordering of Slice is correct");
  else
    FAILMSG("ordering of Slice is NOT correct");

  Slice<vector<unsigned>::iterator>::iterator i = s.begin();
  if (*i == 1)
    PASSMSG("deref of begin is correct");
  else
    FAILMSG("deref of begin is NOT correct");
  if (i[0] == 1)
    PASSMSG("operator[] of begin is correct");
  else
    FAILMSG("operator[] of begin is NOT correct");
  if (*(i + 1) == 3)
    PASSMSG("operator+ is correct");
  else
    FAILMSG("operator+ is NOT correct");
  ++i;
  if (*i == 3)
    PASSMSG("deref of ++begin is correct");
  else
    FAILMSG("deref of ++begin is NOT correct");
  if (i - s.begin() == 1)
    PASSMSG("operator- is correct");
  else
    FAILMSG("operator- is NOT correct");
  i++;
  if (i != s.end())
    FAILMSG("++ past end is NOT correct");
  else
    PASSMSG("++ past end is correct");

  Slice<vector<unsigned>::iterator>::const_iterator ci = s.begin();
  if (ci.first() == v.begin() + 1)
    PASSMSG("first is correct");
  else
    FAILMSG("first is NOT correct");
  if (ci.offset() == 0)
    PASSMSG("offset is correct");
  else
    FAILMSG("offset is NOT correct");
  if (ci.stride() == 2)
    PASSMSG("stride is correct");
  else
    FAILMSG("stride is NOT correct");
  if (*ci == 1)
    PASSMSG("deref of begin is correct");
  else
    FAILMSG("deref of begin is NOT correct");
  if (ci[0] == 1)
    PASSMSG("operator[] of begin is correct");
  else
    FAILMSG("operator[] of begin is NOT correct");
  if (*(ci + 1) == 3)
    PASSMSG("operator+ is correct");
  else
    FAILMSG("operator+ is NOT correct");
  ++ci;
  if (*ci == 3)
    PASSMSG("deref of ++begin is correct");
  else
    FAILMSG("deref of ++begin is NOT correct");
  if (s.begin() < ci)
    PASSMSG("correct const iterator ordering");
  else
    FAILMSG("NOT correct const iterator ordering");
  if (ci - s.begin() == 1)
    PASSMSG("operator- is correct");
  else
    FAILMSG("operator- is NOT correct");
  ci++;
  if (ci != s.end())
    FAILMSG("++ past end is NOT correct");
  else
    PASSMSG("++ past end is correct");

  Slice<vector<unsigned>::iterator> s2(v.begin(), 3, 2);
  if (s2.size() == 3)
    PASSMSG("size of Slice is correct");
  else
    FAILMSG("size of Slice is NOT correct");
  if (s2[1] == 2)
    PASSMSG("indexing of Slice is correct");
  else
    FAILMSG("size of Slice is NOT correct");

  array<double, 6> da;
  Slice<double *> das(da.data(), 2, 3);
  if (das.size() == 2)
    PASSMSG("size of Slice is correct");
  else
    FAILMSG("size of Slice is NOT correct");

  vector<double> db_vector(6);
  double *const db = &db_vector[0];
  db[0] = 0;
  Slice<vector<double>::iterator> dbs(db_vector.begin(), 2, 3);
  if (dbs.size() == 2)
    PASSMSG("size of Slice is correct");
  else
    FAILMSG("size of Slice is NOT correct");

  Slice<vector<unsigned>::iterator> const cs = s;
  if (cs[1] == 3)
    PASSMSG("indexing of const Slice is correct");
  else
    FAILMSG("size of const Slice is NOT correct");
  if (cs.front() == 1)
    PASSMSG("front of const Slice is correct");
  else
    FAILMSG("front of const Slice is NOT correct");
  if (cs.back() == 3)
    PASSMSG("back of const Slice is correct");
  else
    FAILMSG("back of const Slice is NOT correct");
  if (cs.begin() < cs.end())
    PASSMSG("ordering of const Slice is correct");
  else
    FAILMSG("ordering of const Slice is NOT correct");
  if (!(cs.end() < cs.begin()))
    PASSMSG("ordering of const Slice is correct");
  else
    FAILMSG("ordering of const Slice is NOT correct");
  if (!cs.empty())
    PASSMSG("emptiness of const Slice is correct");
  else
    FAILMSG("emptiness of const Slice is NOT correct");
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
