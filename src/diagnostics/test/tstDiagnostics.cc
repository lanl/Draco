//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/test/tstDiagnostics.cc
 * \author Thomas M. Evans
 * \date   Fri Dec  9 16:16:27 2005
 * \brief  Diagnostics test.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "diagnostics/Diagnostics.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <array>

using namespace std;
using namespace rtt_diagnostics;
using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_ints(rtt_dsxx::UnitTest &ut) {
  // add an integer quantity
  Diagnostics::integers["A"];
  if (Diagnostics::integers["A"] != 0)
    ITFAILS;

  // now update some values
  Diagnostics::integers["A"] = 21;
  Diagnostics::integers["A"]++;

  if (Diagnostics::integers["A"] != 22)
    ITFAILS;

  // add another
  Diagnostics::integers["B"] = 51;

  if (ut.numFails == 0)
    PASSMSG("Diagnostics integers ok.");
}

//------------------------------------------------------------------------------------------------//

void test_floats(rtt_dsxx::UnitTest &ut) {
  if (Diagnostics::integers["A"] != 22)
    ITFAILS;
  if (Diagnostics::integers["B"] != 51)
    ITFAILS;

  // make a fraction entry
  Diagnostics::doubles["A_of_B"] =
      Diagnostics::integers["A"] / static_cast<double>(Diagnostics::integers["B"]);

  // check it
  if (!soft_equiv(Diagnostics::doubles["A_of_B"], 22.0 / 51.0))
    ITFAILS;

  // erase and check
  if (Diagnostics::doubles.erase("A_of_B") != 1)
    ITFAILS;
  if (Diagnostics::doubles.count("A_of_B") != 0)
    ITFAILS;

  Diagnostics::integers.erase("A");

  if (ut.numFails == 0)
    PASSMSG("Diagnostics doubles ok.");
}

//------------------------------------------------------------------------------------------------//

void test_vectors(rtt_dsxx::UnitTest &ut) {
  Diagnostics::vec_integers["A"];
  Diagnostics::vec_integers["B"];

  if (!Diagnostics::vec_integers["A"].empty())
    ITFAILS;
  if (!Diagnostics::vec_integers["B"].empty())
    ITFAILS;

  Diagnostics::vec_integers["A"].resize(10);
  if (Diagnostics::vec_integers["A"].size() != 10)
    ITFAILS;

  Diagnostics::vec_doubles["B"].resize(2);
  Diagnostics::vec_doubles["B"][0] = 1.1;
  Diagnostics::vec_doubles["B"][1] = 2.4;

  if (Diagnostics::vec_doubles["B"].size() != 2)
    ITFAILS;

  vector<double> ref(2, 1.1);
  ref[1] = 2.4;
  if (!soft_equiv(Diagnostics::vec_doubles["B"].begin(), Diagnostics::vec_doubles["B"].end(),
                  ref.begin(), ref.end()))
    ITFAILS;

  Diagnostics::vec_doubles["B"].clear();
  if (!Diagnostics::vec_integers["B"].empty())
    ITFAILS;

  if (Diagnostics::integers["A"] != 0)
    ITFAILS;
  if (Diagnostics::integers["B"] != 51)
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("Diagnostics vectors ok.");
}

//------------------------------------------------------------------------------------------------//

void test_macro(rtt_dsxx::UnitTest &ut) {
  cout << endl;
  array<int, 4> level = {1, 0, 0, 0};

#ifdef DRACO_DIAGNOSTICS_LEVEL_1
  cout << ">>> Testing Level 1 Block diagnostics." << endl;
  level[1] = 1;
  level[0] = 0;
#endif

#ifdef DRACO_DIAGNOSTICS_LEVEL_2
  cout << ">>> Testing Level 2 Block diagnostics." << endl;
  level[2] = 1;
  level[0] = 0;
#endif

#ifdef DRACO_DIAGNOSTICS_LEVEL_3
  cout << ">>> Testing Level 3 Block diagnostics." << endl;
  level[3] = 1;
  level[0] = 0;
#endif

  cout << endl;

  DIAGNOSTICS_ONE(integers["L1"] = 1);
  DIAGNOSTICS_TWO(integers["L2"] = 1);
  DIAGNOSTICS_THREE(integers["L3"] = 1);

  if (level[0] == 1) {
    if (Diagnostics::integers.count("L1") != 0)
      ITFAILS;
    if (Diagnostics::integers.count("L2") != 0)
      ITFAILS;
    if (Diagnostics::integers.count("L3") != 0)
      ITFAILS;
  }

  if (level[1] == 1)
    if (Diagnostics::integers["L1"] != 1)
      ITFAILS;

  if (level[2] == 1)
    if (Diagnostics::integers["L2"] != 1)
      ITFAILS;

  if (level[3] == 1)
    if (Diagnostics::integers["L3"] != 1)
      ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("Diagnostics macro ok.");
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_ints(ut);
    test_floats(ut);
    test_vectors(ut);
    test_macro(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstDiagnostics.cc
//------------------------------------------------------------------------------------------------//
