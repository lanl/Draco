//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstData_Table.cc
 * \author Paul Henning
 * \brief  DBC_Ptr tests.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Data_Table.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include <array>

using namespace std;
using rtt_dsxx::Data_Table;

//------------------------------------------------------------------------------------------------//
void test_array(rtt_dsxx::UnitTest &ut) {
  constexpr array<int, 3> myarray = {3, 100, 12};

  bool caught = false;
  try {
    Data_Table<int> dt(&(myarray.front()), &(myarray.back()));
  } catch (rtt_dsxx::assertion & /* error */) {
    caught = true;
  }
  FAIL_IF(caught);

  caught = false;
  try {
    Data_Table<int> dt(&(myarray.back()), &(myarray.front()));
  } catch (rtt_dsxx::assertion & /* error */) {
    caught = true;
  }
  FAIL_IF_NOT(caught);

  Data_Table<int> dt(&(myarray.front()), &(myarray.back()));

  caught = false;
  try {
    FAIL_IF_NOT(dt.size() == 3);
    FAIL_IF_NOT(dt[0] == myarray[0]);
    FAIL_IF_NOT(dt[1] == myarray[1]);
    FAIL_IF_NOT(dt[2] == myarray[2]);
    FAIL_IF_NOT(dt.front() == myarray[0]);
    FAIL_IF_NOT(dt.back() == myarray[2]);
    FAIL_IF_NOT(dt.begin() == &(myarray.front()));
    FAIL_IF_NOT(dt.end() == &(myarray.back()));
    FAIL_IF_NOT(dt.access() == &dt[0]);

    {
      // Suppress clang-tidy warning since we don't care about performance and really want to test
      // this form of the ctor.
      Data_Table<int> dt3(dt); // NOLINT(performance-unnecessary-copy-initialization)
      FAIL_IF_NOT(dt3.size() == dt.size());
      FAIL_IF_NOT(dt3.begin() == dt.begin());
      FAIL_IF_NOT(dt3.end() == dt.end());
      FAIL_IF_NOT(dt3.front() == dt.front());
      FAIL_IF_NOT(dt3.back() == dt.back());
    }

    {
      Data_Table<int> dt3;
      FAIL_IF_NOT(dt3.size() == 0);
      dt3 = dt;
      FAIL_IF_NOT(dt3.size() == dt.size());
      FAIL_IF_NOT(dt3.begin() == dt.begin());
      FAIL_IF_NOT(dt3.end() == dt.end());
      FAIL_IF_NOT(dt3.front() == dt.front());
      FAIL_IF_NOT(dt3.back() == dt.back());
    }
  } catch (rtt_dsxx::assertion & /* error */) {
    caught = true;
  }
  FAIL_IF(caught);

#ifdef DEBUG
  /*
  GCC will issue a warning at compile time for a Release build (with -ftree-vrp, which is enabled by
  default with -O2 or higher).  The warning appears because the size of dt is known at compile time.
*/
  caught = false;
  try {
    std::cout << dt[3];
  } catch (rtt_dsxx::assertion & /* error */) {
    caught = true;
  }
  FAIL_IF_NOT(caught);
#endif

  if (ut.numFails == 0)
    PASSMSG("test_array");
  else
    FAILMSG("test_array FAILED!");

  return;
}

//------------------------------------------------------------------------------------------------//
void test_scalar(rtt_dsxx::UnitTest &ut) {
  Data_Table<int> dt(32);

  bool caught = false;
  try {
    FAIL_IF_NOT(dt.size() == 1);
    FAIL_IF_NOT(dt[0] == 32);

    Data_Table<int> dt2;
    dt2 = dt;
    FAIL_IF_NOT(dt[0] == dt2[0]);
    FAIL_IF(&(dt[0]) == &(dt2[0]));
    FAIL_IF_NOT(dt.front() == 32);
    FAIL_IF_NOT(dt.back() == 32);
    FAIL_IF_NOT(*(dt.begin()) == 32);
    FAIL_IF_NOT(dt2.size() == 1);

    Data_Table<int> dt3(dt2);
    FAIL_IF_NOT(dt[0] == dt3[0]);
    FAIL_IF(&(dt[0]) == &(dt3[0]));

    dt = dt3;
    FAIL_IF_NOT(dt[0] == dt3[0]);
    FAIL_IF(&(dt[0]) == &(dt3[0]));
  } catch (rtt_dsxx::assertion & /* error */) {
    caught = true;
  }
  FAIL_IF(caught);

#ifdef DEBUG
  // GCC in RELEASE mode will refuse to compile this code.
  caught = false;
  try {
    // access one past end of valid range.
    std::cout << dt[1];
  } catch (rtt_dsxx::assertion & /* error */) {
    caught = true;
  }
  FAIL_IF_NOT(caught);

#endif

  if (ut.numFails == 0)
    PASSMSG("test_scalar");
  else
    FAILMSG("test_scalar FAILED!");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  if (ut.dbcOn() && !ut.dbcNothrow()) {
    try {
      // >>> UNIT TESTS
      test_array(ut);
      test_scalar(ut);
    } catch (rtt_dsxx::assertion &error) {
      cout << "ERROR: While testing tstData_Table_Ptr, " << error.what() << endl;
      ut.numFails++;
    }

    catch (...) {
      cout << "ERROR: While testing " << argv[0] << ", An unknown exception was thrown" << endl;
      ut.numFails++;
    }
  } else {
    PASSMSG("Unit tests only works if DBC is on and the DBC nothrow option is off.");
  }
  return ut.numFails;
}

//------------------------------------------------------------------------------------------------//
// end of tstData_Table.cc
//------------------------------------------------------------------------------------------------//
