//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   units/test/tstFundUnit.cc
 * \author Kelly Thompson
 * \date   Wed Oct  8 13:50:19 2003
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "units/FundUnit.hh"
#include <sstream>

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_construction(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::FundUnit;
  using rtt_units::L_cf;
  using rtt_units::L_labels;
  using rtt_units::L_null;
  using rtt_units::Ltype;
  using std::cout;
  using std::endl;
  using std::ostringstream;
  using std::string;

  FundUnit<Ltype> myLength(L_null, L_cf.data(), L_labels);

  // Is the enum value correct?
  if (myLength.enumVal() == L_null) {
    PASSMSG("Enum Value was set correctly.");
  } else {
    ostringstream msg;
    msg << "Enum Value was not set correctly "
        << "(" << myLength.enumVal() << " != L_null (0) )" << endl;
    FAILMSG(msg.str());
  }

  // Is the conversion factor correct?
  if (soft_equiv(myLength.cf(), 0.0)) {
    PASSMSG("Conversion Factor was set correctly.");
  } else {
    ostringstream msg;
    msg << "Conversion Factor was not set correctly "
        << "(" << myLength.cf() << " != 0.0)" << endl;
    FAILMSG(msg.str());
  }

  // Is the label correct?
  if (myLength.label() == string("NA")) {
    PASSMSG("Unit label was set correctly.");
  } else {
    ostringstream msg;
    msg << "Unit label was not set correctly "
        << "(\"" << myLength.label() << R"(" != "NA"))" << endl; // " fix emacs
    FAILMSG(msg.str());
  }

  //-----------------------------------//

  // Use "long labels" string...

  FundUnit<Ltype> myLength2(rtt_units::L_cm, rtt_units::L_cf.data(), rtt_units::L_long_labels);

  // Is the label correct?
  std::string sentinelValue("centimeter");
  if (myLength2.label() == sentinelValue) {
    PASSMSG("Unit long label was set correctly.");
  } else {
    ostringstream msg;
    msg << "Unit long label was not set correctly "
        << "(\"" << myLength2.label() << "\" != \"" << sentinelValue << "\")" << endl;
    FAILMSG(msg.str());
  }
  // Is the conversion factor correct?
  if (soft_equiv(myLength2.cf(), 100.0)) {
    PASSMSG("Conversion Factor was set correctly.");
  } else {
    ostringstream msg;
    msg << "Conversion Factor was not set correctly "
        << "(" << myLength2.cf() << " != 100.0)" << endl;
    FAILMSG(msg.str());
  }

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_construction(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstFundUnit.cc
//------------------------------------------------------------------------------------------------//
