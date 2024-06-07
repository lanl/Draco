//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   src/units/test/tstUnitSystem.cc
 * \brief  test the UnitSystem class.
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "units/UnitSystem.hh"
#include <sstream>

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_ctor(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::UnitSystem;
  using rtt_units::UnitSystemType;
  using std::cout;
  using std::endl;
  using std::ostringstream;
  using std::string;

  {
    UnitSystem uX4(UnitSystemType().X4());

    // Test the length conversion factor for uX4
    double expVal(100.0);
    if (soft_equiv(uX4.L(), expVal)) {
      ostringstream msg;
      msg << "Conversion Value for X4 UnitSystem is correct." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Conversion Value for X4 UnitSystem is incorrect." << endl
          << "\tThe value returned was " << uX4.L() << " != " << expVal << endl;
      FAILMSG(msg.str());
    }

    // Test the length lable for uX4
    string expValS("cm");
    if (uX4.Lname() == expValS) {
      ostringstream msg;
      msg << "Unit label for X4 UnitSystem is correct." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Unit label for X4 UnitSystem is incorrect." << endl
          << "\tThe value returned was " << uX4.Lname() << " != " << expValS << endl;
      FAILMSG(msg.str());
    }
  }

  //--------------------------------------------------------------------//

  {
    UnitSystem uSI(UnitSystemType().SI());

    // Test the length conversion factor for uSI
    double expVal(1.0);
    if (soft_equiv(uSI.M(), expVal)) {
      ostringstream msg;
      msg << "Conversion Value for SI UnitSystem is correct." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Conversion Value for SI UnitSystem is incorrect." << endl
          << "\tThe value returned was " << uSI.M() << " != " << expVal << endl;
      FAILMSG(msg.str());
    }

    // Test the length lable for uSI
    string expValS("kg");
    if (uSI.Mname() == expValS) {
      ostringstream msg;
      msg << "Unit label for SI UnitSystem is correct." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Unit label for SI UnitSystem is incorrect." << endl
          << "\tThe value returned was " << uSI.Mname() << " != " << expValS << endl;
      FAILMSG(msg.str());
    }
  }

  return;
}

//------------------------------------------------------------------------------------------------//

void test_def_ctor(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::UnitSystem;
  using rtt_units::UnitSystemType;
  using std::cout;
  using std::endl;
  using std::ostringstream;
  using std::string;

  {
    UnitSystem uDef;
    UnitSystem uSI(UnitSystemType().SI());

    // Test the time conversion factor for uDef
    double expVal(uSI.t());
    if (soft_equiv(uDef.t(), expVal)) {
      ostringstream msg;
      msg << "Time Conversion Value for Default UnitSystem is correct." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Time Conversion Value for Default UnitSystem is incorrect." << endl
          << "\tThe value returned was " << uDef.t() << " != " << expVal << endl;
      FAILMSG(msg.str());
    }

    // Test the time lable for uDef
    string expValS(uSI.tname());
    if (uDef.tname() == expValS) {
      ostringstream msg;
      msg << "Time Unit label for Default UnitSystem is correct." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Time Unit label for Default UnitSystem is incorrect." << endl
          << "\tThe value returned was " << uDef.tname() << " != " << expValS << endl;
      FAILMSG(msg.str());
    }

    //------------------------------------------------------------------------------------------------//

    // Test the Temperature conversion factor for uDef
    expVal = uSI.T();
    if (soft_equiv(uDef.T(), expVal)) {
      ostringstream msg;
      msg << "Temperature Conversion Value for Default UnitSystem is correct." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Temperature Conversion Value for Default UnitSystem is incorrect." << endl
          << "\tThe value returned was " << uDef.T() << " != " << expVal << endl;
      FAILMSG(msg.str());
    }

    // Test the tetmperature lable for uDef
    expValS = uSI.Tname();
    if (uDef.Tname() == expValS) {
      ostringstream msg;
      msg << "Temperature Unit label for Default UnitSystem is correct." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Temperature Unit label for Default UnitSystem is incorrect." << endl
          << "\tThe value returned was " << uDef.Tname() << " != " << expValS << endl;
      FAILMSG(msg.str());
    }
  }
}

//--------------------------------------------------------------------//

void test_more_accessors(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::UnitSystem;
  using rtt_units::UnitSystemType;
  using std::cout;
  using std::endl;
  using std::ostringstream;
  using std::string;

  // Create a custom UnitSystem
  UnitSystem myus(UnitSystemType().I(rtt_units::I_amp).A(rtt_units::A_deg).Q(rtt_units::Q_mol));

  // Test the electric current conversion factor.
  double expVal(1.0);
  if (soft_equiv(myus.I(), expVal)) {
    ostringstream msg;
    msg << "Electric Current Conversion Value for Default UnitSystem is "
           "correct."
        << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Electric Current Conversion Value for Default UnitSystem is "
           "incorrect."
        << endl
        << "\tThe value returned was " << myus.I() << " != " << expVal << endl;
    FAILMSG(msg.str());
  }

  // Test the electric current label
  string expValS("Amp");
  if (myus.Iname() == expValS) {
    ostringstream msg;
    msg << " Electric Current Unit label for UnitSystem is correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Electric Current Unit label for UnitSystem is incorrect." << endl
        << "\tThe value returned was " << myus.Iname() << " != " << expValS << endl;
    FAILMSG(msg.str());
  }

  //------------------------------------------------------------------------------------------------//

  // Test the angle conversion factor
  expVal = 57.295779512896171;
  if (soft_equiv(myus.A(), expVal)) {
    ostringstream msg;
    msg << "Angle Conversion Value for Default UnitSystem is correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Angle Conversion Value for Default UnitSystem is incorrect." << endl
        << "\tThe value returned was " << myus.A() << " != " << expVal << endl;
    FAILMSG(msg.str());
  }

  // Test the Angle label
  expValS = "deg";
  if (myus.Aname() == expValS) {
    ostringstream msg;
    msg << " Angle Unit label for UnitSystem is correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Angle Unit label for UnitSystem is incorrect." << endl
        << "\tThe value returned was " << myus.Aname() << " != " << expValS << endl;
    FAILMSG(msg.str());
  }

  //------------------------------------------------------------------------------------------------//

  // Test the quantity conversion factor.
  expVal = 1.0;
  if (soft_equiv(myus.Q(), expVal)) {
    ostringstream msg;
    msg << "Quantity Conversion Value for Default UnitSystem is correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Quantity Conversion Value for Default UnitSystem is incorrect." << endl
        << "\tThe value returned was " << myus.Q() << " != " << expVal << endl;
    FAILMSG(msg.str());
  }

  // Test the quantity label
  expValS = "mol";
  if (myus.Qname() == expValS) {
    ostringstream msg;
    msg << "Quantity Unit label for UnitSystem is correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Quantity Unit label for UnitSystem is incorrect." << endl
        << "\tThe value returned was " << myus.Qname() << " != " << expValS << endl;
    FAILMSG(msg.str());
  }

  return;
}

//------------------------------------------------------------------------------------------------//

void test_aux_accessors(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::UnitSystem;
  using rtt_units::UnitSystemType;
  using std::cout;
  using std::endl;
  using std::ostringstream;
  using std::string;

  // Create SI and X4 UnitSystem
  UnitSystem uX4(UnitSystemType().X4());

  // velocity
  {
    double vel_cf = uX4.L() / uX4.t();
    if (soft_equiv(vel_cf, uX4.v())) {
      ostringstream msg;
      msg << "velocity conversion succeeded." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Velocity conversion failed. Found uX4.v()"
          << " = " << uX4.v() << " != " << vel_cf << endl;
      FAILMSG(msg.str());
    }
  }

  return;
}
//------------------------------------------------------------------------------------------------//

void logic_test(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::UnitSystem;
  using rtt_units::UnitSystemType;
  using std::cout;
  using std::endl;
  using std::ostringstream;
  using std::string;

  // Create SI and X4 UnitSystem
  UnitSystem si;
  UnitSystem uX4(UnitSystemType().X4());

  // Look at auxillary conversion factors.
  {
    double myL = 1;                        // 1.0 cm
    double myLsi = myL * si.L() / uX4.L(); // 0.01 m
    if (soft_equiv(myLsi, 0.01)) {
      ostringstream msg;
      msg << "Length conversion succeeded." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Length conversion failed. Found " << myL << " " << uX4.Lname() << " = " << myLsi
          << " " << si.Lname() << endl;
      FAILMSG(msg.str());
    }
  }

  // Mass
  {
    double myM = 1;                        // 1.0 g
    double myMsi = myM * si.M() / uX4.M(); // 0.001 kg
    if (soft_equiv(myMsi, 0.001)) {
      ostringstream msg;
      msg << "Mass conversion succeeded." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Mass conversion failed. Found " << myM << " " << uX4.Mname() << " = " << myMsi << " "
          << si.Mname() << endl;
      FAILMSG(msg.str());
    }
  }

  // velocity
  {
    double myV = 1.0;                      // 1.0 cm/sh
    double myVsi = myV * si.v() / uX4.v(); // 1e6 m/s
    if (soft_equiv(myVsi, 1.0e+6)) {
      ostringstream msg;
      msg << "Velocity conversion succeeded." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Velocity conversion failed. Found " << myV << " " << uX4.Lname() << "/" << uX4.tname()
          << " = " << myVsi << " " << si.Lname() << "/" << si.tname() << endl;
      FAILMSG(msg.str());
    }
  }

  // acceleration
  {
    double mya = 1.0;                      // 1.0 cm/sh/sh
    double myasi = mya * si.a() / uX4.a(); // 1e14 m/s/s
    if (soft_equiv(myasi, 1.0e+14)) {
      ostringstream msg;
      msg << "Acceleration conversion succeeded." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Acceleration conversion failed. Found " << mya << " " << uX4.Lname() << "/"
          << uX4.tname() << "^2 = " << myasi << " " << si.Lname() << "/" << si.tname() << "^2"
          << endl;
      FAILMSG(msg.str());
    }
  }

  // force
  {
    double myF = 1.0;                      // 1.0 g*cm/sh/sh
    double myFsi = myF * si.f() / uX4.f(); // 1e11 kg*m/s/s
    if (soft_equiv(myFsi, 1.0e+11)) {
      ostringstream msg;
      msg << "Force conversion succeeded." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Force conversion failed. Found " << myF << " " << uX4.Mname() << "*" << uX4.Lname()
          << "/" << uX4.tname() << "^2 = " << myFsi << " " << si.Mname() << "*" << si.Lname() << "/"
          << si.tname() << "^2" << endl;
      FAILMSG(msg.str());
    }
  }

  // energy
  {
    double myE = 1.0;                      // 1.0 g*cm*cm/sh/sh
    double myEsi = myE * si.e() / uX4.e(); // 1e11 kg*m*m/s/s
    if (soft_equiv(myEsi, 1.0e+9)) {
      ostringstream msg;
      msg << "Energy conversion succeeded." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Energy conversion failed. Found " << myE << " " << uX4.Mname() << "*" << uX4.Lname()
          << "^2/" << uX4.tname() << "^2 = " << myEsi << " " << si.Mname() << "*" << si.Lname()
          << "^2/" << si.tname() << "^2" << endl;
      FAILMSG(msg.str());
    }
  }

  // power
  {
    double myP = 1.0;                      // 1.0 g*cm*cm/sh/sh/sh
    double myPsi = myP * si.p() / uX4.p(); // 1e11 kg*m*m/s/s/s
    if (soft_equiv(myPsi, 1.0e+17)) {
      ostringstream msg;
      msg << "Power conversion succeeded." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Power conversion failed. Found " << myP << " " << uX4.Mname() << "*" << uX4.Lname()
          << "^2/" << uX4.tname() << "^3 = " << myPsi << " " << si.Mname() << "*" << si.Lname()
          << "^2/" << si.tname() << "^3" << endl;
      FAILMSG(msg.str());
    }
  }

  return;
}

//------------------------------------------------------------------------------------------------//

void test_eq_op(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::UnitSystem;
  using rtt_units::UnitSystemType;
  using std::cout;
  using std::endl;
  using std::ostringstream;
  using std::string;

  // Create SI and X4 UnitSystem
  UnitSystem si;
  UnitSystem uX4(UnitSystemType().X4());

  if (si == uX4) {
    ostringstream msg;
    msg << "Equality operator failed. " << endl;
    FAILMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Equality operator succeeded." << endl;
    PASSMSG(msg.str());
  }

  if (si != uX4) {
    ostringstream msg;
    msg << "Inequality operator succeeded." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Inequality operator failed. " << endl;
    FAILMSG(msg.str());
  }

  // More thorough branch coverage for equality operator.
  {
    UnitSystem u1(UnitSystemType().L(rtt_units::L_m));
    if (si == u1) {
      ostringstream msg;
      msg << "Equality operator failed for M. " << endl;
      FAILMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Equality operator succeeded for M." << endl;
      PASSMSG(msg.str());
    }
  }
  // More thorough branch coverage for equality operator.
  {
    UnitSystem u1(UnitSystemType().L(rtt_units::L_m).M(rtt_units::M_kg));
    if (si == u1) {
      ostringstream msg;
      msg << "Equality operator failed for t. " << endl;
      FAILMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Equality operator succeeded for t." << endl;
      PASSMSG(msg.str());
    }
  }

  // More thorough branch coverage for equality operator.
  {
    UnitSystem u1(UnitSystemType().L(rtt_units::L_m).M(rtt_units::M_kg).t(rtt_units::t_s));
    if (si == u1) {
      ostringstream msg;
      msg << "Equality operator failed for T. " << endl;
      FAILMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Equality operator succeeded for T." << endl;
      PASSMSG(msg.str());
    }
  }
  // More thorough branch coverage for equality operator.
  {
    UnitSystem u1(
        UnitSystemType().L(rtt_units::L_m).M(rtt_units::M_kg).t(rtt_units::t_s).T(rtt_units::T_K));
    if (si == u1) {
      ostringstream msg;
      msg << "Equality operator failed for I." << endl;
      FAILMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Equality operator succeeded for I." << endl;
      PASSMSG(msg.str());
    }
  }
  // More thorough branch coverage for equality operator.
  {
    UnitSystem u1(UnitSystemType()
                      .L(rtt_units::L_m)
                      .M(rtt_units::M_kg)
                      .t(rtt_units::t_s)
                      .T(rtt_units::T_K)
                      .I(rtt_units::I_amp));
    if (si == u1) {
      ostringstream msg;
      msg << "Equality operator failed for A." << endl;
      FAILMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Equality operator succeeded for A." << endl;
      PASSMSG(msg.str());
    }
  }
  // More thorough branch coverage for equality operator.
  {
    UnitSystem u1(UnitSystemType()
                      .L(rtt_units::L_m)
                      .M(rtt_units::M_kg)
                      .t(rtt_units::t_s)
                      .T(rtt_units::T_K)
                      .I(rtt_units::I_amp)
                      .A(rtt_units::A_rad));
    if (si == u1) {
      ostringstream msg;
      msg << "Equality operator failed for Q." << endl;
      FAILMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Equality operator succeeded for Q." << endl;
      PASSMSG(msg.str());
    }
  }
  // More thorough branch coverage for equality operator.
  {
    UnitSystem u1(UnitSystemType()
                      .L(rtt_units::L_m)
                      .M(rtt_units::M_kg)
                      .t(rtt_units::t_s)
                      .T(rtt_units::T_K)
                      .I(rtt_units::I_amp)
                      .A(rtt_units::A_rad)
                      .Q(rtt_units::Q_mol));
    if (si == u1) {
      ostringstream msg;
      msg << "Equality operator succeeded for SI." << endl;
      PASSMSG(msg.str());
    } else {
      ostringstream msg;
      msg << "Equality operator failed for SI." << endl;
      FAILMSG(msg.str());
    }
  }

  return;
}

//------------------------------------------------------------------------------------------------//

void test_valid_units(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::UnitSystem;
  using rtt_units::UnitSystemType;
  using std::cout;
  using std::endl;
  using std::ostringstream;
  using std::string;

  // Create X4 UnitSystem
  UnitSystem uX4(UnitSystemType().X4());

  if (uX4.validUnits()) {
    ostringstream msg;
    msg << "validUnits() is working for uX4" << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "validUnits() is not working got uX4" << endl;
    FAILMSG(msg.str());
  }

  // Create empty UnitSystem
  UnitSystemType empty;
  UnitSystem uempty(empty);

  if (uempty.validUnits()) {
    ostringstream msg;
    msg << "validUnits() is working for uempty" << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "validUnits() is not working for uempty" << endl;
    FAILMSG(msg.str());
  }

  // Create a faulty UnitSystem (negative cf)
  constexpr std::array<double, 2> my_cf = {0.0, -1.0};

  // Fail with L.
  {
    bool with_dbc(true);
    bool found_assert(false);
    // If DRACO_DBC_LEVE=7, a rtt_dsxx:assertion will be fired when myus is constructed.  If not, we
    // must test the validUnits() member function explicitly.
    try {
      UnitSystem myus(UnitSystemType().L(rtt_units::L_m, my_cf.data()));
      // If DRACO_DBC_LEVE=0 we will get here so call validUnits by manually.
      with_dbc = false;
      if (!myus.validUnits()) {
        ostringstream msg;
        msg << "validUnits() is working for myus (L)." << endl;
        PASSMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "validUnits() is not working for myus (L)." << endl;
        FAILMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion & /*assert*/) {
      ostringstream msg;
      msg << "Expected assertion caught with invalid cf for L." << endl;
      //  << "The message was: " << assert.what() << endl;
      PASSMSG(msg.str());
      found_assert = true;
    } catch (...) {
      FAILMSG("Should not get here: Fail with L.");
    }

    if (with_dbc && !found_assert) {
      ostringstream msg;
      msg << "Expected assertion was not caught with invalid cf for L." << endl;
      FAILMSG(msg.str());
    }
  }

  // Fail with M.
  {
    bool with_dbc(true);
    bool found_assert(false);
    // If DRACO_DBC_LEVEL=7 a rtt_dsxx:assertion will be fired when myus is constructed.  If not, we
    // must test the validUnits() member function explicitly.
    try {
      UnitSystem myus(UnitSystemType().M(rtt_units::M_kg, my_cf.data()));
      // If DRACO_DBC_LEVEL == 0 we will get here so call validUnits by manually.
      with_dbc = false;
      if (myus.validUnits()) {
        ostringstream msg;
        msg << "validUnits() is not working for myus (M)." << endl;
        FAILMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "validUnits() is working for myus (M)." << endl;
        PASSMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion & /*assert*/) {
      ostringstream msg;
      msg << "Expected assertion caught with invalid cf for M." << endl;
      PASSMSG(msg.str());
      found_assert = true;
    } catch (...) {
      FAILMSG("Should not get here: Fail with L.");
    }

    if (with_dbc && !found_assert) {
      ostringstream msg;
      msg << "Expected assertion was not caught with invalid cf for M." << endl;
      FAILMSG(msg.str());
    }
  }

  // Fail with t.
  {
    bool with_dbc(true);
    bool found_assert(false);
    // If DRACO_DBC_LEVEL=7 a rtt_dsxx:assertion will be fired when myus is constructed.  If not, we
    // must test the validUnits() member function explicitly.
    try {
      UnitSystem myus(UnitSystemType().t(rtt_units::t_s, my_cf.data()));
      // If DRACO_DBC_LEVEL == 0 we will get here so call validUnits by manually.
      with_dbc = false;
      if (myus.validUnits()) {
        ostringstream msg;
        msg << "validUnits() is not working for myus (t)." << endl;
        FAILMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "validUnits() is working for myus (t)." << endl;
        PASSMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion & /*assert*/) {
      ostringstream msg;
      msg << "Expected assertion caught with invalid cf for t." << endl;
      PASSMSG(msg.str());
      found_assert = true;
    }
    if (with_dbc && !found_assert) {
      ostringstream msg;
      msg << "Expected assertion was not caught with invalid cf for t." << endl;
      FAILMSG(msg.str());
    }
  }

  // Fail with T
  {
    bool with_dbc(true);
    bool found_assert(false);
    // If DRACO_DBC_LEVEL=7 a rtt_dsxx:assertion will be fired when myus is constructed.  If not, we
    // must test the validUnits() member function explicitly.
    try {
      UnitSystem myus(UnitSystemType().T(rtt_units::T_K, my_cf.data()));
      // If DRACO_DBC_LEVEL == 0 we will get here so call validUnits by manually.
      with_dbc = false;
      if (myus.validUnits()) {
        ostringstream msg;
        msg << "validUnits() is not working for myus (T)." << endl;
        FAILMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "validUnits() is working for myus (T)." << endl;
        PASSMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion & /*assert*/) {
      ostringstream msg;
      msg << "Expected assertion caught with invalid cf for T." << endl;
      PASSMSG(msg.str());
      found_assert = true;
    }
    if (with_dbc && !found_assert) {
      ostringstream msg;
      msg << "Expected assertion was not caught with invalid cf for T." << endl;
      FAILMSG(msg.str());
    }
  }

  // Fail with I
  {
    bool with_dbc(true);
    bool found_assert(false);
    // If DRACO_DBC_LEVEL=7 a rtt_dsxx:assertion will be fired when myus is constructed.  If not, we
    // must test the validUnits() member function explicitly.
    try {
      UnitSystem myus(UnitSystemType().I(rtt_units::I_amp, my_cf.data()));
      // If DRACO_DBC_LEVEL == 0 we will get here so call validUnits by manually.
      with_dbc = false;
      if (myus.validUnits()) {
        ostringstream msg;
        msg << "validUnits() is not working for myus (I)." << endl;
        FAILMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "validUnits() is working for myus (I)." << endl;
        PASSMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion & /*assert*/) {
      ostringstream msg;
      msg << "Expected assertion caught with invalid cf for I." << endl;
      PASSMSG(msg.str());
      found_assert = true;
    }
    if (with_dbc && !found_assert) {
      ostringstream msg;
      msg << "Expected assertion was not caught with invalid cf for I." << endl;
      FAILMSG(msg.str());
    }
  }

  // Fail with A
  {
    bool with_dbc(true);
    bool found_assert(false);
    // If DRACO_DBC_LEVEL=7 a rtt_dsxx:assertion will be fired when myus is constructed.  If not, we
    // must test the validUnits() member function explicitly.
    try {
      UnitSystem myus(UnitSystemType().A(rtt_units::A_rad, my_cf.data()));
      // If DRACO_DBC_LEVEL == 0 we will get here so call validUnits by manually.
      with_dbc = false;
      if (myus.validUnits()) {
        ostringstream msg;
        msg << "validUnits() is not working for myus (A)." << endl;
        FAILMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "validUnits() is working for myus (A)." << endl;
        PASSMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion & /*assert*/) {
      ostringstream msg;
      msg << "Expected assertion caught with invalid cf for A." << endl;
      PASSMSG(msg.str());
      found_assert = true;
    }
    if (with_dbc && !found_assert) {
      ostringstream msg;
      msg << "Expected assertion was not caught with invalid cf for A." << endl;
      FAILMSG(msg.str());
    }
  }

  // Fail with Q
  {
    bool with_dbc(true);
    bool found_assert(false);
    // If DRACO_DBC_LEVEL=7 a rtt_dsxx:assertion will be fired when myus is constructed.  If not, we
    // must test the validUnits() member function explicitly.
    try {
      UnitSystem myus(UnitSystemType().Q(rtt_units::Q_mol, my_cf.data()));
      // If DRACO_DBC_LEVEL == 0 we will get here so call validUnits by manually.
      with_dbc = false;
      if (myus.validUnits()) {
        ostringstream msg;
        msg << "validUnits() is not working for myus (Q)." << endl;
        FAILMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "validUnits() is working for myus (Q)." << endl;
        PASSMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion & /*assert*/) {
      ostringstream msg;
      msg << "Expected assertion caught with invalid cf for Q." << endl;
      PASSMSG(msg.str());
      found_assert = true;
    }
    if (with_dbc && !found_assert) {
      ostringstream msg;
      msg << "Expected assertion was not caught with invalid cf for Q." << endl;
      FAILMSG(msg.str());
    }
  }

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_ctor(ut);
    test_def_ctor(ut);
    test_more_accessors(ut);
    test_aux_accessors(ut);
    logic_test(ut);
    test_eq_op(ut);
    test_valid_units(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstUnitSystem.cc
//------------------------------------------------------------------------------------------------//
