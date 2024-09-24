//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   src/units/test/tstUnitSystem.cc
 * \brief  test the UnitSystem class.
 *\note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Assert.hh"
#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "units/UnitSystem.hh"
#include "units/UnitSystemEnums.hh"
#include "units/UnitSystemType.hh"
#include <array>
#include <iostream>
#include <sstream>
#include <string>

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_ctor(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::UnitSystem;
  using rtt_units::UnitSystemType;
  using std::cout;
  using std::ostringstream;
  using std::string;

  {
    UnitSystem const uX4(UnitSystemType().X4());

    // Test the length conversion factor for uX4
    double const expVal(100.0);
    if (soft_equiv(uX4.L(), expVal)) {
      PASSMSG("Conversion Value for X4 UnitSystem is correct.\n");
    } else {
      ostringstream msg;
      msg << "Conversion Value for X4 UnitSystem is incorrect.\n\tThe value returned was "
          << uX4.L() << " != " << expVal << '\n';
      FAILMSG(msg.str());
    }

    // Test the length lable for uX4
    string const expValS("cm");
    if (uX4.Lname() == expValS) {
      PASSMSG("Unit label for X4 UnitSystem is correct.\n");
    } else {
      ostringstream msg;
      msg << "Unit label for X4 UnitSystem is incorrect.\n\tThe value returned was " << uX4.Lname()
          << " != " << expValS << '\n';
      FAILMSG(msg.str());
    }
  }

  //----------------------------------------------------------------------------------------------//

  {
    UnitSystem const uSI(UnitSystemType().SI());

    // Test the length conversion factor for uSI
    double const expVal(1.0);
    if (soft_equiv(uSI.M(), expVal)) {
      PASSMSG("Conversion Value for SI UnitSystem is correct.\n");
    } else {
      ostringstream msg;
      msg << "Conversion Value for SI UnitSystem is incorrect.\n\tThe value returned was "
          << uSI.M() << " != " << expVal << '\n';
      FAILMSG(msg.str());
    }

    // Test the length lable for uSI
    string const expValS("kg");
    if (uSI.Mname() == expValS) {
      PASSMSG("Unit label for SI UnitSystem is correct.\n");
    } else {
      ostringstream msg;
      msg << "Unit label for SI UnitSystem is incorrect.\n\tThe value returned was " << uSI.Mname()
          << " != " << expValS << '\n';
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
  using std::ostringstream;
  using std::string;

  {
    UnitSystem const uDef;
    UnitSystem const uSI(UnitSystemType().SI());

    // Test the time conversion factor for uDef
    double expVal(uSI.t());
    if (soft_equiv(uDef.t(), expVal)) {
      PASSMSG("Time Conversion Value for Default UnitSystem is correct.\n");
    } else {
      ostringstream msg;
      msg << "Time Conversion Value for Default UnitSystem is incorrect.\n"
          << "\tThe value returned was " << uDef.t() << " != " << expVal << '\n';
      FAILMSG(msg.str());
    }

    // Test the time lable for uDef
    string expValS(uSI.tname());
    if (uDef.tname() == expValS) {
      PASSMSG("Time Unit label for Default UnitSystem is correct.\n");
    } else {
      ostringstream msg;
      msg << "Time Unit label for Default UnitSystem is incorrect.\n"
          << "\tThe value returned was " << uDef.tname() << " != " << expValS << '\n';
      FAILMSG(msg.str());
    }

    //--------------------------------------------------------------------------------------------//

    // Test the Temperature conversion factor for uDef
    expVal = uSI.T();
    if (soft_equiv(uDef.T(), expVal)) {
      PASSMSG("Temperature Conversion Value for Default UnitSystem is correct.\n");
    } else {
      ostringstream msg;
      msg << "Temperature Conversion Value for Default UnitSystem is incorrect.\n"
          << "\tThe value returned was " << uDef.T() << " != " << expVal << '\n';
      FAILMSG(msg.str());
    }

    // Test the tetmperature lable for uDef
    expValS = uSI.Tname();
    if (uDef.Tname() == expValS) {
      PASSMSG("Temperature Unit label for Default UnitSystem is correct.\n");
    } else {
      ostringstream msg;
      msg << "Temperature Unit label for Default UnitSystem is incorrect.\n"
          << "\tThe value returned was " << uDef.Tname() << " != " << expValS << '\n';
      FAILMSG(msg.str());
    }
  }
}

//------------------------------------------------------------------------------------------------//

void test_more_accessors(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::UnitSystem;
  using rtt_units::UnitSystemType;
  using std::cout;
  using std::ostringstream;
  using std::string;

  // Create a custom UnitSystem
  UnitSystem const myus(
      UnitSystemType().I(rtt_units::I_amp).A(rtt_units::A_deg).Q(rtt_units::Q_mol));

  // Test the electric current conversion factor.
  double expVal(1.0);
  if (soft_equiv(myus.I(), expVal)) {
    PASSMSG("Electric Current Conversion Value for Default UnitSystem is correct.\n");
  } else {
    ostringstream msg;
    msg << "Electric Current Conversion Value for Default UnitSystem is incorrect.\n"
        << "\tThe value returned was " << myus.I() << " != " << expVal << '\n';
    FAILMSG(msg.str());
  }

  // Test the electric current label
  string expValS("Amp");
  if (myus.Iname() == expValS) {
    PASSMSG(" Electric Current Unit label for UnitSystem is correct.\n");
  } else {
    ostringstream msg;
    msg << "Electric Current Unit label for UnitSystem is incorrect.\n"
        << "\tThe value returned was " << myus.Iname() << " != " << expValS << '\n';
    FAILMSG(msg.str());
  }

  //----------------------------------------------------------------------------------------------//

  // Test the angle conversion factor
  expVal = 57.295779512896171;
  if (soft_equiv(myus.A(), expVal)) {
    PASSMSG("Angle Conversion Value for Default UnitSystem is correct.\n");
  } else {
    ostringstream msg;
    msg << "Angle Conversion Value for Default UnitSystem is incorrect.\n"
        << "\tThe value returned was " << myus.A() << " != " << expVal << '\n';
    FAILMSG(msg.str());
  }

  // Test the Angle label
  expValS = "deg";
  if (myus.Aname() == expValS) {
    PASSMSG(" Angle Unit label for UnitSystem is correct.\n");
  } else {
    ostringstream msg;
    msg << "Angle Unit label for UnitSystem is incorrect.\n"
        << "\tThe value returned was " << myus.Aname() << " != " << expValS << '\n';
    FAILMSG(msg.str());
  }

  //----------------------------------------------------------------------------------------------//

  // Test the quantity conversion factor.
  expVal = 1.0;
  if (soft_equiv(myus.Q(), expVal)) {
    PASSMSG("Quantity Conversion Value for Default UnitSystem is correct.\n");
  } else {
    ostringstream msg;
    msg << "Quantity Conversion Value for Default UnitSystem is incorrect.\n"
        << "\tThe value returned was " << myus.Q() << " != " << expVal << '\n';
    FAILMSG(msg.str());
  }

  // Test the quantity label
  expValS = "mol";
  if (myus.Qname() == expValS) {
    PASSMSG("Quantity Unit label for UnitSystem is correct.\n");
  } else {
    ostringstream msg;
    msg << "Quantity Unit label for UnitSystem is incorrect.\n"
        << "\tThe value returned was " << myus.Qname() << " != " << expValS << '\n';
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
  using std::ostringstream;
  using std::string;

  // Create SI and X4 UnitSystem
  UnitSystem const uX4(UnitSystemType().X4());

  // velocity
  {
    double const vel_cf = uX4.L() / uX4.t();
    if (soft_equiv(vel_cf, uX4.v())) {
      PASSMSG("velocity conversion succeeded.\n");
    } else {
      ostringstream msg;
      msg << "Velocity conversion failed. Found uX4.v()"
          << " = " << uX4.v() << " != " << vel_cf << '\n';
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
  using std::ostringstream;
  using std::string;

  // Create SI and X4 UnitSystem
  UnitSystem const si;
  UnitSystem const uX4(UnitSystemType().X4());

  // Look at auxillary conversion factors.
  {
    double const myL = 1;                        // 1.0 cm
    double const myLsi = myL * si.L() / uX4.L(); // 0.01 m
    if (soft_equiv(myLsi, 0.01)) {
      PASSMSG("Length conversion succeeded.\n");
    } else {
      ostringstream msg;
      msg << "Length conversion failed. Found " << myL << " " << uX4.Lname() << " = " << myLsi
          << " " << si.Lname() << '\n';
      FAILMSG(msg.str());
    }
  }

  // Mass
  {
    double const myM = 1;                        // 1.0 g
    double const myMsi = myM * si.M() / uX4.M(); // 0.001 kg
    if (soft_equiv(myMsi, 0.001)) {
      PASSMSG("Mass conversion succeeded.\n");
    } else {
      ostringstream msg;
      msg << "Mass conversion failed. Found " << myM << " " << uX4.Mname() << " = " << myMsi << " "
          << si.Mname() << '\n';
      FAILMSG(msg.str());
    }
  }

  // velocity
  {
    double const myV = 1.0;                      // 1.0 cm/sh
    double const myVsi = myV * si.v() / uX4.v(); // 1e6 m/s
    if (soft_equiv(myVsi, 1.0e+6)) {
      PASSMSG("Velocity conversion succeeded.\n");
    } else {
      ostringstream msg;
      msg << "Velocity conversion failed. Found " << myV << " " << uX4.Lname() << "/" << uX4.tname()
          << " = " << myVsi << " " << si.Lname() << "/" << si.tname() << '\n';
      FAILMSG(msg.str());
    }
  }

  // acceleration
  {
    double const mya = 1.0;                      // 1.0 cm/sh/sh
    double const myasi = mya * si.a() / uX4.a(); // 1e14 m/s/s
    if (soft_equiv(myasi, 1.0e+14)) {
      PASSMSG("Acceleration conversion succeeded.\n");
    } else {
      ostringstream msg;
      msg << "Acceleration conversion failed. Found " << mya << " " << uX4.Lname() << "/"
          << uX4.tname() << "^2 = " << myasi << " " << si.Lname() << "/" << si.tname() << "^2"
          << '\n';
      FAILMSG(msg.str());
    }
  }

  // force
  {
    double const myF = 1.0;                      // 1.0 g*cm/sh/sh
    double const myFsi = myF * si.f() / uX4.f(); // 1e11 kg*m/s/s
    if (soft_equiv(myFsi, 1.0e+11)) {
      PASSMSG("Force conversion succeeded.\n");
    } else {
      ostringstream msg;
      msg << "Force conversion failed. Found " << myF << " " << uX4.Mname() << "*" << uX4.Lname()
          << "/" << uX4.tname() << "^2 = " << myFsi << " " << si.Mname() << "*" << si.Lname() << "/"
          << si.tname() << "^2\n";
      FAILMSG(msg.str());
    }
  }

  // energy
  {
    double const myE = 1.0;                      // 1.0 g*cm*cm/sh/sh
    double const myEsi = myE * si.e() / uX4.e(); // 1e11 kg*m*m/s/s
    if (soft_equiv(myEsi, 1.0e+9)) {
      PASSMSG("Energy conversion succeeded.\n");
    } else {
      ostringstream msg;
      msg << "Energy conversion failed. Found " << myE << " " << uX4.Mname() << "*" << uX4.Lname()
          << "^2/" << uX4.tname() << "^2 = " << myEsi << " " << si.Mname() << "*" << si.Lname()
          << "^2/" << si.tname() << "^2\n";
      FAILMSG(msg.str());
    }
  }

  // power
  {
    double const myP = 1.0;                      // 1.0 g*cm*cm/sh/sh/sh
    double const myPsi = myP * si.p() / uX4.p(); // 1e11 kg*m*m/s/s/s
    if (soft_equiv(myPsi, 1.0e+17)) {
      PASSMSG("Power conversion succeeded.\n");
    } else {
      ostringstream msg;
      msg << "Power conversion failed. Found " << myP << " " << uX4.Mname() << "*" << uX4.Lname()
          << "^2/" << uX4.tname() << "^3 = " << myPsi << " " << si.Mname() << "*" << si.Lname()
          << "^2/" << si.tname() << "^3\n";
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
  using std::ostringstream;
  using std::string;

  // Create SI and X4 UnitSystem
  UnitSystem const si;
  UnitSystem const uX4(UnitSystemType().X4());

  if (si == uX4) {
    FAILMSG("Equality operator failed.\n");
  } else {
    PASSMSG("Equality operator succeeded.\n");
  }

  if (si != uX4) {
    PASSMSG("Inequality operator succeeded.\n");
  } else {
    FAILMSG("Inequality operator failed.\n");
  }

  // More thorough branch coverage for equality operator.
  {
    UnitSystem const u1(UnitSystemType().L(rtt_units::L_m));
    if (si == u1) {
      FAILMSG("Equality operator failed for M.\n");
    } else {
      PASSMSG("Equality operator succeeded for M.\n");
    }
  }
  // More thorough branch coverage for equality operator.
  {
    UnitSystem const u1(UnitSystemType().L(rtt_units::L_m).M(rtt_units::M_kg));
    if (si == u1) {
      FAILMSG("Equality operator failed for t.\n");
    } else {
      PASSMSG("Equality operator succeeded for t.\n");
    }
  }

  // More thorough branch coverage for equality operator.
  {
    UnitSystem const u1(UnitSystemType().L(rtt_units::L_m).M(rtt_units::M_kg).t(rtt_units::t_s));
    if (si == u1) {
      FAILMSG("Equality operator failed for T.\n");
    } else {
      PASSMSG("Equality operator succeeded for T.\n");
    }
  }
  // More thorough branch coverage for equality operator.
  {
    UnitSystem const u1(
        UnitSystemType().L(rtt_units::L_m).M(rtt_units::M_kg).t(rtt_units::t_s).T(rtt_units::T_K));
    if (si == u1) {
      FAILMSG("Equality operator failed for I.\n");
    } else {
      PASSMSG("Equality operator succeeded for I.\n");
    }
  }
  // More thorough branch coverage for equality operator.
  {
    UnitSystem const u1(UnitSystemType()
                            .L(rtt_units::L_m)
                            .M(rtt_units::M_kg)
                            .t(rtt_units::t_s)
                            .T(rtt_units::T_K)
                            .I(rtt_units::I_amp));
    if (si == u1) {
      FAILMSG("Equality operator failed for A.\n");
    } else {
      PASSMSG("Equality operator succeeded for A.\n");
    }
  }
  // More thorough branch coverage for equality operator.
  {
    UnitSystem const u1(UnitSystemType()
                            .L(rtt_units::L_m)
                            .M(rtt_units::M_kg)
                            .t(rtt_units::t_s)
                            .T(rtt_units::T_K)
                            .I(rtt_units::I_amp)
                            .A(rtt_units::A_rad));
    if (si == u1) {
      FAILMSG("Equality operator failed for Q.\n");
    } else {
      PASSMSG("Equality operator succeeded for Q.\n");
    }
  }
  // More thorough branch coverage for equality operator.
  {
    UnitSystem const u1(UnitSystemType()
                            .L(rtt_units::L_m)
                            .M(rtt_units::M_kg)
                            .t(rtt_units::t_s)
                            .T(rtt_units::T_K)
                            .I(rtt_units::I_amp)
                            .A(rtt_units::A_rad)
                            .Q(rtt_units::Q_mol));
    if (si == u1) {
      PASSMSG("Equality operator succeeded for SI.\n");
    } else {
      FAILMSG("Equality operator failed for SI.\n");
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
  using std::ostringstream;
  using std::string;

  // Create X4 UnitSystem
  UnitSystem const uX4(UnitSystemType().X4());

  if (uX4.validUnits()) {
    PASSMSG("validUnits() is working for uX4\n");
  } else {
    FAILMSG("validUnits() is not working got uX4\n");
  }

  // Create empty UnitSystem
  UnitSystemType const empty;
  UnitSystem const uempty(empty);

  if (uempty.validUnits()) {
    PASSMSG("validUnits() is working for uempty\n");
  } else {
    FAILMSG("validUnits() is not working for uempty\n");
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
      UnitSystem const myus(UnitSystemType().L(rtt_units::L_m, my_cf.data()));
      // If DRACO_DBC_LEVE=0 we will get here so call validUnits by manually.
      with_dbc = false;
      if (!myus.validUnits()) {
        PASSMSG("validUnits() is working for myus (L).\n");
      } else {
        FAILMSG("validUnits() is not working for myus (L).\n");
      }
    } catch (rtt_dsxx::assertion & /*assert*/) {
      ostringstream msg;
      msg << "Expected assertion caught with invalid cf for L.\n";
      //  << "The message was: " << assert.what() << '\n';
      PASSMSG(msg.str());
      found_assert = true;
    } catch (...) {
      FAILMSG("Should not get here: Fail with L.");
    }

    if (with_dbc && !found_assert) {
      FAILMSG("Expected assertion was not caught with invalid cf for L.\n");
    }
  }

  // Fail with M.
  {
    bool with_dbc(true);
    bool found_assert(false);
    // If DRACO_DBC_LEVEL=7 a rtt_dsxx:assertion will be fired when myus is constructed.  If not, we
    // must test the validUnits() member function explicitly.
    try {
      UnitSystem const myus(UnitSystemType().M(rtt_units::M_kg, my_cf.data()));
      // If DRACO_DBC_LEVEL == 0 we will get here so call validUnits by manually.
      with_dbc = false;
      if (myus.validUnits()) {
        FAILMSG("validUnits() is not working for myus (M).\n");
      } else {
        PASSMSG("validUnits() is working for myus (M).\n");
      }
    } catch (rtt_dsxx::assertion & /*assert*/) {
      PASSMSG("Expected assertion caught with invalid cf for M.\n");
      found_assert = true;
    } catch (...) {
      FAILMSG("Should not get here: Fail with L.");
    }

    if (with_dbc && !found_assert) {
      ostringstream msg;
      msg << "Expected assertion was not caught with invalid cf for M.\n";
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
      UnitSystem const myus(UnitSystemType().t(rtt_units::t_s, my_cf.data()));
      // If DRACO_DBC_LEVEL == 0 we will get here so call validUnits by manually.
      with_dbc = false;
      if (myus.validUnits()) {
        FAILMSG("validUnits() is not working for myus (t).\n");
      } else {
        PASSMSG("validUnits() is working for myus (t).\n");
      }
    } catch (rtt_dsxx::assertion & /*assert*/) {
      PASSMSG("Expected assertion caught with invalid cf for t.\n");
      found_assert = true;
    }
    if (with_dbc && !found_assert) {
      FAILMSG("Expected assertion was not caught with invalid cf for t.\n");
    }
  }

  // Fail with T
  {
    bool with_dbc(true);
    bool found_assert(false);
    // If DRACO_DBC_LEVEL=7 a rtt_dsxx:assertion will be fired when myus is constructed.  If not, we
    // must test the validUnits() member function explicitly.
    try {
      UnitSystem const myus(UnitSystemType().T(rtt_units::T_K, my_cf.data()));
      // If DRACO_DBC_LEVEL == 0 we will get here so call validUnits by manually.
      with_dbc = false;
      if (myus.validUnits()) {
        FAILMSG("validUnits() is not working for myus (T).\n");
      } else {
        PASSMSG("validUnits() is working for myus (T).\n");
      }
    } catch (rtt_dsxx::assertion & /*assert*/) {
      PASSMSG("Expected assertion caught with invalid cf for T.\n");
      found_assert = true;
    }
    if (with_dbc && !found_assert) {
      FAILMSG("Expected assertion was not caught with invalid cf for T.\n");
    }
  }

  // Fail with I
  {
    bool with_dbc(true);
    bool found_assert(false);
    // If DRACO_DBC_LEVEL=7 a rtt_dsxx:assertion will be fired when myus is constructed.  If not, we
    // must test the validUnits() member function explicitly.
    try {
      UnitSystem const myus(UnitSystemType().I(rtt_units::I_amp, my_cf.data()));
      // If DRACO_DBC_LEVEL == 0 we will get here so call validUnits by manually.
      with_dbc = false;
      if (myus.validUnits()) {
        FAILMSG("validUnits() is not working for myus (I).\n");
      } else {
        PASSMSG("validUnits() is working for myus (I).\n");
      }
    } catch (rtt_dsxx::assertion & /*assert*/) {
      PASSMSG("Expected assertion caught with invalid cf for I.\n");
      found_assert = true;
    }
    if (with_dbc && !found_assert) {
      FAILMSG("Expected assertion was not caught with invalid cf for I.\n");
    }
  }

  // Fail with A
  {
    bool with_dbc(true);
    bool found_assert(false);
    // If DRACO_DBC_LEVEL=7 a rtt_dsxx:assertion will be fired when myus is constructed.  If not, we
    // must test the validUnits() member function explicitly.
    try {
      UnitSystem const myus(UnitSystemType().A(rtt_units::A_rad, my_cf.data()));
      // If DRACO_DBC_LEVEL == 0 we will get here so call validUnits by manually.
      with_dbc = false;
      if (myus.validUnits()) {
        FAILMSG("validUnits() is not working for myus (A).\n");
      } else {
        PASSMSG("validUnits() is working for myus (A).\n");
      }
    } catch (rtt_dsxx::assertion & /*assert*/) {
      PASSMSG("Expected assertion caught with invalid cf for A.\n");
      found_assert = true;
    }
    if (with_dbc && !found_assert) {
      FAILMSG("Expected assertion was not caught with invalid cf for A.\n");
    }
  }

  // Fail with Q
  {
    bool with_dbc(true);
    bool found_assert(false);
    // If DRACO_DBC_LEVEL=7 a rtt_dsxx:assertion will be fired when myus is constructed.  If not, we
    // must test the validUnits() member function explicitly.
    try {
      UnitSystem const myus(UnitSystemType().Q(rtt_units::Q_mol, my_cf.data()));
      // If DRACO_DBC_LEVEL == 0 we will get here so call validUnits by manually.
      with_dbc = false;
      if (myus.validUnits()) {
        FAILMSG("validUnits() is not working for myus (Q).\n");
      } else {
        PASSMSG("validUnits() is working for myus (Q).\n");
      }
    } catch (rtt_dsxx::assertion & /*assert*/) {
      PASSMSG("Expected assertion caught with invalid cf for Q.\n");
      found_assert = true;
    }
    if (with_dbc && !found_assert) {
      FAILMSG("Expected assertion was not caught with invalid cf for Q.\n");
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
