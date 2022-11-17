//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   units/test/tstUnitSystemType.cc
 * \author Kelly Thompson
 * \date   Wed Oct  8 13:50:19 2003
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "units/UnitSystemType.hh"
#include <sstream>

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_copy_ctor(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::FundUnit;
  using rtt_units::Ltype;
  using rtt_units::UnitSystemType;
  using std::endl;
  using std::ostringstream;

  //Let's use the default SI units to test

  UnitSystemType si(UnitSystemType().SI());
  UnitSystemType si2(si);

  //Only testing one element
  if (si.L().enumVal() == si2.L().enumVal()) {
    ostringstream msg;
    msg << "Copied UnitSystemType container has "
        << "orig.L() == L_copy.L()." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Copied UnitSystemType container does not have "
        << "orig.L() == L_copy.L()." << endl;
    FAILMSG(msg.str());
  }

  if (&si != &si2) {
    ostringstream msg;
    msg << "Orig and copied object are seperate. " << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Orig and copied object point to same memory. " << endl;
    FAILMSG(msg.str());
  }
}

void test_move_ctor(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::FundUnit;
  using rtt_units::Ltype;
  using rtt_units::UnitSystemType;
  using std::endl;
  using std::ostringstream;

  //Let's use the default SI units to test

  UnitSystemType si(UnitSystemType().SI());

  auto Ltest = si.L().enumVal();

  UnitSystemType si2(std::move(si));

  //Only testing one element
  if (si2.L().enumVal() == Ltest) {
    ostringstream msg;
    msg << "Moved UnitSystemType container has "
        << "orig.L() == L_copy.L()." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Moved UnitSystemType container does not have "
        << "orig.L() == L_copy.L()." << endl;
    FAILMSG(msg.str());
  }
}

void test_default_ctor(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::FundUnit;
  using rtt_units::Ltype;
  using rtt_units::UnitSystemType;
  using std::endl;
  using std::ostringstream;

  // Create a container with no data...
  UnitSystemType ust;

  // ensure that the container is empty but valid...
  if (ust.L().enumVal() == 0) {
    ostringstream msg;
    msg << "Empty UnitSystemType container has "
        << "L().enumVal() == L_null." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Empty UnitSystemType container does not have "
        << "L().enumVal() == L_null." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(ust.L().cf(), 0.0)) {
    PASSMSG("Empty UnitSystemType container has L().cf() == 0.0).");
  } else {
    ostringstream msg;
    msg << "Empty UnitSystemType container does not have"
        << "L().cf() == 0.0)." << endl;
    FAILMSG(msg.str());
  }
  if (ust.L().label() == "NA") {
    PASSMSG("Empty UnitSystemType container has L().label() == NA).");
  } else {
    ostringstream msg;
    msg << "Empty UnitSystemType container does not have "
        << "L().label() == NA)." << endl;
    FAILMSG(msg.str());
  }

  // Establish a fundamental unit for the length, mass and Temperature portions of UnitSystemType...
  ust.L(rtt_units::L_m).M(rtt_units::M_kg).T(rtt_units::T_K);

  // Check that the object was updated...
  if (ust.L().enumVal() == rtt_units::L_m) {
    PASSMSG("UnitSystemType container has L().enumVal() == L_m).");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "L().enumVal() == L_m)." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(ust.L().cf(), 1.0))
    PASSMSG("UnitSystemType container has L().cf() == 1.0).");
  else
    FAILMSG("UnitSystemType container does not have L().cf() == 1.0).");
  if (ust.L().label() == "m") {
    PASSMSG("UnitSystemType container has L().label() == \"m\").");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "L().label() == \"m\")." << endl;
    FAILMSG(msg.str());
  }

  // Ensure that other empty parts of ust are valid...
  if (ust.Q().enumVal() == rtt_units::Q_null) {
    PASSMSG("UnitSystemType container has Q().enumVal() == Q_null).");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "Q().enumVal() == Q_null)." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(ust.Q().cf(), 0.0))
    PASSMSG("UnitSystemType container has Q().cf() == 0.0).");
  else
    FAILMSG("UnitSystemType container does not have Q().cf() == 0.0).");
  if (ust.Q().label() == "NA") {
    PASSMSG("UnitSystemType container has Q().label() == \"NA\").");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "Q().label() == \"NA\")." << endl;
    FAILMSG(msg.str());
  }
  //---------------------------------------------------------------------//
  if (ust.A().enumVal() == rtt_units::A_null) {
    PASSMSG("UnitSystemType container has A().enumVal() == A_null).");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "A().enumVal() == A_null)." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(ust.A().cf(), 0.0))
    PASSMSG("UnitSystemType container has A().cf() == 0.0).");
  else
    FAILMSG("UnitSystemType container does not have A().cf() == 0.0).");
  if (ust.A().label() == "NA") {
    PASSMSG("UnitSystemType container has A().label() == \"NA\").");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "A().label() == \"NA\"." << endl;
    FAILMSG(msg.str());
  }
  //---------------------------------------------------------------------//

  // Add more FundUnits to ust using a cascading set of manipulators...
  ust.Q(rtt_units::Q_mol)
      .A(rtt_units::A_rad)
      .I(rtt_units::I_amp)
      .t(rtt_units::t_sh, rtt_units::t_cf.data(), rtt_units::t_long_labels);

  //---------------------------------------------------------------------//

  // Ensure that the the updated values are in place...
  if (ust.Q().enumVal() == rtt_units::Q_mol) {
    PASSMSG("UnitSystemType container has Q().enumVal() == Q_mol).");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "Q().enumVal() == Q_mol.  Instead the valued retuned was: " << ust.Q().enumVal() << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(ust.Q().cf(), 1.0)) {
    PASSMSG("UnitSystemType container has Q().cf() == 1.0).");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "Q().cf() == 1.0.  Instead the value returned was: " << ust.Q().cf() << endl;
    FAILMSG(msg.str());
  }
  if (ust.Q().label() == "mol") {
    PASSMSG("UnitSystemType container has Q().label() == \"mol\").");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << R"(Q().label() == "mol".  Instead the value returned was: ")" << ust.Q().label() << "\""
        << endl;
    // " fix emacs font-lock issue with odd number of quotes.
    FAILMSG(msg.str());
  }
  //---------------------------------------------------------------------//
  if (ust.A().enumVal() == rtt_units::A_rad) {
    PASSMSG("UnitSystemType container has A().enumVal() == A_rad).");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "A().enumVal() == A_rad." << endl
        << "\tInstead the returned value was: " << ust.A().enumVal() << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(ust.A().cf(), 1.0)) {
    PASSMSG("UnitSystemType container has A().cf() == 1.0).");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "A().cf() == 1.0." << endl
        << "\tInstead the returned value was:" << ust.A().cf() << endl;
    FAILMSG(msg.str());
  }
  if (ust.A().label() == "rad") {
    PASSMSG("UnitSystemType container has A().label() == \"rad\").");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "A().label() == \"rad\"." << endl
        << "\tInstead the returned value was: \"" << ust.A().label() << "\"" << endl;
    FAILMSG(msg.str());
  }
  //---------------------------------------------------------------------//
  if (ust.I().enumVal() == rtt_units::I_amp) {
    PASSMSG("UnitSystemType container has I().enumVal() == I_amp).");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "I().enumVal() == I_amp." << endl
        << "\tInstead the returned value was: " << ust.I().enumVal() << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(ust.I().cf(), 1.0)) {
    PASSMSG("UnitSystemType container has I().cf() == 1.0).");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "I().cf() == 1.0." << endl
        << "\tInstead the returned value was:" << ust.I().cf() << endl;
    FAILMSG(msg.str());
  }
  if (ust.I().label() == "Amp") {
    PASSMSG("UnitSystemType container has I().label() == \"Amp\").");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "I().label() == \"Amp\"." << endl
        << "\tInstead the returned value was: \"" << ust.I().label() << "\"" << endl;
    FAILMSG(msg.str());
  }
  //---------------------------------------------------------------------//
  if (ust.t().enumVal() == rtt_units::t_sh) {
    PASSMSG("UnitSystemType container has t().enumVal() == t_sh.");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "t().enumVal() == t_rad." << endl
        << "\tInstead the returned value was: " << ust.t().enumVal() << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(ust.t().cf(), 1.0e8)) {
    PASSMSG("UnitSystemType container has t().cf() == 1.0e8.");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "t().cf() == 1.0e8." << endl
        << "\tInstead the returned value was:" << ust.t().cf() << endl;
    FAILMSG(msg.str());
  }
  if (ust.t().label() == "shake") {
    PASSMSG("UnitSystemType container has t().label() == \"shake\".");
  } else {
    ostringstream msg;
    msg << "UnitSystemType container does not have "
        << "t().label() == \"shake\"." << endl
        << "\tInstead the returned value was: \"" << ust.t().label() << "\"" << endl;
    FAILMSG(msg.str());
  }

  return;
}

//------------------------------------------------------------------------//
void test_qualified_ctor(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::UnitSystemType;
  using std::endl;
  using std::ostringstream;

  UnitSystemType ust(rtt_units::L_cm, rtt_units::M_g, rtt_units::t_us, rtt_units::T_K,
                     rtt_units::I_amp, rtt_units::A_rad, rtt_units::Q_mol);

  //----------------------------------------//

  if (ust.M().enumVal() == rtt_units::M_g) {
    ostringstream msg;
    msg << "UnitSystemType has "
        << "ust.M().enumVal() == rtt_units::M_g" << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "UnitSystemType does not have "
        << "ust.M().enumVal() == rtt_units::M_g" << endl
        << "The value reported was " << ust.M().enumVal() << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(ust.M().cf(), 1000.0)) {
    ostringstream msg;
    msg << "UnitSystemType has "
        << "ust.M().cf() == 1000.0" << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "UnitSystemType does not have "
        << "ust.M().cf() == 1000.0" << endl
        << "The value reported was " << ust.M().cf() << endl;
    FAILMSG(msg.str());
  }
  if (ust.M().label() == "g") {
    ostringstream msg;
    msg << "UnitSystemType has "
        << "ust.M().label() == \"g\"" << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "UnitSystemType did not return the expected value " << endl
        << "( ust.M().label() == \"g\" ).  Instead it returned the " << endl
        << "value: \"" << ust.M().label() << "\"" << endl;
    FAILMSG(msg.str());
  }

  //----------------------------------------//

  if (ust.t().enumVal() == rtt_units::t_us) {
    ostringstream msg;
    msg << "UnitSystemType has "
        << "ust.t().enumVal() == rtt_units::t_us" << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "UnitSystemType does not have "
        << "ust.t().enumVal() == rtt_units::t_us" << endl
        << "\tThe value reported was " << ust.t().enumVal() << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(ust.t().cf(), 1.0e6)) {
    ostringstream msg;
    msg << "UnitSystemType has "
        << "ust.t().cf() == 1.0e6" << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "UnitSystemType does not have "
        << "ust.t().cf() == 1.0e6" << endl
        << "\tThe value reported was " << ust.t().cf() << endl;
    FAILMSG(msg.str());
  }
  if (ust.t().label() == "us") {
    ostringstream msg;
    msg << "UnitSystemType has "
        << "ust.t().label() == \"us\"" << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "UnitSystemType did not return the expected value " << endl
        << "\t( ust.t().label() == \"us\" ).  Instead it returned the " << endl
        << "\tvalue: \"" << ust.t().label() << "\"" << endl;
    FAILMSG(msg.str());
  }

  //----------------------------------------//

  if (ust.T().enumVal() == rtt_units::T_K) {
    ostringstream msg;
    msg << "UnitSystemType has "
        << "ust.T().enumVal() == rtt_units::t_us" << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "UnitSystemType does not have "
        << "ust.T().enumVal() == rtt_units::T_K" << endl
        << "\tThe value reported was " << ust.T().enumVal() << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(ust.T().cf(), 1.0)) {
    ostringstream msg;
    msg << "UnitSystemType has "
        << "ust.T().cf() == 1.0" << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "UnitSystemType does not have "
        << "ust.T().cf() == 1.0" << endl
        << "\tThe value reported was " << ust.T().cf() << endl;
    FAILMSG(msg.str());
  }
  if (ust.T().label() == "K") {
    ostringstream msg;
    msg << "UnitSystemType has "
        << "ust.T().label() == \"K\"" << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "UnitSystemType did not return the expected value " << endl
        << "\t( ust.T().label() == \"K\" ).  Instead it returned the " << endl
        << "\tvalue: \"" << ust.T().label() << "\"" << endl;
    FAILMSG(msg.str());
  }

  return;
}

//------------------------------------------------------------------------//
/*!
 * \brief Test special unit system types.
 */
void test_sust_ctor(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::UnitSystemType;
  using std::endl;
  using std::ostringstream;

  // Create the two "speical" unit system types.

  UnitSystemType si(UnitSystemType().SI());
  UnitSystemType x4(UnitSystemType().X4());
  UnitSystemType cgs(UnitSystemType().CGS());

  // Check Length data...
  if (si.L().enumVal() == rtt_units::L_m) {
    PASSMSG("si has L().enumVal() == L_m).");
  } else {
    ostringstream msg;
    msg << "si does not have "
        << "L().enumVal() == L_m)." << endl;
    FAILMSG(msg.str());
  }
  if (x4.L().enumVal() == rtt_units::L_cm) {
    PASSMSG("x4 has L().enumVal() == L_cm).");
  } else {
    ostringstream msg;
    msg << "x4 does not have "
        << "L().enumVal() == L_cm." << endl;
    FAILMSG(msg.str());
  }
  if (cgs.L().enumVal() == rtt_units::L_cm) {
    PASSMSG("cgs has L().enumVal() == L_cm).");
  } else {
    ostringstream msg;
    msg << "cgs does not have "
        << "L().enumVal() == L_cm." << endl;
    FAILMSG(msg.str());
  }

  if (soft_equiv(si.L().cf(), 1.0)) {
    PASSMSG("si container has L().cf() == 1.0).");
  } else {
    FAILMSG("si container does not have L().cf() == 1.0).");
  }
  if (si.L().label() == "m") {
    PASSMSG("si container has L().label() == \"m\").");
  } else {
    ostringstream msg;
    msg << "si container does not have "
        << "L().label() == \"m\")." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(x4.L().cf(), 100.0)) {
    PASSMSG("x4 units container has L().cf() == 100.0.");
  } else {
    FAILMSG("x4 units container does not have L().cf() == 100.0.");
  }
  if (x4.L().label() == "cm") {
    PASSMSG("x4 units container has L().label() == \"cm\").");
  } else {
    ostringstream msg;
    msg << "x4 units container does not have "
        << "L().label() == \"cm\")." << endl;
    FAILMSG(msg.str());
  }

  // Check Mass data...
  // ----------------------------------------

  if (si.M().enumVal() == rtt_units::M_kg) {
    PASSMSG("si has M().enumVal() == L_m).");
  } else {
    ostringstream msg;
    msg << "si does not have "
        << "M().enumVal() == L_m)." << endl;
    FAILMSG(msg.str());
  }
  if (x4.M().enumVal() == rtt_units::M_g) {
    PASSMSG("x4 has M().enumVal() == M_g.");
  } else {
    ostringstream msg;
    msg << "x4 does not have "
        << "M().enumVal() == M_g." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(si.M().cf(), 1.0)) {
    PASSMSG("si container has M().cf() == 1.0.");
  } else {
    FAILMSG("si container does not have M().cf() == 1.0.");
  }
  if (si.M().label() == "kg") {
    PASSMSG("si container has M().label() == \"kg\").");
  } else {
    ostringstream msg;
    msg << "si container does not have "
        << "M().label() == \"kg\")." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(x4.M().cf(), 1000.0)) {
    PASSMSG("x4 units container has M().cf() == 1000.0.");
  } else {
    FAILMSG("x4 units container does not have M().cf() == 1000.0.");
  }
  if (x4.M().label() == "g") {
    PASSMSG("x4 units container has M().label() == \"g\".");
  } else {
    ostringstream msg;
    msg << "x4 units container does not have "
        << "M().label() == \"g\"." << endl;
    FAILMSG(msg.str());
  }

  // Check time data...
  // ----------------------------------------

  if (si.t().enumVal() == rtt_units::t_s) {
    PASSMSG("si has t().enumVal() == t_s.");
  } else {
    ostringstream msg;
    msg << "si does not have "
        << "t().enumVal() == t_s." << endl;
    FAILMSG(msg.str());
  }
  if (x4.t().enumVal() == rtt_units::t_sh) {
    PASSMSG("x4 has t().enumVal() == t_sh.");
  } else {
    ostringstream msg;
    msg << "x4 does not have "
        << "t().enumVal() == t_sh." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(si.t().cf(), 1.0)) {
    PASSMSG("si container has t().cf() == 1.0.");
  } else {
    FAILMSG("si container does not have t().cf() == 1.0.");
  }
  if (si.t().label() == "s") {
    PASSMSG("si container has t().label() == \"s\").");
  } else {
    ostringstream msg;
    msg << "si container does not have "
        << "t().label() == \"s\"." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(x4.t().cf(), 1.0e8)) {
    PASSMSG("x4 units container has t().cf() == 1.0e8.");
  } else {
    FAILMSG("x4 units container does not have t().cf() == 1.0e8.");
  }
  if (x4.t().label() == "sh") {
    PASSMSG("x4 units container has t().label() == \"sh\".");
  } else {
    ostringstream msg;
    msg << "x4 units container does not have "
        << "t().label() == \"sh\"." << endl;
    FAILMSG(msg.str());
  }

  // Check Temperature data...
  // ----------------------------------------

  if (si.T().enumVal() == rtt_units::T_K) {
    PASSMSG("si has T().enumVal() == T_K).");
  } else {
    ostringstream msg;
    msg << "si does not have "
        << "T().enumVal() == T_K)." << endl;
    FAILMSG(msg.str());
  }
  if (x4.T().enumVal() == rtt_units::T_keV) {
    PASSMSG("x4 has T().enumVal() == T_kev.");
  } else {
    ostringstream msg;
    msg << "x4 does not have "
        << "T().enumVal() == T_kev." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(si.T().cf(), 1.0)) {
    PASSMSG("si container has T().cf() == 1.0).");
  } else {
    FAILMSG("si container does not have T().cf() == 1.0).");
  }
  if (si.T().label() == "K") {
    PASSMSG("si container has T().label() == \"K\".");
  } else {
    ostringstream msg;
    msg << "si container does not have "
        << "T().label() == \"K\"." << endl;
    FAILMSG(msg.str());
  }

  double const keV2K(1.16045193028089e7);
  if (soft_equiv(x4.T().cf(), 1.0 / keV2K)) {
    std::ostringstream msg;
    msg << "x4 units container has T().cf() == " << x4.T().cf();
    PASSMSG(msg.str());
  } else {
    std::ostringstream msg;
    msg << "x4 units container does not have T().cf() == " << 1.0 / keV2K
        << "\n\t found T().cf() = " << x4.T().cf();
    FAILMSG(msg.str());
  }
  if (x4.T().label() == "keV") {
    PASSMSG("x4 units container has T().label() == \"keV\").");
  } else {
    ostringstream msg;
    msg << "x4 units container does not have "
        << "T().label() == \"keV\")." << endl;
    FAILMSG(msg.str());
  }

  // Check Angle data...
  // ----------------------------------------

  if (si.A().enumVal() == rtt_units::A_rad) {
    PASSMSG("si has A().enumVal() == A_rad.");
  } else {
    ostringstream msg;
    msg << "si does not have "
        << "A().enumVal() == A_rad." << endl;
    FAILMSG(msg.str());
  }
  if (x4.A().enumVal() == rtt_units::A_rad) {
    PASSMSG("x4 has A().enumVal() == A_rad.");
  } else {
    ostringstream msg;
    msg << "x4 does not have "
        << "A().enumVal() == A_rad." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(si.A().cf(), 1.0)) {
    PASSMSG("si container has A().cf() == 1.0.");
  } else {
    FAILMSG("si container does not have A().cf() == 1.0.");
  }
  if (si.A().label() == "rad") {
    PASSMSG("si container has A().label() == \"rad\".");
  } else {
    ostringstream msg;
    msg << "si container does not have "
        << "A().label() == \"rad\"." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(x4.A().cf(), 1.0)) {
    PASSMSG("x4 units container has A().cf() == 1.0.");
  } else {
    FAILMSG("x4 units container does not have A().cf() == 1.0.");
  }
  if (x4.A().label() == "rad") {
    PASSMSG("x4 units container has A().label() == \"rad\").");
  } else {
    ostringstream msg;
    msg << "x4 units container does not have "
        << "A().label() == \"rad\")." << endl;
    FAILMSG(msg.str());
  }

  // Check Current data...
  // ----------------------------------------

  if (si.I().enumVal() == rtt_units::I_amp) {
    PASSMSG("si has I().enumVal() == I_amp).");
  } else {
    ostringstream msg;
    msg << "si does not have "
        << "I().enumVal() == I_amp)." << endl;
    FAILMSG(msg.str());
  }
  if (x4.I().enumVal() == rtt_units::I_amp) {
    PASSMSG("x4 has I().enumVal() == I_amp.");
  } else {
    ostringstream msg;
    msg << "x4 does not have "
        << "I().enumVal() == I_amp." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(si.I().cf(), 1.0)) {
    PASSMSG("si container has I().cf() == 1.0.");
  } else {
    FAILMSG("si container does not have I().cf() == 1.0.");
  }
  if (si.I().label() == "Amp") {
    PASSMSG("si container has I().label() == \"Amp\".");
  } else {
    ostringstream msg;
    msg << "si container does not have "
        << "I().label() == \"Amp\"." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(x4.I().cf(), 1.0)) {
    PASSMSG("x4 units container has I().cf() == 1.0");
  } else {
    FAILMSG("x4 units container does not have I().cf() == 1.0.");
  }
  if (x4.I().label() == "Amp") {
    PASSMSG("x4 units container has I().label() == \"Amp\").");
  } else {
    ostringstream msg;
    msg << "x4 units container does not have "
        << "I().label() == \"Amp\")." << endl;
    FAILMSG(msg.str());
  }

  // Check Quantity data...
  // ----------------------------------------

  if (si.Q().enumVal() == rtt_units::Q_mol) {
    PASSMSG("si has Q().enumVal() == Q_mol.");
  } else {
    ostringstream msg;
    msg << "si does not have "
        << "Q().enumVal() == Q_mol." << endl;
    FAILMSG(msg.str());
  }
  if (x4.Q().enumVal() == rtt_units::Q_mol) {
    PASSMSG("x4 has Q().enumVal() == Q_mol.");
  } else {
    ostringstream msg;
    msg << "x4 does not have "
        << "Q().enumVal() == Q_mol." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(si.Q().cf(), 1.0)) {
    PASSMSG("si container has Q().cf() == 1.0.");
  } else {
    FAILMSG("si container does not have Q().cf() == 1.0.");
  }
  if (si.Q().label() == "mol") {
    PASSMSG("si container has Q().label() == \"mol\".");
  } else {
    ostringstream msg;
    msg << "si container does not have "
        << "Q().label() == \"mol\"." << endl;
    FAILMSG(msg.str());
  }
  if (soft_equiv(x4.Q().cf(), 1.0)) {
    PASSMSG("x4 units container has Q().cf() == 1.0.");
  } else {
    FAILMSG("x4 units container does not have Q().cf() == 1.0.");
  }
  if (x4.Q().label() == "mol") {
    PASSMSG("x4 units container has Q().label() == \"mol\").");
  } else {
    ostringstream msg;
    msg << "x4 units container does not have "
        << "Q().label() == \"mol\")." << endl;
    FAILMSG(msg.str());
  }

  return;
}

//------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_default_ctor(ut);
    test_qualified_ctor(ut);
    test_sust_ctor(ut); // specific UnitSystemType ctors
    test_copy_ctor(ut);
    test_move_ctor(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstUnitSystemTypes.cc
//------------------------------------------------------------------------------------------------//
