//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   units/test/tstUnitSystemEnums.cc
 * \author Kelly Thompson
 * \date   Wed Oct  8 13:50:19 2003
 * \brief  Unit tests for units/UnitSystemEnum class
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"
#include <iomanip>
#include <sstream>

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_enumValues(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using std::string;

  {
    int constexpr iSenVal(3);
    constexpr std::array<double, iSenVal> adSenVal = {0.0, 1.0, 100.0};
    string const sSenVal("NA,m,cm");
    string const sSenVal2("no length unit specified,meter,centimeter");

    constexpr bool is_iSenVal = rtt_units::num_Ltype == iSenVal;
    UT_MSG(is_iSenVal, "num_Ltype has the expected value.");

    UT_MSG(soft_equiv(rtt_units::L_cf.begin(), rtt_units::L_cf.end(), adSenVal.begin(),
                      adSenVal.end()),
           "L_cf has the expected values.");

    UT_MSG(rtt_units::L_labels == sSenVal, "L_labels has the expected values.");

    UT_MSG(rtt_units::L_long_labels == sSenVal2, "L_long_labels has the expected values.");
  }

  //------------------------------------------------------------------------------------------------//
  {
    int constexpr iSenVal(3);
    constexpr std::array<double, iSenVal> adSenVal = {0.0, 1.0, 1000.0};
    string const sSenVal("NA,kg,g");
    string const sSenVal2("no mass unit specified,kilogram,gram");

    constexpr bool is_iSenVal = rtt_units::num_Mtype == iSenVal;
    UT_MSG(is_iSenVal, "num_Mtype has the expected value.");

    if (soft_equiv(rtt_units::M_cf.begin(), rtt_units::M_cf.end(), adSenVal.begin(),
                   adSenVal.end()))
      PASSMSG("M_cf has the expected values.");
    else
      FAILMSG("M_cf does not have the expected values.");

    if (rtt_units::M_labels == sSenVal)
      PASSMSG("M_labels has the expected values.");
    else
      PASSMSG("M_labels does not have the expected values.");

    if (rtt_units::M_long_labels == sSenVal2)
      PASSMSG("M_long_labels has the expected values.");
    else
      PASSMSG("M_long_labels does not have the expected values.");
  }
  //------------------------------------------------------------------------------------------------//
  {
    int constexpr iSenVal(6);
    constexpr std::array<double, iSenVal> adSenVal = {0.0, 1.0, 1000.0, 1.0e6, 1.0e8, 1.0e9};
    string const sSenVal("NA,s,ms,us,sh,ns");
    string const sSenVal2("no time unit "
                          "specified,second,milisecond,microsecond,shake,"
                          "nanosecond");

    constexpr bool is_iSenVal = rtt_units::num_ttype == iSenVal;
    UT_MSG(is_iSenVal, "num_ttype has the expected value.");

    if (soft_equiv(rtt_units::t_cf.begin(), rtt_units::t_cf.end(), adSenVal.begin(),
                   adSenVal.end()))
      PASSMSG("t_cf has the expected values.");
    else
      FAILMSG("t_cf does not have the expected values.");

    if (rtt_units::t_labels == sSenVal)
      PASSMSG("t_labels has the expected values.");
    else
      PASSMSG("t_labels does not have the expected values.");

    if (rtt_units::t_long_labels == sSenVal2)
      PASSMSG("t_long_labels has the expected values.");
    else
      PASSMSG("t_long_labels does not have the expected values.");
  }
  //------------------------------------------------------------------------------------------------//
  {
    double constexpr keV2K(1.16045193028089e7);
    int constexpr iSenVal(4);
    constexpr std::array<double, iSenVal> adSenVal = {0.0, 1.0, 1.0 / keV2K, 1.0e3 / keV2K};
    string const sSenVal("NA,K,keV,eV");

    constexpr bool is_iSenVal = rtt_units::num_Ttype == iSenVal;
    UT_MSG(is_iSenVal, "num_Ttype has the expected value.");

    if (soft_equiv(rtt_units::T_cf.begin(), rtt_units::T_cf.end(), adSenVal.begin(),
                   adSenVal.end())) {
      PASSMSG("T_cf has the expected values.");
    } else {
      std::ostringstream msg;
      msg << "T_cf does not have the expected values.";
      for (size_t i = 0; i < iSenVal; ++i)
        msg << "\n\ti = " << i << ", T_cf = " << std::setprecision(16) << rtt_units::T_cf[i]
            << " =? adSenVal = " << adSenVal[i];
      FAILMSG(msg.str());
    }

    if (rtt_units::T_labels == sSenVal)
      PASSMSG("T_labels has the expected values.");
    else
      PASSMSG("T_labels does not have the expected values.");
  }
  //------------------------------------------------------------------------------------------------//
  {
    int constexpr iSenVal(2);
    constexpr std::array<double, iSenVal> adSenVal = {0.0, 1.0};
    string const sSenVal("NA,Amp");

    constexpr bool is_iSenVal = rtt_units::num_Itype == iSenVal;
    UT_MSG(is_iSenVal, "num_Itype has the expected value.");

    if (soft_equiv(rtt_units::I_cf.begin(), rtt_units::I_cf.end(), adSenVal.begin(),
                   adSenVal.end()))
      PASSMSG("I_cf has the expected values.");
    else
      FAILMSG("I_cf does not have the expected values.");

    if (rtt_units::I_labels == sSenVal)
      PASSMSG("I_labels has the expected values.");
    else
      PASSMSG("I_labels does not have the expected values.");
  }
  //------------------------------------------------------------------------------------------------//
  {
    int constexpr iSenVal(3);
    constexpr std::array<double, iSenVal> adSenVal = {0.0, 1.0, 57.295779512896171};
    string const sSenVal("NA,rad,deg");

    constexpr bool is_iSenVal = rtt_units::num_Atype == iSenVal;
    UT_MSG(is_iSenVal, "num_Atype has the expected value.");

    if (soft_equiv(rtt_units::A_cf.begin(), rtt_units::A_cf.end(), adSenVal.begin(),
                   adSenVal.end()))
      PASSMSG("A_cf has the expected values.");
    else
      FAILMSG("A_cf does not have the expected values.");

    if (rtt_units::A_labels == sSenVal)
      PASSMSG("A_labels has the expected values.");
    else
      PASSMSG("A_labels does not have the expected values.");
  }
  //------------------------------------------------------------------------------------------------//
  {
    int constexpr iSenVal(2);
    constexpr std::array<double, iSenVal> adSenVal = {0.0, 1.0};
    string const sSenVal("NA,mol");

    constexpr bool is_iSenVal = rtt_units::num_Qtype == iSenVal;
    UT_MSG(is_iSenVal, "num_Qtype has the expected value.");

    if (soft_equiv(rtt_units::Q_cf.begin(), rtt_units::Q_cf.end(), adSenVal.begin(),
                   adSenVal.end()))
      PASSMSG("Q_cf has the expected values.");
    else
      FAILMSG("Q_cf does not have the expected values.");

    if (rtt_units::Q_labels == sSenVal)
      PASSMSG("Q_labels has the expected values.");
    else
      PASSMSG("Q_labels does not have the expected values.");
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_enumValues(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstUnitSystemEnums.cc
//------------------------------------------------------------------------------------------------//
