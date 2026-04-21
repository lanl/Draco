//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   src/units/test/tstPhysicalConstants.cc
 * \author Kelly Thompson
 * \brief  test the PhysicalConstants class
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "units/PhysicalConstants.hh"
#include "units/PhysicalConstantsSI.hh"
#include "units/UnitSystem.hh"
#include "units/UnitSystemType.hh"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <sstream>
#include <string>

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_static_access(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::EV2K;
  using std::cout;
  using std::ostringstream;
  using std::string;

  // eV2K

  double constexpr dev = 11604.51930280894;
  if (soft_equiv(dev, EV2K)) {
    PASSMSG("Found expected value for eV2K (static).\n");
  } else {
    ostringstream msg;
    msg << "Did not find expected value for eV2K (static).\n\tThe value returned was "
        << std::setprecision(16) << EV2K << " != " << std::setprecision(16) << dev << ".\n";
    FAILMSG(msg.str());
  }

  return;
}

//------------------------------------------------------------------------------------------------//

void test_ctor(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::PhysicalConstants;
  using rtt_units::UnitSystem;
  using rtt_units::UnitSystemType;
  using std::ostringstream;

  UnitSystem const us(UnitSystemType().SI());
  PhysicalConstants const pc_def;
  PhysicalConstants const pc_si(us);

  // Check equality of objects built with different constructors.

  if (soft_equiv(pc_def.h(), pc_si.h())) {
    PASSMSG("Default and SI objects have same h.\n");
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same h.\n\tdefault.h() = " << std::setprecision(16)
        << pc_def.h() << " != si.h() = " << std::setprecision(16) << pc_si.h() << ".\n";
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.R(), pc_si.R())) {
    PASSMSG("Default and SI objects have same R.\n");
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same R.\n\tdefault.R() = " << std::setprecision(16)
        << pc_def.R() << " != si.R() = " << std::setprecision(16) << pc_si.R() << ".\n";
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.k(), pc_si.k())) {
    PASSMSG("Default and SI objects have same k.\n");
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same k.\n\tdefault.k() = " << std::setprecision(16)
        << pc_def.k() << " != si.k() = " << std::setprecision(16) << pc_si.k() << ".\n";
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.e(), pc_si.e())) {
    PASSMSG("Default and SI objects have same e.\n");
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same e.\n\tdefault.e() = " << std::setprecision(16)
        << pc_def.e() << " != si.e() = " << std::setprecision(16) << pc_si.e() << ".\n";
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.c(), pc_si.c())) {
    PASSMSG("Default and SI objects have same c.\n");
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same c.\n\tdefault.c() = " << std::setprecision(16)
        << pc_def.c() << " != si.c() = " << std::setprecision(16) << pc_si.c() << ".\n";
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.sigma(), pc_si.sigma())) {
    PASSMSG("Default and SI objects have same sigma.\n");
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same sigma.\n\tdefault.sigma() = "
        << std::setprecision(16) << pc_def.sigma() << " != si.sigma() = " << std::setprecision(16)
        << pc_si.sigma() << ".\n";
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.G(), pc_si.G())) {
    PASSMSG("Default and SI objects have same G.\n");
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same G.\n\tdefault.G() = " << std::setprecision(16)
        << pc_def.G() << " != si.G() = " << std::setprecision(16) << pc_si.G() << ".\n";
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.g(), pc_si.g())) {
    PASSMSG("Default and SI objects have same g.\n");
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same g.\n\tdefault.g() = " << std::setprecision(16)
        << pc_def.g() << " != si.g() = " << std::setprecision(16) << pc_si.g() << ".\n";
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.F(), pc_si.F())) {
    PASSMSG("Default and SI objects have same F.\n");
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same F.\n\tdefault.F() = " << std::setprecision(16)
        << pc_def.F() << " != si.F() = " << std::setprecision(16) << pc_si.F() << ".\n";
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.mu0(), pc_si.mu0())) {
    PASSMSG("Default and SI objects have same mu0.\n");
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same mu0.\n\tdefault.mu0() = "
        << std::setprecision(16) << pc_def.mu0() << " != si.mu0() = " << std::setprecision(16)
        << pc_si.mu0() << ".\n";
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.epsi0(), pc_si.epsi0())) {
    PASSMSG("Default and SI objects have same epsi0.\n");
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same epsi0.\n\tdefault.epsi0() = "
        << std::setprecision(16) << pc_def.epsi0() << " != si.epsi0() = " << std::setprecision(16)
        << pc_si.epsi0() << ".\n";
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.me(), pc_si.me())) {
    PASSMSG("Default and SI objects have same me.\n");
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same me.\n\tdefault.me() = " << std::setprecision(16)
        << pc_def.me() << " != si.me() = " << std::setprecision(16) << pc_si.me() << ".\n";
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.mp(), pc_si.mp())) {
    PASSMSG("Default and SI objects have same mp.\n");
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same mp.\n\tdefault.mp() = " << std::setprecision(16)
        << pc_def.mp() << " != si.mp() = " << std::setprecision(16) << pc_si.mp() << ".\n";
    FAILMSG(msg.str());
  }

  return;
}

//------------------------------------------------------------------------------------------------//

void test_scaled_values(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::PhysicalConstants;
  using rtt_units::UnitSystem;
  using rtt_units::UnitSystemType;
  using std::ostringstream;
  using std::pow;

  // test scaled values against expected values

  UnitSystem const us(UnitSystemType().X4());
  UnitSystem const si(UnitSystemType().SI());
  PhysicalConstants pc(us);

  double dev;

  // Avogadro

  dev = 6.02214129e+23;
  if (soft_equiv(dev, pc.avogadro()) && soft_equiv(dev, pc.Na())) {
    PASSMSG("Found expected value for Avogadro's number.\n");
  } else {
    ostringstream msg;
    msg << "Did not find expected value for Avogadro's number.\n\tThe value returned was "
        << std::setprecision(16) << pc.avogadro() << " != " << std::setprecision(16) << dev
        << ".\n";
    FAILMSG(msg.str());
  }

  // Planck

  dev = 6.6260755e-34; // J * s
  dev = dev * us.e() / si.e() * us.t() / si.t();
  if (soft_equiv(pc.planck(), dev)) {
    PASSMSG("Scaled Planck constant looks correct.\n");
  } else {
    ostringstream msg;
    msg << "Scaled Planck constant is not correct.\n\tvalue =  " << std::setprecision(16)
        << pc.planck() << " != " << std::setprecision(16) << dev << ".\n";
    FAILMSG(msg.str());
  }

  // gasConstant

  dev = 8.3144621; // J/mol/K
  dev = dev * us.e() / si.e() / (us.Q() / si.Q() * us.T() / si.T());
  if (soft_equiv(pc.gasConstant(), dev, 1.0e-8)) {
    PASSMSG("Scaled Gas constant looks correct.\n");
  } else {
    ostringstream msg;
    msg << "Scaled Gas constant is not correct.\n\tvalue =  " << std::setprecision(16)
        << pc.gasConstant() << " != " << std::setprecision(16) << dev << ".\n";
    FAILMSG(msg.str());
  }

  // Boltzmann

  dev = pc.gasConstant() / pc.avogadro();
  if (soft_equiv(pc.boltzmann(), dev)) {
    PASSMSG("Scaled Boltzmann constant looks correct.\n");
  } else {
    ostringstream msg;
    msg << "Scaled Boltzmann constant is not correct.\n\tvalue =  " << std::setprecision(16)
        << pc.boltzmann() << " != " << std::setprecision(16) << dev << ".\n";
    FAILMSG(msg.str());
  }

  // Electron Charge

  dev = 1.60217733e-19; // Amp / sec (charge)
  dev = dev / (us.t() / us.T());
  if (soft_equiv(pc.electronCharge(), dev)) {
    PASSMSG("Scaled electron charge looks correct.\n");
  } else {
    ostringstream msg;
    msg << "Scaled electron charge is not correct.\n\tvalue =  " << std::setprecision(16)
        << pc.electronCharge() << " != " << std::setprecision(16) << dev << ".\n";
    FAILMSG(msg.str());
  }

  // Speed of Light

  dev = 2.99792458e8; // m/s
  dev = dev * us.v() / si.v();
  if (soft_equiv(pc.speedOfLight(), dev)) {
    PASSMSG("Scaled speed of light looks correct.\n");
  } else {
    ostringstream msg;
    msg << "Scaled speed of light is not correct.\n\tvalue =  " << std::setprecision(16)
        << pc.speedOfLight() << " != " << std::setprecision(16) << dev << ".\n";
    FAILMSG(msg.str());
  }

  // Stefan-Boltzmann constant

  dev = 5.6703726e-8; // W / m^2 / K^4
  dev = dev * us.p() / si.p() * pow(si.L() / us.L(), 2) * pow(si.T() / us.T(), 4);

  if (soft_equiv(pc.stefanBoltzmann(), dev, 1.0e-8)) {
    PASSMSG("Scaled Stefan-Boltzmann constant looks correct.\n");
  } else {
    ostringstream msg;
    msg << "Scaled Stefan-Boltzmann constant is not correct.\n\tvalue =  " << std::setprecision(16)
        << pc.stefanBoltzmann() << " != " << std::setprecision(16) << dev << ".\n";
    FAILMSG(msg.str());
  }

  // Gravitational constant

  dev = 6.67259e-11; // N / m^2 / kg^2
  dev = dev * us.f() / si.f() / pow(us.L() / si.L(), 2) / pow(us.M() / si.M(), 2);
  if (soft_equiv(pc.gravitationalConstant(), dev)) {
    PASSMSG("Scaled Gravitational constant looks correct.\n");
  } else {
    ostringstream msg;
    msg << "Scaled Gravitational constant is not correct.\n"
        << "\tvalue =  " << std::setprecision(16) << pc.gravitationalConstant()
        << " != " << std::setprecision(16) << dev << ".\n";
    FAILMSG(msg.str());
  }

  // Acceleration from Gravity

  dev = 9.80665; // m/s^2
  dev = dev * us.a() / si.a();
  if (soft_equiv(pc.accelerationFromGravity(), dev)) {
    PASSMSG("Scaled Acceleration from Gravity looks correct.\n");
  } else {
    ostringstream msg;
    msg << "Scaled Acceleration from Gravity is not correct.\n\tvalue =  " << std::setprecision(16)
        << pc.accelerationFromGravity() << " != " << std::setprecision(16) << dev << ".\n";
    FAILMSG(msg.str());
  }

  // Faraday constant

  dev = pc.avogadro() * pc.electronCharge();
  if (soft_equiv(pc.faradayConstant(), dev)) {
    PASSMSG("Scaled Faraday constant looks correct.\n");
  } else {
    ostringstream msg;
    msg << "Scaled Faraday constant is not correct.\n\tvalue =  " << std::setprecision(16)
        << pc.faradayConstant() << " != " << std::setprecision(16) << dev << ".\n";
    FAILMSG(msg.str());
  }

  // Permeability of vacuum (free space)

  dev = 4.0 * std::numbers::pi * 1.0e-7; // H / m
  dev = dev * si.L() / us.L();
  if (soft_equiv(pc.permeabilityOfVacuum(), dev)) {
    PASSMSG("Scaled Permeability of vacuum (free space) looks correct.\n");
  } else {
    ostringstream msg;
    msg << "Scaled Permeability of vacuum (free space) is not correct.\n\tvalue =  "
        << std::setprecision(16) << pc.permeabilityOfVacuum() << " != " << std::setprecision(16)
        << dev << ".\n";
    FAILMSG(msg.str());
  }

  // PERMITTIVITY OF FREE SPACE (F/M)

  dev = 1.0 / pc.permeabilityOfVacuum() / pow(pc.speedOfLight(), 2);
  if (soft_equiv(pc.permittivityOfFreeSpace(), dev)) {
    PASSMSG("Scaled PERMITTIVITY OF FREE SPACE (F/M) looks correct.\n");
  } else {
    ostringstream msg;
    msg << "Scaled PERMITTIVITY OF FREE SPACE (F/M) is not correct.\n\tvalue =  "
        << std::setprecision(16) << pc.permittivityOfFreeSpace() << " != " << std::setprecision(16)
        << dev << ".\n";
    FAILMSG(msg.str());
  }

  // Classical electron radius
  dev = 2.8179403227e-15; // m (NIST value)
  dev = dev * us.L() / si.L();
  if (soft_equiv(pc.classicalElectronRadius(), dev, 2e-9)) {
    PASSMSG("Scaled classical electron radius looks correct.\n");
  } else {
    ostringstream msg;
    msg << "Scaled classical electron radius is not correct.\n\tvalue =  " << std::setprecision(16)
        << pc.classicalElectronRadius() << " != " << std::setprecision(16) << dev << ".\n";
    FAILMSG(msg.str());
  }
  // test alias
  if (soft_equiv(pc.re(), dev, 2e-9)) {
    PASSMSG("Scaled classical electron radius looks correct.\n");
  } else {
    ostringstream msg;
    msg << "Scaled classical electron radius is not correct.\n\tvalue =  " << std::setprecision(16)
        << pc.re() << " != " << std::setprecision(16) << dev << ".\n";
    FAILMSG(msg.str());
  }

  // Electron mass

  dev = 9.1093897e-31; // kg
  dev = dev * us.M() / si.M();
  if (soft_equiv(pc.electronMass(), dev)) {
    PASSMSG("Scaled Electron mass looks correct.\n");
  } else {
    ostringstream msg;
    msg << "Scaled Electron mass is not correct.\n\tvalue =  " << std::setprecision(16)
        << pc.electronMass() << " != " << std::setprecision(16) << dev << ".\n";
    FAILMSG(msg.str());
  }

  // Proton mass

  dev = 1.6726231e-27; // kg
  dev = dev * us.M() / si.M();
  if (soft_equiv(pc.protonMass(), dev)) {
    PASSMSG("Scaled Proton mass looks correct.\n");
  } else {
    ostringstream msg;
    msg << "Scaled Proton mass is not correct.\n\tvalue =  " << std::setprecision(16)
        << pc.protonMass() << " != " << std::setprecision(16) << dev << ".\n";
    FAILMSG(msg.str());
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_static_access(ut);
    test_ctor(ut);
    test_scaled_values(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstPhysicalConstants.cc
//------------------------------------------------------------------------------------------------//
