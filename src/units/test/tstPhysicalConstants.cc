//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   src/units/test/tstPhysicalConstants.cc
 * \author Kelly Thompson
 * \date   Mon Nov  3 22:35:14 2003
 * \brief  test the PhysicalConstants class
 * \note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"
#include "units/PhysicalConstantsSI.hh"
#include <iomanip>
#include <sstream>

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_static_access(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::EV2K;
  using rtt_units::PI;
  using std::cout;
  using std::endl;
  using std::ostringstream;
  using std::string;

  // PI

  double dev(3.14159265358979324);
  if (soft_equiv(dev, PI)) {
    ostringstream msg;
    msg << "Found expected value for PI (static)." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Did not find expected value for PI (static)." << endl
        << "\tThe value returned was " << std::setprecision(16) << PI
        << " != " << std::setprecision(16) << dev << "." << endl;
    FAILMSG(msg.str());
  }

  // eV2K

  dev = 11604.51930280894;
  if (soft_equiv(dev, EV2K)) {
    ostringstream msg;
    msg << "Found expected value for eV2K (static)." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Did not find expected value for eV2K (static)." << endl
        << "\tThe value returned was " << std::setprecision(16) << EV2K
        << " != " << std::setprecision(16) << dev << "." << endl;
    FAILMSG(msg.str());
  }

  return;
}

//------------------------------------------------------------------------------------------------//

void test_ctor(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::PhysicalConstants;
  using rtt_units::PI;
  using rtt_units::UnitSystem;
  using rtt_units::UnitSystemType;
  using std::endl;
  using std::ostringstream;

  UnitSystem us(UnitSystemType().SI());
  PhysicalConstants pc_def;
  PhysicalConstants pc_si(us);

  // Check equality of objects built with different constructors.

  if (soft_equiv(pc_def.h(), pc_si.h())) {
    ostringstream msg;
    msg << "Default and SI objects have same h." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same h." << endl
        << "\tdefault.h() = " << std::setprecision(16) << pc_def.h()
        << " != si.h() = " << std::setprecision(16) << pc_si.h() << "." << endl;
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.R(), pc_si.R())) {
    ostringstream msg;
    msg << "Default and SI objects have same R." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same R." << endl
        << "\tdefault.R() = " << std::setprecision(16) << pc_def.R()
        << " != si.R() = " << std::setprecision(16) << pc_si.R() << "." << endl;
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.k(), pc_si.k())) {
    ostringstream msg;
    msg << "Default and SI objects have same k." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same k." << endl
        << "\tdefault.k() = " << std::setprecision(16) << pc_def.k()
        << " != si.k() = " << std::setprecision(16) << pc_si.k() << "." << endl;
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.e(), pc_si.e())) {
    ostringstream msg;
    msg << "Default and SI objects have same e." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same e." << endl
        << "\tdefault.e() = " << std::setprecision(16) << pc_def.e()
        << " != si.e() = " << std::setprecision(16) << pc_si.e() << "." << endl;
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.c(), pc_si.c())) {
    ostringstream msg;
    msg << "Default and SI objects have same c." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same c." << endl
        << "\tdefault.c() = " << std::setprecision(16) << pc_def.c()
        << " != si.c() = " << std::setprecision(16) << pc_si.c() << "." << endl;
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.sigma(), pc_si.sigma())) {
    ostringstream msg;
    msg << "Default and SI objects have same sigma." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same sigma." << endl
        << "\tdefault.sigma() = " << std::setprecision(16) << pc_def.sigma()
        << " != si.sigma() = " << std::setprecision(16) << pc_si.sigma() << "." << endl;
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.G(), pc_si.G())) {
    ostringstream msg;
    msg << "Default and SI objects have same G." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same G." << endl
        << "\tdefault.G() = " << std::setprecision(16) << pc_def.G()
        << " != si.G() = " << std::setprecision(16) << pc_si.G() << "." << endl;
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.g(), pc_si.g())) {
    ostringstream msg;
    msg << "Default and SI objects have same g." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same g." << endl
        << "\tdefault.g() = " << std::setprecision(16) << pc_def.g()
        << " != si.g() = " << std::setprecision(16) << pc_si.g() << "." << endl;
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.F(), pc_si.F())) {
    ostringstream msg;
    msg << "Default and SI objects have same F." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same F." << endl
        << "\tdefault.F() = " << std::setprecision(16) << pc_def.F()
        << " != si.F() = " << std::setprecision(16) << pc_si.F() << "." << endl;
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.mu0(), pc_si.mu0())) {
    ostringstream msg;
    msg << "Default and SI objects have same mu0." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same mu0." << endl
        << "\tdefault.mu0() = " << std::setprecision(16) << pc_def.mu0()
        << " != si.mu0() = " << std::setprecision(16) << pc_si.mu0() << "." << endl;
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.epsi0(), pc_si.epsi0())) {
    ostringstream msg;
    msg << "Default and SI objects have same epsi0." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same epsi0." << endl
        << "\tdefault.epsi0() = " << std::setprecision(16) << pc_def.epsi0()
        << " != si.epsi0() = " << std::setprecision(16) << pc_si.epsi0() << "." << endl;
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.me(), pc_si.me())) {
    ostringstream msg;
    msg << "Default and SI objects have same me." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same me." << endl
        << "\tdefault.me() = " << std::setprecision(16) << pc_def.me()
        << " != si.me() = " << std::setprecision(16) << pc_si.me() << "." << endl;
    FAILMSG(msg.str());
  }

  if (soft_equiv(pc_def.mp(), pc_si.mp())) {
    ostringstream msg;
    msg << "Default and SI objects have same mp." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Default and SI objects do not have same mp." << endl
        << "\tdefault.mp() = " << std::setprecision(16) << pc_def.mp()
        << " != si.mp() = " << std::setprecision(16) << pc_si.mp() << "." << endl;
    FAILMSG(msg.str());
  }

  return;
}

//------------------------------------------------------------------------------------------------//

void test_scaled_values(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::PhysicalConstants;
  using rtt_units::PI;
  using rtt_units::UnitSystem;
  using rtt_units::UnitSystemType;
  using std::endl;
  using std::ostringstream;
  using std::pow;

  // test scaled values against expected values

  UnitSystem us(UnitSystemType().X4());
  UnitSystem si(UnitSystemType().SI());
  PhysicalConstants pc(us);

  double dev;

  // Avagadro

  dev = 6.02214129e+23;
  if (soft_equiv(dev, pc.avogadro()) && soft_equiv(dev, pc.Na())) {
    ostringstream msg;
    msg << "Found expected value for Avogadro's number." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Did not find expected value for Avogadro's number." << endl
        << "\tThe value returned was " << std::setprecision(16) << pc.avogadro()
        << " != " << std::setprecision(16) << dev << "." << endl;
    FAILMSG(msg.str());
  }

  // Planck

  dev = 6.6260755e-34; // J * s
  dev = dev * us.e() / si.e() * us.t() / si.t();
  if (soft_equiv(pc.planck(), dev)) {
    ostringstream msg;
    msg << "Scaled Planck constant looks correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Scaled Planck constant is not correct." << endl
        << "\tvalue =  " << std::setprecision(16) << pc.planck() << " != " << std::setprecision(16)
        << dev << "." << endl;
    FAILMSG(msg.str());
  }

  // gasConstant

  dev = 8.3144621; // J/mol/K
  dev = dev * us.e() / si.e() / (us.Q() / si.Q() * us.T() / si.T());
  if (soft_equiv(pc.gasConstant(), dev, 1.0e-8)) {
    ostringstream msg;
    msg << "Scaled Gas constant looks correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Scaled Gas constant is not correct." << endl
        << "\tvalue =  " << std::setprecision(16) << pc.gasConstant()
        << " != " << std::setprecision(16) << dev << "." << endl;
    FAILMSG(msg.str());
  }

  // Boltzmann

  dev = pc.gasConstant() / pc.avogadro();
  if (soft_equiv(pc.boltzmann(), dev)) {
    ostringstream msg;
    msg << "Scaled Boltzmann constant looks correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Scaled Boltzmann constant is not correct." << endl
        << "\tvalue =  " << std::setprecision(16) << pc.boltzmann()
        << " != " << std::setprecision(16) << dev << "." << endl;
    FAILMSG(msg.str());
  }

  // Electron Charge

  dev = 1.60217733e-19; // Amp / sec (charge)
  dev = dev / (us.t() / us.T());
  if (soft_equiv(pc.electronCharge(), dev)) {
    ostringstream msg;
    msg << "Scaled electron charge looks correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Scaled electron charge is not correct." << endl
        << "\tvalue =  " << std::setprecision(16) << pc.electronCharge()
        << " != " << std::setprecision(16) << dev << "." << endl;
    FAILMSG(msg.str());
  }

  // Speed of Light

  dev = 2.99792458e8; // m/s
  dev = dev * us.v() / si.v();
  if (soft_equiv(pc.speedOfLight(), dev)) {
    ostringstream msg;
    msg << "Scaled speed of light looks correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Scaled speed of light is not correct." << endl
        << "\tvalue =  " << std::setprecision(16) << pc.speedOfLight()
        << " != " << std::setprecision(16) << dev << "." << endl;
    FAILMSG(msg.str());
  }

  // Stefan-Boltzmann constant

  dev = 5.6703726e-8; // W / m^2 / K^4
  dev = dev * us.p() / si.p() * pow(si.L() / us.L(), 2) * pow(si.T() / us.T(), 4);

  if (soft_equiv(pc.stefanBoltzmann(), dev, 1.0e-8)) {
    ostringstream msg;
    msg << "Scaled Stefan-Boltzmann constant looks correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Scaled Stefan-Boltzmann constant is not correct." << endl
        << "\tvalue =  " << std::setprecision(16) << pc.stefanBoltzmann()
        << " != " << std::setprecision(16) << dev << "." << endl;
    FAILMSG(msg.str());
  }

  // Gravitational constant

  dev = 6.67259e-11; // N / m^2 / kg^2
  dev = dev * us.f() / si.f() / pow(us.L() / si.L(), 2) / pow(us.M() / si.M(), 2);
  if (soft_equiv(pc.gravitationalConstant(), dev)) {
    ostringstream msg;
    msg << "Scaled Gravitational constant looks correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Scaled Gravitaitonal constant is not correct." << endl
        << "\tvalue =  " << std::setprecision(16) << pc.gravitationalConstant()
        << " != " << std::setprecision(16) << dev << "." << endl;
    FAILMSG(msg.str());
  }

  // Acceleration from Gravity

  dev = 9.80665; // m/s^2
  dev = dev * us.a() / si.a();
  if (soft_equiv(pc.accelerationFromGravity(), dev)) {
    ostringstream msg;
    msg << "Scaled Acceleration from Gravity looks correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Scaled Acceleration from Gravity is not correct." << endl
        << "\tvalue =  " << std::setprecision(16) << pc.accelerationFromGravity()
        << " != " << std::setprecision(16) << dev << "." << endl;
    FAILMSG(msg.str());
  }

  // Faraday constant

  dev = pc.avogadro() * pc.electronCharge();
  if (soft_equiv(pc.faradayConstant(), dev)) {
    ostringstream msg;
    msg << "Scaled Faraday constant looks correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Scaled Faraday constant is not correct." << endl
        << "\tvalue =  " << std::setprecision(16) << pc.faradayConstant()
        << " != " << std::setprecision(16) << dev << "." << endl;
    FAILMSG(msg.str());
  }

  // Permeability of vacuum (free space)

  dev = 4.0 * PI * 1.0e-7; // H / m
  dev = dev * si.L() / us.L();
  if (soft_equiv(pc.permeabilityOfVacuum(), dev)) {
    ostringstream msg;
    msg << "Scaled Permeability of vacuum (free space) looks correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Scaled Permeability of vacuum (free space) is not correct." << endl
        << "\tvalue =  " << std::setprecision(16) << pc.permeabilityOfVacuum()
        << " != " << std::setprecision(16) << dev << "." << endl;
    FAILMSG(msg.str());
  }

  // PERMITTIVITY OF FREE SPACE (F/M)

  dev = 1.0 / pc.permeabilityOfVacuum() / pow(pc.speedOfLight(), 2);
  if (soft_equiv(pc.permittivityOfFreeSpace(), dev)) {
    ostringstream msg;
    msg << "Scaled PERMITTIVITY OF FREE SPACE (F/M) looks correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Scaled PERMITTIVITY OF FREE SPACE (F/M) is not correct." << endl
        << "\tvalue =  " << std::setprecision(16) << pc.permittivityOfFreeSpace()
        << " != " << std::setprecision(16) << dev << "." << endl;
    FAILMSG(msg.str());
  }

  // Classical electron radius
  dev = 2.8179403227e-15; // m (NIST value)
  dev = dev * us.L() / si.L();
  if (soft_equiv(pc.classicalElectronRadius(), dev, 2e-9)) {
    ostringstream msg;
    msg << "Scaled classical electron radius looks correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Scaled classical electron radius is not correct." << endl
        << "\tvalue =  " << std::setprecision(16) << pc.classicalElectronRadius()
        << " != " << std::setprecision(16) << dev << "." << endl;
    FAILMSG(msg.str());
  }
  // test alias
  if (soft_equiv(pc.re(), dev, 2e-9)) {
    ostringstream msg;
    msg << "Scaled classical electron radius looks correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Scaled classical electron radius is not correct." << endl
        << "\tvalue =  " << std::setprecision(16) << pc.re() << " != " << std::setprecision(16)
        << dev << "." << endl;
    FAILMSG(msg.str());
  }

  // Electron mass

  dev = 9.1093897e-31; // kg
  dev = dev * us.M() / si.M();
  if (soft_equiv(pc.electronMass(), dev)) {
    ostringstream msg;
    msg << "Scaled Electron mass looks correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Scaled Electron mass is not correct." << endl
        << "\tvalue =  " << std::setprecision(16) << pc.electronMass()
        << " != " << std::setprecision(16) << dev << "." << endl;
    FAILMSG(msg.str());
  }

  // Proton mass

  dev = 1.6726231e-27; // kg
  dev = dev * us.M() / si.M();
  if (soft_equiv(pc.protonMass(), dev)) {
    ostringstream msg;
    msg << "Scaled Proton mass looks correct." << endl;
    PASSMSG(msg.str());
  } else {
    ostringstream msg;
    msg << "Scaled Proton mass is not correct." << endl
        << "\tvalue =  " << std::setprecision(16) << pc.protonMass()
        << " != " << std::setprecision(16) << dev << "." << endl;
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
