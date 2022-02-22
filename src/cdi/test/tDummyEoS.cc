//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/test/tDummyEoS.cc
 * \author Thomas M. Evans
 * \date   Tue Oct  9 10:52:50 2001
 * \brief  EoS class test.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "DummyEoS.hh"
#include "cdi/EoS.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <memory>
#include <sstream>

using namespace std;

using rtt_cdi::EoS;
using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_EoS(rtt_dsxx::UnitTest &ut) {
  // --------------------- //
  // Create an EoS object. //
  // --------------------- //

  // The smart pointer points to a generic EoS object.
  std::shared_ptr<EoS> spEoS;

  // The actual instatniate is specific (dummyEoS).
  if ((spEoS = std::make_shared<rtt_cdi_test::DummyEoS>()), spEoS) {
    // If we get here then the object was successfully instantiated.
    PASSMSG("Smart Pointer to new EoS object created.");
  } else {
    FAILMSG("Unable to create a Smart Pointer to new EoS object.");
  }

  // --------- //
  // EoS Tests //
  // --------- //

  double temperature = 5800.0;                                                     // Kelvin
  double density = 27.0;                                                           // g/cm^3
  double tabulatedSpecificElectronInternalEnergy = temperature + 1000.0 * density; // kJ/g

  double seie = spEoS->getSpecificElectronInternalEnergy(temperature, density);

  if (soft_equiv(seie, tabulatedSpecificElectronInternalEnergy)) {
    ostringstream message;
    message << "The getSpecificElectronInternalEnergy( dbl, dbl) "
            << "request returned the expected value.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "The getSpecificElectronInternalEnergy( dbl, dbl) "
            << "request returned a value that is out of spec.";
    FAILMSG(message.str());
  }

  // try using a vectors of temps. and densities vtemperature.size() == vdensity.size()

  std::vector<double> vtemperature(3);
  vtemperature[0] = 5000.0; // Kelvin
  vtemperature[1] = 7000.0; // Kelvin
  vtemperature[2] = 3000.0; // Kelvin

  std::vector<double> vdensity(3);
  vdensity[0] = 0.35; // g/cm^3
  vdensity[1] = 1.0;  // g/cm^3
  vdensity[2] = 9.8;  // g/mcm^3

  // Retrieve electron based heat capacities.
  std::vector<double> vRefCve(vtemperature.size());
  for (size_t i = 0; i < vtemperature.size(); ++i)
    vRefCve[i] = vtemperature[i] + vdensity[i] / 1000.0;

  std::vector<double> vCve = spEoS->getElectronHeatCapacity(vtemperature, vdensity);

  if (soft_equiv(vCve.begin(), vCve.end(), vRefCve.begin(), vRefCve.end())) {
    ostringstream message;
    message << "The getElectronHeatCapacity( vect, vect ) request"
            << " returned the expected values.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "The getElectronHeatCapacity( vect, vect ) request"
            << " returned values that are out of spec.";
    FAILMSG(message.str());
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_EoS(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tDummyEoS.cc
//------------------------------------------------------------------------------------------------//
