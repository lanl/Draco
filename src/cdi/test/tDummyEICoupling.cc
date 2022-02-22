//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/test/tDummyEICoupling.cc
 * \author Mathew Cleveland
 * \date   March 2019
 * \brief  EICoupling class test.
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "DummyEICoupling.hh"
#include "cdi/EICoupling.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <memory>
#include <sstream>

using namespace std;

using rtt_cdi::EICoupling;
using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_EICoupling(rtt_dsxx::UnitTest &ut) {
  // ---------------------------- //
  // Create an EICoupling object. //
  // ---------------------------- //

  // The smart pointer points to a generic EICoupling object.
  std::shared_ptr<EICoupling> spEICoupling;

  // The actual instatniate is specific (dummyEoS).
  if ((spEICoupling = std::make_shared<rtt_cdi_test::DummyEICoupling>()), spEICoupling) {
    // If we get here then the object was successfully instantiated.
    PASSMSG("Smart Pointer to new EICoupling object created.");
  } else {
    FAILMSG("Unable to create a Smart Pointer to new EICoupling object.");
  }

  // --------- //
  // EoS Tests //
  // --------- //

  double etemperature = 1.0; // keV
  double itemperature = 2.0; // keV
  double density = 3.0;      // g/cm^3
  double w_e = 4.0;
  double w_i = 5.0;
  double ei_coupling_ref =
      etemperature + 10.0 * itemperature + 100.0 * density + 1000.0 * w_e + 10000.0 * w_i;

  double ei_coupling =
      spEICoupling->getElectronIonCoupling(etemperature, itemperature, density, w_e, w_i);

  if (soft_equiv(ei_coupling, ei_coupling_ref)) {
    ostringstream message;
    message << "The getElectronIonCoupling scalar request"
            << "request returned the expected value.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "The getElectonIonCoupling scalar request"
            << "request returned a value that is out of spec.";
    FAILMSG(message.str());
  }

  std::vector<double> vetemperature{1.0, 3.0, 5.0}; // Kelvin
  std::vector<double> vitemperature{2.0, 4.0, 6.0}; // Kelvin

  std::vector<double> vdensity{1.0, 2.0, 3.8}; // g/cm^3

  std::vector<double> vw_e{1.0, 2.0, 3.8};

  std::vector<double> vw_i{1.0, 2.0, 3.8};

  // Retrieve electron based heat capacities.
  std::vector<double> vRefEICoupling(vetemperature.size());
  for (size_t i = 0; i < vetemperature.size(); ++i)
    vRefEICoupling[i] = vetemperature[i] + 10.0 * vitemperature[i] + 100.0 * vdensity[i] +
                        1000.0 * vw_e[i] + 10000.0 * vw_i[i];

  std::vector<double> vEICoupling =
      spEICoupling->getElectronIonCoupling(vetemperature, vitemperature, vdensity, vw_e, vw_i);

  if (soft_equiv(vEICoupling.begin(), vEICoupling.end(), vRefEICoupling.begin(),
                 vRefEICoupling.end())) {
    ostringstream message;
    message << "The getElectronIonCoupling vector request"
            << " returned the expected values.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "The getElectronIonCoupling vector request"
            << " returned values that are out of spec.";
    FAILMSG(message.str());
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_EICoupling(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tDummyEICoupling.cc
//------------------------------------------------------------------------------------------------//
