//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ipcress/test/tIpcressWithCDI.cc
 * \author Thomas M. Evans
 * \date   Mon Oct 29 17:16:32 2001
 * \brief  Ipcress + CDI test.
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "cdi_ipcress_test.hh"
#include "cdi/CDI.hh" // this includes everything from CDI
#include "cdi_ipcress/IpcressFile.hh"
#include "cdi_ipcress/IpcressGrayOpacity.hh"
#include "cdi_ipcress/IpcressMultigroupOpacity.hh"
#include "ds++/Release.hh"

using namespace std;

using rtt_cdi::CDI;
using rtt_cdi::GrayOpacity;
using rtt_cdi::MultigroupOpacity;
using rtt_cdi_ipcress::IpcressFile;
using rtt_cdi_ipcress::IpcressGrayOpacity;
using rtt_cdi_ipcress::IpcressMultigroupOpacity;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_ipcress_CDI(rtt_dsxx::ScalarUnitTest &ut) {

  // ----------------------------------------------- //
  // Test the data file "analyticOpacities.ipcress"  //
  // ----------------------------------------------- //

  // -----------------------------------------------------------------
  // The Opacities in this file are computed from the following analytic
  // formula:
  //     opacity = rho * T^4,
  // rho is the density and T is the temperature.
  //
  // The grid in this data file has the following structure:
  //    T   = { 0.1, 1.0, 10.0 } keV.
  //    rho = { 0.1, 0.5, 1.0 } g/cm^3
  //    E_bounds = { 0.01, 0.03, 0.07, 0.1, 0.3, 0.7, 1.0, 3.0, 7.0
  //                 10.0, 30.0, 70.0 100.0 } keV.
  //-----------------------------------------------------------------

  // Ipcress data filename (IPCRESS format required)
  string op_data_file = ut.getTestSourcePath() + "analyticOpacities.ipcress";

  // ------------------------- //
  // Create IpcressFile object //
  // ------------------------- //

  // Create a smart pointer to a IpcressFile object
  shared_ptr<const IpcressFile> spGFAnalytic;

  // Try to instantiate the object.
  try {
    spGFAnalytic = std::make_shared<const IpcressFile>(op_data_file);
  } catch (rtt_dsxx::assertion const &excpt) {
    FAILMSG(excpt.what());
    ostringstream message;
    message << "Aborting tests because unable to instantiate IpcressFile "
            << "object";
    FAILMSG(message.str());
    return;
  }

  // If we make it here then spGFAnalytic was successfully instantiated.
  PASSMSG("shared_ptr to new IpcressFile object created (spGFAnalytic).");

  // ----------------------------------- //
  // Create a IpcressGrayOpacity object. //
  // ----------------------------------- //

  // Material identifier.  This data file has two materials: First has mat tag
  // "10001".
  int const matid = 10001;

  // Create a smart pointer to an opacity object.
  shared_ptr<const GrayOpacity> spOp_Analytic_ragray;

  // Try to instantiate the opacity object.
  try {
    spOp_Analytic_ragray = std::make_shared<const IpcressGrayOpacity>(
        spGFAnalytic, matid, rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION);
  } catch (rtt_dsxx::assertion const &excpt) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressGrayOpacity object for "
            << "two-mats.ipcress data.\n\t" << excpt.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we get here then the object was successfully instantiated.
  {
    ostringstream message;
    message << "shared_ptr to new IpcressGrayOpacity object created for "
            << "analyticOpacities.ipcress.";
    PASSMSG(message.str());
  }

  // ----------------- //
  // Create CDI object //
  // ----------------- //

  shared_ptr<CDI> spCDI_Analytic;
  if ((spCDI_Analytic = std::make_shared<CDI>()), spCDI_Analytic) {
    ostringstream message;
    message << "shared_ptr to CDI object created successfully (GrayOpacity).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Failed to create shared_ptr to CDI object (GrayOpacity).";
    FAILMSG(message.str());
  }

  // ------------------ //
  // Gray Opacity Tests //
  // ------------------ //

  // set the gray opacity
  spCDI_Analytic->setGrayOpacity(spOp_Analytic_ragray);

  double temperature = 10.0;                                   // keV
  double density = 1.0;                                        // g/cm^3
  double tabulatedGrayOpacity = density * pow(temperature, 4); // cm^2/g

  rtt_cdi::Model r = rtt_cdi::ROSSELAND;
  rtt_cdi::Reaction a = rtt_cdi::ABSORPTION;

  double opacity = spCDI_Analytic->gray(r, a)->getOpacity(temperature, density);

  if (rtt_dsxx::soft_equiv(opacity, tabulatedGrayOpacity)) {
    ostringstream message;
    message << spCDI_Analytic->gray(r, a)->getDataDescriptor()
            << " getOpacity computation was good.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spCDI_Analytic->gray(r, a)->getDataDescriptor()
            << " getOpacity value is out of spec.";
    FAILMSG(message.str());
  }

  // try using a vector of temps.

  vector<double> vtemperature(2);
  vtemperature[0] = 0.5; // keV
  vtemperature[1] = 0.7; // keV
  density = 0.35;        // g/cm^3
  vector<double> vRefOpacity(vtemperature.size());
  for (size_t i = 0; i < vtemperature.size(); ++i)
    vRefOpacity[i] = density * pow(vtemperature[i], 4);

  vector<double> vOpacity = spCDI_Analytic->gray(r, a)->getOpacity(vtemperature, density);

  if (rtt_dsxx::soft_equiv(vOpacity, vRefOpacity)) {
    ostringstream message;
    message << spCDI_Analytic->gray(r, a)->getDataDescriptor()
            << " getOpacity computation was good for a vector of temps.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spCDI_Analytic->gray(r, a)->getDataDescriptor()
            << " getOpacity value is out of spec. for a vector of temps.";
    FAILMSG(message.str());
  }

  // STL-like accessor
  //
  // The virtual base class does not support STL-like accessors so we don't test
  // this feature.

  // ----------------------------------------- //
  // Create a IpcressMultigorupOpacity object. //
  // ----------------------------------------- //

  // Create a smart pointer to an opacity object.
  shared_ptr<const MultigroupOpacity> spOp_Analytic_ramg;

  // Try to instantiate the opacity object.
  try {
    spOp_Analytic_ramg = std::make_shared<const IpcressMultigroupOpacity>(
        spGFAnalytic, matid, rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION);
  } catch (rtt_dsxx::assertion const &excpt) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressMultigroupOpacity "
            << "object for two-mats.ipcress data.\n\t" << excpt.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we get here then the object was successfully instantiated.
  {
    ostringstream message;
    message << "shared_ptr to new Ipcress multigroup opacity object created"
            << "\n\tfor analyticOpacities.ipcress.";
    PASSMSG(message.str());
  }

  // ----------------------------------------------- //
  // Create a new CDI that has both Gray and MG data //
  // ----------------------------------------------- //

  // Add the multigroup opacity object to this CDI object.

  spCDI_Analytic->setMultigroupOpacity(spOp_Analytic_ramg);

  // --------------- //
  // MG Opacity test //
  // --------------- //

  // Set up the new test problem.

  temperature = 0.3; // keV
  density = 0.7;     // g/cm^3

  // This is the solution we compare against.
  int numGroups = 12;
  vector<double> tabulatedMGOpacity(numGroups);
  for (int i = 0; i < numGroups; ++i)
    tabulatedMGOpacity[i] = density * pow(temperature, 4); // cm^2/gm

  // Request the multigroup opacity vector.
  vector<double> mgOpacity = spCDI_Analytic->mg(r, a)->getOpacity(temperature, density);

  if (rtt_dsxx::soft_equiv(mgOpacity, tabulatedMGOpacity)) {
    ostringstream message;
    message << spCDI_Analytic->mg(r, a)->getDataDescriptor() << " getOpacity computation was good.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spCDI_Analytic->mg(r, a)->getDataDescriptor() << " getOpacity value is out of spec.";
    FAILMSG(message.str());
  }

  // Finally lets check to see if CDI catches some inappropriate accesses

  // multigroup access
  bool caught = false;
  try {
    spCDI_Analytic->mg(r, rtt_cdi::SCATTERING);
  } catch (const rtt_dsxx::assertion & /*excpt*/) {
    PASSMSG("Good, caught illegal accessor to CDI-mg().");
    caught = true;
  }
  if (!caught)
    FAILMSG("Failed to catch illegal accessor to CDI-mg().");

  // gray access
  caught = false;
  try {
    spCDI_Analytic->gray(rtt_cdi::ANALYTIC, a);
  } catch (const rtt_dsxx::assertion & /*excpt*/) {
    PASSMSG("Good, caught illegal accessor to CDI-gray().");
    caught = true;
  }
  if (!caught)
    FAILMSG("Failed to catch illegal accessor to CDI-gray().");
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_ipcress_CDI(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tIpcressWithCDI.cc
//------------------------------------------------------------------------------------------------//
