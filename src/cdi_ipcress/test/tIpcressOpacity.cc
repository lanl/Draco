//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ipcress/test/tIpcressOpacity.cc
 * \author Thomas M. Evans
 * \date   Fri Oct 26 10:50:44 2001
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi_ipcress_test.hh"
#include "cdi/OpacityCommon.hh"
#include "cdi_ipcress/IpcressFile.hh"
#include "cdi_ipcress/IpcressGrayOpacity.hh"
#include "cdi_ipcress/IpcressMultigroupOpacity.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/path.hh"
#include <algorithm>

using namespace std;

using rtt_cdi::GrayOpacity;
using rtt_cdi::MultigroupOpacity;
using rtt_cdi_ipcress::IpcressFile;
using rtt_cdi_ipcress::IpcressGrayOpacity;
using rtt_cdi_ipcress::IpcressMultigroupOpacity;
using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void file_check_two_mats(rtt_dsxx::ScalarUnitTest &ut) {
  // Ipcress data filename (IPCRESS format required)
  string op_data_file = ut.getTestSourcePath() + "two-mats.ipcress";
  FAIL_IF_NOT(rtt_dsxx::fileExists(op_data_file));

  cout << "Starting test \"file_check_two_mats\"...\n"
       << "Examining the ipcress file \"" << op_data_file << "\"...\n\n";

  // ------------------------- //
  // Create IpcressFile object //
  // ------------------------- //

  {
    shared_ptr<IpcressFile> spIF;

    // Attempt to instantiate the object.
    try {
      spIF = std::make_shared<IpcressFile>(op_data_file);
    } catch (rtt_dsxx::assertion const &error) {
      cout << "While testing tIpcressOpacity, " << error.what() << endl;
      return;
    }

    // If we make it here then spIF was successfully instantiated.
    PASSMSG("shared_ptr to new IpcressFile object created for two_mats.ipcress data.");

    // Test the IpcressFile object.
    if (spIF->getDataFilename() == op_data_file) {
      PASSMSG("IpcressFile object is now linked to the two_mats.ipcress data file.");
    } else {
      FAILMSG("IpcressFile object failed to link itself to the two_mats.ipcress data file.");
    }

    if (spIF->getNumMaterials() == 2) {
      PASSMSG("The correct number of materials was found in the two_mats.ipcress data file.");
    } else {
      FAILMSG(
          "spIF did not find the correct number of materials in the two_mats.ipcress data file.");
    }
  }

  // ---------------------- //
  // Create Opacity object. //
  // ---------------------- //

  // Material identifier. This data file has two materials: 10001 and 10002
  size_t const matid(10001);

  // Try to instantiate the Opacity object. (Rosseland, Gray Total for material 10001 in the IPCRESS
  // file pointed to by spIF).
  shared_ptr<GrayOpacity> spOp_twomat_rgt;

  try {
    shared_ptr<IpcressFile> spIF(new IpcressFile(op_data_file));
    spOp_twomat_rgt =
        std::make_shared<IpcressGrayOpacity>(spIF, matid, rtt_cdi::ROSSELAND, rtt_cdi::TOTAL);
    // spIF goes out of scope
  } catch (rtt_dsxx::assertion const &excpt) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressOpacity object for "
            << "two_mats.ipcress data.\n\t" << excpt.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we get here then the object was successfully instantiated.
  PASSMSG("shared_ptr to new Opacity object created for two_mats.ipcress data.");

  // ----------------- //
  // Gray Opacity Test //
  // ----------------- //

  double temperature = 0.1;                        // keV
  double density = 27.0;                           // g/cm^3
  double tabulatedGrayOpacity = 6.157321485417703; // cm^2/g

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_twomat_rgt, temperature, density,
                                                   tabulatedGrayOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // Check accessor functions

  rtt_cdi::OpacityModelType omt(spOp_twomat_rgt->getOpacityModelType());
  if (omt == rtt_cdi::IPCRESS_TYPE)
    PASSMSG("OpacityModelType() returned expected value.");
  else
    FAILMSG("OpacityModelType() did not return the expected value.");

  string edp(spOp_twomat_rgt->getEnergyPolicyDescriptor());
  if (edp == string("gray"))
    PASSMSG("EDP = gray");
  else
    FAILMSG("EDP != gray");

  if (!spOp_twomat_rgt->data_in_tabular_form())
    ITFAILS;

  size_t nd(spOp_twomat_rgt->getNumDensities());
  size_t nt(spOp_twomat_rgt->getNumTemperatures());
  if (nd != 5)
    FAILMSG("Found wrong number of density values.");
  if (nt != 5)
    FAILMSG("Found wrong number of temperature values.");

  vector<double> densGrid(spOp_twomat_rgt->getDensityGrid());
  vector<double> tempGrid(spOp_twomat_rgt->getTemperatureGrid());
  if (densGrid.size() != nd)
    ITFAILS;
  if (tempGrid.size() != nt)
    ITFAILS;

  vector<double> expected_densGrid = {0.01, 0.1, 1.0, 10.0, 100.0};
  vector<double> expected_tempGrid = {0.01, 0.2575, 0.505, 0.7525, 1.0};

  FAIL_IF_NOT(soft_equiv(densGrid.begin(), densGrid.end(), expected_densGrid.begin(),
                         expected_densGrid.end()));
  FAIL_IF_NOT(soft_equiv(tempGrid.begin(), tempGrid.end(), expected_tempGrid.begin(),
                         expected_tempGrid.end()));

  // --------------- //
  // MG Opacity test //
  // --------------- //

  // Create a Multigroup Rosseland Total Opacity object.
  shared_ptr<MultigroupOpacity> spOp_twomat_rtmg;

  // Try to instantiate the Opacity object.
  try {
    shared_ptr<IpcressFile> spIF(new IpcressFile(op_data_file));
    spOp_twomat_rtmg =
        std::make_shared<IpcressMultigroupOpacity>(spIF, matid, rtt_cdi::ROSSELAND, rtt_cdi::TOTAL);
  } catch (rtt_dsxx::assertion const &excpt) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressOpacity object for two_mats.ipcress data."
            << "\n\t" << excpt.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // Setup the test point.
  temperature = 0.01; // keV
  density = 2.0;      // g/cm^3

  // Check accessor functions

  omt = spOp_twomat_rtmg->getOpacityModelType();
  if (omt == rtt_cdi::IPCRESS_TYPE)
    PASSMSG("OpacityModelType() returned expected value.");
  else
    FAILMSG("OpacityModelType() did not return the expected value.");

  edp = spOp_twomat_rtmg->getEnergyPolicyDescriptor();
  if (edp == string("mg"))
    PASSMSG("EDP = mg");
  else
    FAILMSG("EDP != mg");

  FAIL_IF_NOT(spOp_twomat_rtmg->data_in_tabular_form());

  size_t numGroups = 33;
  FAIL_IF_NOT(spOp_twomat_rtmg->getNumGroups() == numGroups);

  // The solution to compare against:
  vector<double> tabulatedMGOpacity = {
      1.65474413066534,  1.25678363987902,  0.969710642123251, 0.809634824544866, 0.724209752535695,
      0.657174133970242, 0.60418067688849,  0.56102837685666,  0.526678715472818, 0.498691712866007,
      0.476029518181352, 0.458069153077808, 0.443645292733958, 0.432168463782148, 0.422976580403301,
      0.415262388915466, 0.408791600852843, 0.402940726234063, 0.397310792049752, 0.391533709574117,
      0.385148585211532, 0.377799945947677, 0.368969609795487, 0.358190516168003, 0.345144753971726,
      0.329280962243614, 0.310427778423081, 0.288524502929352, 0.263718520343973, 0.236783648295996,
      0.208267355739943, 0.179391299272131, 0.15400996177649}; // KeV, numGroups entries.

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_twomat_rtmg, temperature, density,
                                                   tabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }
}

//------------------------------------------------------------------------------------------------//
void file_check_analytic(rtt_dsxx::ScalarUnitTest &ut) {

  cout << "\nStarting test \"file_check_analytic\"...\n";

  // ----------------------------------------------------------------
  // The Opacities in this file are computed from the following analytic formula:
  //     opacity = rho * T^4,
  // rho is the density and T is the temperature.
  //
  // The grid in this data file has the following structure:
  //    T   = { 0.1, 1.0, 10.0 } keV.
  //    rho = { 0.1, 0.5, 1.0 } g/cm^3
  //    E_bounds = { 0.01, 0.03, 0.07, 0.1, 0.3, 0.7, 1.0, 3.0, 7.0, 10.0, 30.0, 70.0, 100.0 } keV.
  //-----------------------------------------------------------------

  // Ipcress data filename (IPCRESS format required)
  string op_data_file = ut.getTestSourcePath() + "analyticOpacities.ipcress";
  FAIL_IF_NOT(rtt_dsxx::fileExists(op_data_file));

  // ------------------------- //
  // Create IpcressFile object //
  // ------------------------- //

  // Create a smart pointer to a IpcressFile object
  shared_ptr<IpcressFile> spGFAnalytic;

  // Try to instantiate the object.
  try {
    spGFAnalytic = std::make_shared<rtt_cdi_ipcress::IpcressFile>(op_data_file);
  } catch (rtt_dsxx::assertion const &error) {
    FAILMSG(error.what());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we make it here then spGFAnalytic was successfully instantiated.
  PASSMSG("shared_ptr to new IpcressFile object created (spGFAnalytic).");

  // Test the IpcressFile object.
  if (spGFAnalytic->getDataFilename() == op_data_file) {
    PASSMSG("IpcressFile object is now linked to the data file.");
  } else {
    FAILMSG("IpcressFile object failed to link itself to the data file.");
  }

  if (spGFAnalytic->getNumMaterials() == 1) {
    PASSMSG("The correct number of materials was found in the data file.");
  } else {
    FAILMSG("spGFAnalytic did not find the correct number of materials in the data file.");
  }

  // --------------------- //
  // Create Opacity object //
  // --------------------- //

  // Create a smart pointer to an Opacity object.
  shared_ptr<GrayOpacity> spOp_Analytic_ragray;

  // material ID
  const int matid = 10001;

  // Try to instantiate the Opacity object.
  try {
    spOp_Analytic_ragray = std::make_shared<IpcressGrayOpacity>(
        spGFAnalytic, matid, rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION);
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressOpacity object for "
            << "analyticOpacities.ipcress data.\n\t" << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we get here then the object was successfully instantiated.
  PASSMSG("shared_ptr to new Opacity object created for analyticOpacities.ipcress.");

  // ----------------- //
  // Gray Opacity Test //
  // ----------------- //

  double temperature = 10.0; // keV
  double density = 1.0;      // g/cm^3
  double tabulatedGrayOpacity = density * pow(temperature, 4);

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_ragray, temperature, density,
                                                   tabulatedGrayOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  //---------------- //
  // MG Opacity test //
  //---------------- //

  // Create a smart pointer to an Opacity object.
  shared_ptr<MultigroupOpacity> spOp_Analytic_ramg;

  // Try to instantiate the Opacity object.
  try {
    spOp_Analytic_ramg = std::make_shared<IpcressMultigroupOpacity>(
        spGFAnalytic, matid, rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION);
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressOpacity object for "
            << "analyticOpacities.ipcress data.\n\t" << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we get here then the object was successfully instantiated.
  PASSMSG("shared_ptr to new Opacity object created for analyticOpacities.ipcress.");

  // Set up the new test problem.

  temperature = 0.3; // keV
  density = 0.7;     // g/cm^3

  // This is the solution we compare against.
  int numGroups = 12;
  vector<double> tabulatedMGOpacity(numGroups);
  for (int i = 0; i < numGroups; ++i)
    tabulatedMGOpacity[i] = density * pow(temperature, 4); // cm^2/gm

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_ramg, temperature, density,
                                                   tabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ------------------------------------------------------------ //
  // Test the Planck routines using analyticOpacities.ipcress data //
  // ------------------------------------------------------------ //

  // The Opacities in this file are computed from the following analytic formula:
  //     opacity = rho * T^4,
  // rho is the density and T is the temperature.

  // spGFAnalytic already points to the correct file so we don't repeat the coding.

  // Ditto for spOpacityAnalytic.

  // ----------------- //
  // Gray Opacity Test //
  // ----------------- //

  // Create a smart pointer to an Opacity object.
  shared_ptr<GrayOpacity> spOp_Analytic_pgray;

  // Try to instantiate the Opacity object.
  try {
    spOp_Analytic_pgray = std::make_shared<IpcressGrayOpacity>(spGFAnalytic, matid, rtt_cdi::PLANCK,
                                                               rtt_cdi::ABSORPTION);
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressOpacity object for "
            << "analyticOpacities.ipcress data.\n\t" << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we get here then the object was successfully instantiated.
  PASSMSG(
      "shared_ptr to new Gray Plank Total Opacity object created for analyticOpacities.ipcress.");

  // Setup the test problem.

  temperature = 3.0;                                     // keV
  density = 0.7;                                         // g/cm^3
  double tabulatedValue = density * pow(temperature, 4); // cm^2/g

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_pgray, temperature, density,
                                                   tabulatedValue)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // --------------- //
  // MG Opacity test //
  // --------------- //

  // Create a smart pointer to an Opacity object.
  shared_ptr<MultigroupOpacity> spOp_Analytic_pmg;

  // Try to instantiate the Opacity object.
  try {
    spOp_Analytic_pmg = std::make_shared<IpcressMultigroupOpacity>(
        spGFAnalytic, matid, rtt_cdi::PLANCK, rtt_cdi::ABSORPTION);
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressOpacity object for "
            << "analyticOpacities.ipcress data.\n\t" << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we get here then the object was successfully instantiated.
  {
    ostringstream message;
    message << "shared_ptr to new Multigroup Planck Total Opacity object created \n\t for "
            << "\"analyticOpacities.ipcress.\"";
    PASSMSG(message.str());
  }

  // Setup the test problem.

  size_t ng = 12;
  tabulatedMGOpacity.resize(ng);
  temperature = 0.4; // keV
  density = 0.22;    // g/cm^3
  for (size_t ig = 0; ig < ng; ++ig)
    tabulatedMGOpacity[ig] = density * pow(temperature, 4); // cm^2/g

  // If this test fails then stop testing.
  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_pmg, temperature, density,
                                                   tabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ------------------------ //
  // Access temperature grid. //
  // ------------------------ //

  rtt_cdi_ipcress_test::testTemperatureGridAccessor(ut, spOp_Analytic_pmg);

  // ------------------------ //
  // Access the density grid. //
  // ------------------------ //

  rtt_cdi_ipcress_test::testDensityGridAccessor(ut, spOp_Analytic_pmg);

  // ----------------------------- //
  // Access the energy boundaries. //
  // ----------------------------- //

  rtt_cdi_ipcress_test::testEnergyBoundaryAccessor(ut, spOp_Analytic_pmg);

  // ------------------------------------------------------------ //
  // Test alternate (vector-based) accessors for getGrayRosseland //
  // ------------------------------------------------------------ //

  // ---------------------- //
  // Vector of temperatures //
  // ---------------------- //

  vector<double> vtemperature(2);
  vtemperature[0] = 0.5; // keV
  vtemperature[1] = 0.7; // keV
  density = 0.35;        // g/cm^3

  vector<double> vtabulatedGrayOpacity(vtemperature.size());
  for (size_t i = 0; i < vtabulatedGrayOpacity.size(); ++i)
    vtabulatedGrayOpacity[i] = density * pow(vtemperature[i], 4);

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_ragray, vtemperature, density,
                                                   vtabulatedGrayOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ---------------------- //
  // Vector of densities    //
  // ---------------------- //

  vector<double> vdensity(3);
  temperature = 0.3; //keV
  vdensity[0] = 0.2; // g/cm^3
  vdensity[1] = 0.4; // g/cm^3
  vdensity[2] = 0.6; // g/cm^3

  vtabulatedGrayOpacity.resize(vdensity.size());
  for (size_t i = 0; i < vtabulatedGrayOpacity.size(); ++i)
    vtabulatedGrayOpacity[i] = vdensity[i] * pow(temperature, 4);

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_ragray, temperature, vdensity,
                                                   vtabulatedGrayOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // -------------------------------------------------------- //
  // Test alternate (vector-based) accessors for getGray Planck //
  // -------------------------------------------------------- //

  // ---------------------- //
  // Vector of temperatures //
  // ---------------------- //

  vtemperature.resize(2);
  vtemperature[0] = 0.5; // keV
  vtemperature[1] = 0.7; // keV
  density = 0.35;        // g/cm^3q

  vtabulatedGrayOpacity.resize(vtemperature.size());
  for (size_t i = 0; i < vtabulatedGrayOpacity.size(); ++i)
    vtabulatedGrayOpacity[i] = density * pow(vtemperature[i], 4);

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_pgray, vtemperature, density,
                                                   vtabulatedGrayOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ------------------- //
  // Vector of densities //
  // ------------------- //

  vdensity.resize(3);
  temperature = 0.3; //keV
  vdensity[0] = 0.2; // g/cm^3
  vdensity[1] = 0.4; // g/cm^3
  vdensity[2] = 0.6; // g/cm^3

  vtabulatedGrayOpacity.resize(vdensity.size());
  for (size_t i = 0; i < vtabulatedGrayOpacity.size(); ++i)
    vtabulatedGrayOpacity[i] = vdensity[i] * pow(temperature, 4);

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_pgray, temperature, vdensity,
                                                   vtabulatedGrayOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ---------------------------------------------------------- //
  // Test alternate (vector-based) accessors for getMGRosseland //
  // ---------------------------------------------------------- //

  // ---------------------- //
  // Vector of temperatures //
  // ---------------------- //

  vtemperature.resize(2);
  vtemperature[0] = 0.5; // keV
  vtemperature[1] = 0.7; // keV
  density = 0.35;        // g/cm^3
  ng = spOp_Analytic_ramg->getNumGroupBoundaries() - 1;

  vector<vector<double>> vtabulatedMGOpacity(vtemperature.size());
  for (size_t i = 0; i < vtemperature.size(); ++i) {
    vtabulatedMGOpacity[i].resize(ng);
    for (size_t ig = 0; ig < ng; ++ig)
      vtabulatedMGOpacity[i][ig] = density * pow(vtemperature[i], 4);
  }

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_ramg, vtemperature, density,
                                                   vtabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ------------------- //
  // Vector of densities //
  // ------------------- //

  vdensity.resize(2);
  vdensity[0] = 0.3; // g/cm^3
  vdensity[1] = 0.7; // g/cm^3
  temperature = 7.0; // keV
  ng = spOp_Analytic_ramg->getNumGroupBoundaries() - 1;

  vtabulatedMGOpacity.resize(vdensity.size());
  for (size_t i = 0; i < vdensity.size(); ++i) {
    vtabulatedMGOpacity[i].resize(ng);
    for (size_t ig = 0; ig < ng; ++ig)
      vtabulatedMGOpacity[i][ig] = vdensity[i] * pow(temperature, 4);
  }

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_ramg, temperature, vdensity,
                                                   vtabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ------------------------------------------------------ //
  // Test alternate (vector-based) accessors for getMG Planck //
  // ------------------------------------------------------ //

  // ---------------------- //
  // Vector of temperatures //
  // ---------------------- //

  vtemperature.resize(2);
  vtemperature[0] = 0.5; // keV
  vtemperature[1] = 0.7; // keV
  density = 0.35;        // g/cm^3
  ng = spOp_Analytic_pmg->getNumGroupBoundaries() - 1;

  vtabulatedMGOpacity.resize(vtemperature.size());
  for (size_t i = 0; i < vtemperature.size(); ++i) {
    vtabulatedMGOpacity[i].resize(ng);
    for (size_t ig = 0; ig < ng; ++ig)
      vtabulatedMGOpacity[i][ig] = density * pow(vtemperature[i], 4);
  }

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_pmg, vtemperature, density,
                                                   vtabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ------------------- //
  // Vector of densities //
  // ------------------- //

  vdensity.resize(2);
  vdensity[0] = 0.3; // g/cm^3
  vdensity[1] = 0.7; // g/cm^3
  temperature = 7.0; // keV
  ng = spOp_Analytic_pmg->getNumGroupBoundaries() - 1;

  vtabulatedMGOpacity.resize(vdensity.size());
  for (size_t i = 0; i < vdensity.size(); ++i) {
    vtabulatedMGOpacity[i].resize(ng);
    for (size_t ig = 0; ig < ng; ++ig)
      vtabulatedMGOpacity[i][ig] = vdensity[i] * pow(temperature, 4);
  }

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_pmg, temperature, vdensity,
                                                   vtabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }
}

//------------------------------------------------------------------------------------------------//
void check_ipcress_stl_accessors(rtt_dsxx::ScalarUnitTest &ut) {

  cout << "\nStarting test \"check_ipcress_stl_accessors\"...\n";

  // Ipcress data filename (IPCRESS format required)
  string op_data_file = ut.getTestSourcePath() + "analyticOpacities.ipcress";

  // ------------------------- //
  // Create IpcressFile object //
  // ------------------------- //

  // Create a smart pointer to a IpcressFile object
  shared_ptr<IpcressFile> spGFAnalytic;

  // Try to instantiate the object.
  try {
    spGFAnalytic = std::make_shared<rtt_cdi_ipcress::IpcressFile>(op_data_file);
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    FAILMSG(error.what());
    FAILMSG("Aborting tests.");
    return;
  }

  // material ID
  const int matid = 10001;

  // -------------------------------------- //
  // Test the STL-like getOpacity accessor  //
  // Using const iterators for Gray objects //
  // -------------------------------------- //

  // These accessors are only available in IpcressOpacity objects so the shared_ptr must be
  // templated on IpcressGrayOpacity and not on cdi/GrayOpacity.

  // Create a new smart pointer to a IpcressGrayOpacity object.
  shared_ptr<IpcressGrayOpacity> spGGOp_Analytic_ra;

  // try to instantiate the Opacity object.
  try {
    spGGOp_Analytic_ra = std::make_shared<IpcressGrayOpacity>(
        spGFAnalytic, matid, rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION);
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressGrayOpacity object for\n\t"
            << "analyticOpacities.ipcress data (shared_ptr not templated on cdi/GrayOpacity).\n\t"
            << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  vector<double> vdensity;
  vector<double> vtemperature;
  vector<double> vtabulatedGrayOpacity;

  double density = 0.35;    // g/cm^3
  double temperature = 7.0; // kev

  // Setup the temperature and density parameters for this test.
  vdensity.resize(6);
  vtemperature.resize(6);

  // (temperature,density) tuples.

  vtemperature[0] = 0.5; // keV
  vdensity[0] = 0.2;     // g/cm^3

  vtemperature[1] = 0.7; // keV
  vdensity[1] = 0.2;     // g/cm^3

  vtemperature[2] = 0.5; // keV
  vdensity[2] = 0.4;     // g/cm^3

  vtemperature[3] = 0.7; // keV
  vdensity[3] = 0.4;     // g/cm^3

  vtemperature[4] = 0.5; // keV
  vdensity[4] = 0.6;     // g/cm^3

  vtemperature[5] = 0.7; // keV
  vdensity[5] = 0.6;     // g/cm^3

  // we want to test the const_iterator version of getOpacity() so we need to create const vectors
  // with the tuple data.
  const vector<double> cvdensity = vdensity;
  const vector<double> cvtemperature = vtemperature;

  size_t nt = cvtemperature.size();
  size_t nd = cvdensity.size();

  // Here is the reference solution
  vtabulatedGrayOpacity.resize(nt);
  for (size_t i = 0; i < nt; ++i)
    vtabulatedGrayOpacity[i] = cvdensity[i] * pow(cvtemperature[i], 4);

  // Here is the solution from Ipcress
  vector<double> graOp(nt);
  spGGOp_Analytic_ra->getOpacity(cvtemperature.begin(), cvtemperature.end(), cvdensity.begin(),
                                 cvdensity.end(), graOp.begin());
  if (rtt_dsxx::soft_equiv(graOp, vtabulatedGrayOpacity)) {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor() << " opacity computation was good for \n\t"
            << spGGOp_Analytic_ra->getDataFilename()
            << " (const-iterator accessor, temp x density).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor() << " opacity value is out of spec. for \n\t"
            << spGGOp_Analytic_ra->getDataFilename()
            << " (non-const-iterator accessor, temp x density).";
    FAILMSG(message.str());
  }

  // ------------------------------------- //
  // Test the STL-like getOpacity accessor //
  // Using non-const iterator              //
  // ------------------------------------- //

  spGGOp_Analytic_ra->getOpacity(vtemperature.begin(), vtemperature.end(), vdensity.begin(),
                                 vdensity.end(), graOp.begin());
  if (rtt_dsxx::soft_equiv(graOp, vtabulatedGrayOpacity)) {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor() << " opacity computation was good for \n\t"
            << spGGOp_Analytic_ra->getDataFilename()
            << " (non-const-iterator accessor, temp x density).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor() << " opacity value is out of spec. for \n\t"
            << spGGOp_Analytic_ra->getDataFilename()
            << " (non-const-iterator accessor, temp x density).";
    FAILMSG(message.str());
  }

  // ------------------------------------- //
  // Test the STL-like getOpacity accessor //
  // const iterator (temperature only)     //
  // ------------------------------------- //

  graOp.resize(nt);
  vtabulatedGrayOpacity.resize(nt);
  for (size_t it = 0; it < nt; ++it)
    vtabulatedGrayOpacity[it] = density * pow(vtemperature[it], 4);

  spGGOp_Analytic_ra->getOpacity(cvtemperature.begin(), cvtemperature.end(), density,
                                 graOp.begin());
  if (rtt_dsxx::soft_equiv(graOp, vtabulatedGrayOpacity)) {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor() << " opacity computation was good for \n\t"
            << spGGOp_Analytic_ra->getDataFilename() << " (const iterator accessor, vtemps).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor() << " opacity value is out of spec. for \n\t"
            << spGGOp_Analytic_ra->getDataFilename() << " (const iterator accessor, vtemps).";
    FAILMSG(message.str());
  }

  // ------------------------------------- //
  // Test the STL-like getOpacity accessor //
  // const iterator ( density only)        //
  // ------------------------------------- //

  graOp.resize(nd);
  vtabulatedGrayOpacity.resize(nd);
  for (size_t id = 0; id < nd; ++id)
    vtabulatedGrayOpacity[id] = vdensity[id] * pow(temperature, 4);

  spGGOp_Analytic_ra->getOpacity(temperature, cvdensity.begin(), cvdensity.end(), graOp.begin());
  if (rtt_dsxx::soft_equiv(graOp, vtabulatedGrayOpacity)) {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor() << " opacity computation was good for \n\t"
            << spGGOp_Analytic_ra->getDataFilename() << " (const iterator accessor, vdensity).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor() << " opacity value is out of spec. for \n\t"
            << spGGOp_Analytic_ra->getDataFilename() << " (const iterator accessor, vdensity).";
    FAILMSG(message.str());
  }

  // -------------------------------------- //
  // Test the STL-like getOpacity accessor  //
  // Using const iterators for MG objects   //
  // -------------------------------------- //

  // These accessors are only available in IpcressOpacity objects so the shared_ptr must be
  // templated on IpcressMultigroupOpacity and not on cdi/MultigroupOpacity.

  // Create a new smart pointer to a IpcressGrayOpacity object.
  shared_ptr<IpcressMultigroupOpacity> spGMGOp_Analytic_ra;

  // try to instantiate the Opacity object.
  try {
    spGMGOp_Analytic_ra = std::make_shared<IpcressMultigroupOpacity>(
        spGFAnalytic, matid, rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION);
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressGrayOpacity object for\n\t"
            << "analyticOpacities.ipcress data (shared_ptr not templated on cdi/GrayOpacity).\n\t"
            << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we get here then the object was successfully instantiated.
  PASSMSG("shared_ptr to new Opacity object created for analyticOpacities.ipcress.");

  // Here is the reference solution
  size_t ng = spGMGOp_Analytic_ra->getNumGroupBoundaries() - 1;
  vector<double> vtabulatedOpacity(ng * nt);

  for (size_t i = 0; i < nt; ++i)
    for (size_t ig = 0; ig < ng; ++ig)
      vtabulatedOpacity[i * ng + ig] = cvdensity[i] * pow(cvtemperature[i], 4);

  // Here is the solution from Ipcress
  vector<double> mgOp(nt * ng);
  spGMGOp_Analytic_ra->getOpacity(cvtemperature.begin(), cvtemperature.end(), cvdensity.begin(),
                                  cvdensity.end(), mgOp.begin());

  if (rtt_dsxx::soft_equiv(mgOp, vtabulatedOpacity)) {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor() << " opacity computation was good for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename()
            << " (const-iterator accessor, temp x density).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor() << " opacity value is out of spec. for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename()
            << " (non-const-iterator accessor, temp x density).";
    FAILMSG(message.str());
  }

  // --------------------------------------- //
  // Test the STL-like getOpacity accessor   //
  // Using non-const iterator for MG objects //
  // --------------------------------------- //

  // clear old data
  for (size_t i = 0; i < nt * ng; ++i)
    mgOp[i] = 0.0;

  // use Ipcress to obtain new data
  spGMGOp_Analytic_ra->getOpacity(vtemperature.begin(), vtemperature.end(), vdensity.begin(),
                                  vdensity.end(), mgOp.begin());

  // compare the results to the reference solution and report our findings.
  if (rtt_dsxx::soft_equiv(mgOp, vtabulatedOpacity)) {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor() << " opacity computation was good for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename()
            << " (non-const-iterator accessor, temp x density).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor() << " opacity value is out of spec. for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename()
            << " (non-const-iterator accessor, temp x density).";
    FAILMSG(message.str());
  }

  // ------------------------------------------------ //
  // Test the STL-like getOpacity accessor            //
  // const iterator (temperature only) for MG data    //
  // ------------------------------------------------ //

  // clear old data
  for (size_t i = 0; i < nt * ng; ++i)
    mgOp[i] = 0.0;

  // Calculate the reference solution.
  for (size_t it = 0; it < nt; ++it)
    for (size_t ig = 0; ig < ng; ++ig)
      vtabulatedOpacity[it * ng + ig] = density * pow(vtemperature[it], 4);

  // Obtain new solution
  spGMGOp_Analytic_ra->getOpacity(cvtemperature.begin(), cvtemperature.end(), density,
                                  mgOp.begin());

  // Compare solutions and report the results.
  if (rtt_dsxx::soft_equiv(mgOp, vtabulatedOpacity)) {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor() << " opacity computation was good for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename() << " (const iterator accessor, vtemps).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor() << " opacity value is out of spec. for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename() << " (const iterator accessor, vtemps).";
    FAILMSG(message.str());
  }

  // ------------------------------------------ //
  // Test the STL-like getOpacity accessor      //
  // const iterator ( density only) for MG data //
  // ------------------------------------------ //

  // clear old data
  for (size_t i = 0; i < nd * ng; ++i)
    mgOp[i] = 0.0;

  // Calculate the reference solution.
  for (size_t id = 0; id < nd; ++id)
    for (size_t ig = 0; ig < ng; ++ig)
      vtabulatedOpacity[id * ng + ig] = vdensity[id] * pow(temperature, 4);

  // Obtain new solution
  spGMGOp_Analytic_ra->getOpacity(temperature, cvdensity.begin(), cvdensity.end(), mgOp.begin());

  // Compare solutions and report the results.
  if (rtt_dsxx::soft_equiv(mgOp, vtabulatedOpacity)) {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor() << " opacity computation was good for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename() << " (const iterator accessor, vdensity).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor() << " opacity value is out of spec. for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename() << " (const iterator accessor, vdensity).";
    FAILMSG(message.str());
  }
}

//------------------------------------------------------------------------------------------------//
void gray_opacity_packing_test(rtt_dsxx::ScalarUnitTest &ut) {

  cout << "\nStarting test \"gray_opacity_packing_test\"...\n";

  // Ipcress data filename (IPCRESS format required)
  string const op_data_file = ut.getTestSourcePath() + "analyticOpacities.ipcress";
  FAIL_IF_NOT(rtt_dsxx::fileExists(op_data_file));

  // ------------------------- //
  // Create IpcressFile object //
  // ------------------------- //

  // Create a smart pointer to a IpcressFile object
  shared_ptr<IpcressFile> spGFAnalytic;

  // Try to instantiate the object.
  try {
    spGFAnalytic = std::make_shared<rtt_cdi_ipcress::IpcressFile>(op_data_file);
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    FAILMSG(error.what());
    FAILMSG("Aborting tests.");
    return;
  }

  // Create a smart pointer to an Opacity object.
  shared_ptr<GrayOpacity> spOp_Analytic_ragray;

  // material ID
  int const matid = 10001;

  spOp_Analytic_ragray = std::make_shared<IpcressGrayOpacity>(
      spGFAnalytic, matid, rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION);
  // pack up the opacity
  vector<char> packed = spOp_Analytic_ragray->pack();

  // make a new IpcressGrayOpacity from packed data
  shared_ptr<GrayOpacity> unpacked_opacity;

  // Try to instantiate the Opacity object.
  try {
    unpacked_opacity = std::make_shared<IpcressGrayOpacity>(packed);
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    message << "Failed to create shared_ptr to unpacked IpcressOpacity object for "
            << "analyticOpacities.ipcress data.\n\t" << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // some simple tests
  FAIL_IF_NOT(unpacked_opacity->getDataFilename() == op_data_file);
  FAIL_IF_NOT(unpacked_opacity->getReactionType() == rtt_cdi::ABSORPTION);
  FAIL_IF_NOT(unpacked_opacity->getModelType() == rtt_cdi::ROSSELAND);

  // ----------------- //
  // Gray Opacity Test //
  // ----------------- //

  double temperature = 10.0; // keV
  double density = 1.0;      // g/cm^3
  double tabulatedGrayOpacity = density * pow(temperature, 4);

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, unpacked_opacity, temperature, density,
                                                   tabulatedGrayOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // try to unpack gray opacity as multigroup opacity
  bool caught = false;
  try {
    shared_ptr<MultigroupOpacity> opacity(new IpcressMultigroupOpacity(packed));
  } catch (rtt_dsxx::assertion const &error) {
    caught = true;
    ostringstream message;
    message << "Good, we caught the following assertion, \n" << error.what();
    PASSMSG(message.str());
  }
  if (!caught) {
    FAILMSG("Failed to catch an illegal packing asserion.");
  }
}

//------------------------------------------------------------------------------------------------//
void mg_opacity_packing_test(rtt_dsxx::ScalarUnitTest &ut) {

  cout << "\nStarting test \"mg_opacity_packing_test\"...\n";

  // Ipcress data filename (IPCRESS format required)
  string const op_data_file = ut.getTestSourcePath() + "analyticOpacities.ipcress";
  FAIL_IF_NOT(rtt_dsxx::fileExists(op_data_file));

  // ------------------------- //
  // Create IpcressFile object //
  // ------------------------- //

  // Create a smart pointer to a IpcressFile object
  shared_ptr<IpcressFile> spGFAnalytic;

  // Try to instantiate the object.
  try {
    spGFAnalytic = std::make_shared<rtt_cdi_ipcress::IpcressFile>(op_data_file);
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    FAILMSG(error.what());
    FAILMSG("Aborting tests.");
    return;
  }

  // material ID
  const int matid = 10001;

  //---------------- //
  // MG Opacity test //
  //---------------- //

  // Create a smart pointer to an Opacity object.
  shared_ptr<MultigroupOpacity> spOp_Analytic_pmg;

  spOp_Analytic_pmg = std::make_shared<IpcressMultigroupOpacity>(
      spGFAnalytic, matid, rtt_cdi::PLANCK, rtt_cdi::ABSORPTION);
  vector<char> packed = spOp_Analytic_pmg->pack();

  // make a new IpcressGrayOpacity from packed data
  shared_ptr<MultigroupOpacity> unpacked_opacity;

  // Try to instantiate the Opacity object.
  try {
    unpacked_opacity = std::make_shared<IpcressMultigroupOpacity>(packed);
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    message << "Failed to create shared_ptr to unpacked IpcressOpacity object for "
            << "analyticOpacities.ipcress data.\n\t" << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // some simple tests
  FAIL_IF_NOT(unpacked_opacity->getDataFilename() == op_data_file);
  FAIL_IF_NOT(unpacked_opacity->getReactionType() == rtt_cdi::ABSORPTION);
  FAIL_IF_NOT(unpacked_opacity->getModelType() == rtt_cdi::PLANCK);

  // Setup the test problem.

  size_t ng(12);
  vector<double> tabulatedMGOpacity(ng);
  double temperature = 0.4; // keV
  double density = 0.22;    // g/cm^3
  for (size_t ig = 0; ig < ng; ++ig)
    tabulatedMGOpacity[ig] = density * pow(temperature, 4); // cm^2/g

  // If this test fails then stop testing.
  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, unpacked_opacity, temperature, density,
                                                   tabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ------------------------ //
  // Access temperature grid. //
  // ------------------------ //

  rtt_cdi_ipcress_test::testTemperatureGridAccessor(ut, unpacked_opacity);

  // ------------------------ //
  // Access the density grid. //
  // ------------------------ //

  rtt_cdi_ipcress_test::testDensityGridAccessor(ut, unpacked_opacity);

  // ----------------------------- //
  // Access the energy boundaries. //
  // ----------------------------- //

  rtt_cdi_ipcress_test::testEnergyBoundaryAccessor(ut, unpacked_opacity);

  // try to unpack multigroup as gray opacity
  bool caught = false;
  try {
    shared_ptr<GrayOpacity> opacity(new IpcressGrayOpacity(packed));
  } catch (rtt_dsxx::assertion const &error) {
    caught = true;
    ostringstream message;
    message << "Good, we caught the following assertion, \n" << error.what();
    PASSMSG(message.str());
  }
  if (!caught) {
    FAILMSG("Failed to catch an illegal packing asserion.");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    file_check_two_mats(ut);
    file_check_analytic(ut);
    check_ipcress_stl_accessors(ut);

    gray_opacity_packing_test(ut);
    mg_opacity_packing_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tIpcressOpacity.cc
//------------------------------------------------------------------------------------------------//
