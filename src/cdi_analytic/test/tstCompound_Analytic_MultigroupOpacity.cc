//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/test/tstCompound_Analytic_MultigroupOpacity.cc
 * \author Thomas M. Evans
 * \date   Tue Nov 13 17:24:12 2001
 * \brief  Compound_Analytic_MultigroupOpacity test.
 * \note   Copyright (C) 2019-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi_analytic_test.hh"
#include "dsxx/Assert.hh"
#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"
#include "cdi/CDI.hh"
#include "cdi/MultigroupOpacity.hh"
#include "cdi/OpacityCommon.hh"
#include "cdi_analytic/Analytic_Models.hh"
#include "cdi_analytic/Compound_Analytic_MultigroupOpacity.hh"
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

using rtt_cdi::CDI;
using rtt_cdi::MultigroupOpacity;
using rtt_cdi_analytic::Analytic_Opacity_Model;
using rtt_cdi_analytic::Compound_Analytic_MultigroupOpacity;
using rtt_cdi_analytic::Constant_Analytic_Opacity_Model;
using rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model;
using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void multigroup_test(rtt_dsxx::UnitTest &ut) {
  // group structure
  vector<double> groups = {0.05, 0.5, 5.0, 50.0};
  vector<shared_ptr<Analytic_Opacity_Model>> models(3);

  // make a Marshak (user-defined) model for the first group
  models[0] = std::make_shared<rtt_cdi_analytic_test::Marshak_Model>(100.0);

  // make a Polynomial model for the second group
  models[1] =
      std::make_shared<rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model>(1.5, 0.0, 0.0, 0.0);

  // make a Constant model for the third group
  models[2] = std::make_shared<rtt_cdi_analytic::Constant_Analytic_Opacity_Model>(3.0);

  // make an analytic multigroup opacity object for absorption
  Compound_Analytic_MultigroupOpacity const opacity(groups, models, rtt_cdi::ABSORPTION);

  // check the interface to multigroup opacity
  {
    string const desc = "Compound Multigroup Absorption";

    FAIL_IF(opacity.data_in_tabular_form());
    FAIL_IF(opacity.getReactionType() != rtt_cdi::ABSORPTION);
    FAIL_IF(opacity.getModelType() != rtt_cdi::ANALYTIC);
    FAIL_IF(opacity.getNumTemperatures() != 0);
    FAIL_IF(opacity.getNumDensities() != 0);
    FAIL_IF(opacity.getTemperatureGrid() != vector<double>());
    FAIL_IF(opacity.getDensityGrid() != vector<double>());
    FAIL_IF(opacity.getNumGroups() != 3);
    FAIL_IF(opacity.getNumGroupBoundaries() != 4);
    FAIL_IF(opacity.getEnergyPolicyDescriptor() != "mg");
    FAIL_IF(opacity.getDataDescriptor() != desc);
    FAIL_IF(opacity.getDataFilename() != string());
    FAIL_IF(opacity.getOpacityModelType() != rtt_cdi::ANALYTIC_TYPE);
  }
  {
    Compound_Analytic_MultigroupOpacity const anal_opacity(groups, models, rtt_cdi::SCATTERING);
    FAIL_IF(anal_opacity.getDataDescriptor() != "Compound Multigroup Scattering");
  }
  {
    Compound_Analytic_MultigroupOpacity const anal_opacity(groups, models, rtt_cdi::TOTAL);
    FAIL_IF(anal_opacity.getDataDescriptor() != "Compound Multigroup Total");
  }

  // check the group structure
  vector<double> mg_groups = opacity.getGroupBoundaries();

  if (soft_equiv(mg_groups.begin(), mg_groups.end(), groups.begin(), groups.end())) {
    PASSMSG("Group boundaries match.");
  } else {
    FAILMSG("Group boundaries do not match.");
  }

  // >>> get opacities

  // scalar density and temperature
  vector<double> sigma = opacity.getOpacity(2.0, 3.0);
  vector<double> ref(3, 0.0);
  {
    ref[0] = 100.0 / 8.0;
    ref[1] = 1.5;
    ref[2] = 3.0;
  }
  if (soft_equiv(sigma.begin(), sigma.end(), ref.begin(), ref.end())) {
    ostringstream message;
    message << "Analytic multigroup opacities are correct for scalar temperature and scalar "
            << "density.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Analytic multigroup opacities are NOT correct for scalar temperature and scalar "
            << "density.";
    FAILMSG(message.str());
  }

  // scalar density/temperature + vector density/temperature
  vector<double> const data_field(3, 2.0);
  vector<vector<double>> sigma_total = opacity.getOpacity(data_field, 3.0);
  vector<vector<double>> sig_rho = opacity.getOpacity(2.0, data_field);

  for (int i = 0; i < 3; i++) {
    vector<double> &test = sigma_total[i];

    if (soft_equiv(test.begin(), test.end(), ref.begin(), ref.end())) {
      ostringstream message;
      message << "Analytic multigroup opacities are correct for temperature field and scalar "
              << "density.";
      PASSMSG(message.str());
    } else {
      ostringstream message;
      message << "Analytic multigroup opacities are NOT correct for temperature field and scalar "
              << "density.";
      FAILMSG(message.str());
    }

    test = sig_rho[i];

    if (soft_equiv(test.begin(), test.end(), ref.begin(), ref.end())) {
      ostringstream message;
      message << "Analytic multigroup opacities are correct for density field and scalar "
              << "temperature.";
      PASSMSG(message.str());
    } else {
      ostringstream message;
      message << "Analytic multigroup opacities are NOT correct for density field and scalar "
              << "temperature.";
      FAILMSG(message.str());
    }
  }

  // Test the get_Analytic_Model() member function.
  {
    shared_ptr<Analytic_Opacity_Model const> const my_mg_opacity_model =
        opacity.get_Analytic_Model(1);
    shared_ptr<Analytic_Opacity_Model const> const expected_model(models[0]);

    if (expected_model == my_mg_opacity_model)
      PASSMSG("get_Analytic_Model() returned the expected MG Opacity model.");
    else
      FAILMSG("get_Analytic_Model() did not return the expected MG Opacity model.");
  }

  return;
}

//------------------------------------------------------------------------------------------------//

void test_CDI(rtt_dsxx::UnitTest &ut) {
  // group structure
  vector<double> groups = {0.05, 0.5, 5.0, 50.0};
  vector<shared_ptr<Analytic_Opacity_Model>> models(3);

  // make a Marshak (user-defined) model for the first group
  models[0] = std::make_shared<rtt_cdi_analytic_test::Marshak_Model>(100.0);

  // make a Polynomial model for the second group
  models[1] =
      std::make_shared<rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model>(1.5, 0.0, 0.0, 0.0);

  // make a Constant model for the third group
  models[2] = std::make_shared<rtt_cdi_analytic::Constant_Analytic_Opacity_Model>(3.0);

  // make an analytic multigroup opacity object for absorption
  shared_ptr<const MultigroupOpacity> const mg(
      new Compound_Analytic_MultigroupOpacity(groups, models, rtt_cdi::ABSORPTION));

  // make a CDI object
  CDI cdi;

  // set the multigroup opacity
  cdi.setMultigroupOpacity(mg);

  // check the energy groups from CDI
  vector<double> mg_groups = CDI::getFrequencyGroupBoundaries();

  if (soft_equiv(mg_groups.begin(), mg_groups.end(), groups.begin(), groups.end())) {
    PASSMSG("CDI Group boundaries match.");
  } else {
    FAILMSG("CDI Group boundaries do not match.");
  }

  // do a quick access test for getOpacity

  // scalar density and temperature
  vector<double> sigma = cdi.mg(rtt_cdi::ANALYTIC, rtt_cdi::ABSORPTION)->getOpacity(2.0, 3.0);
  vector<double> ref(3, 0.0);
  {
    ref[0] = 100.0 / 8.0;
    ref[1] = 1.5;
    ref[2] = 3.0;
  }
  if (soft_equiv(sigma.begin(), sigma.end(), ref.begin(), ref.end())) {
    PASSMSG("CDI multigroup opacities are correct for scalar temperature and scalar density.");
  } else {
    FAILMSG("CDI multigroup opacities are NOT correct for scalar temperature and scalar density.");
  }
}

//------------------------------------------------------------------------------------------------//

void packing_test(rtt_dsxx::UnitTest &ut) {

  // group structure
  vector<double> groups = {0.05, 0.5, 5.0, 50.0};
  vector<shared_ptr<Analytic_Opacity_Model>> models(3);

  // make a Polynomial model for the first group
  models[0] =
      std::make_shared<rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model>(0.0, 100.0, -3.0, 0.0);

  // make a Polynomial model for the second group
  models[1] =
      std::make_shared<rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model>(1.5, 0.0, 0.0, 0.0);

  // make a Constant model for the third group
  models[2] = std::make_shared<rtt_cdi_analytic::Constant_Analytic_Opacity_Model>(3.0);

  // make an analytic multigroup opacity object for absorption
  shared_ptr<const MultigroupOpacity> const mg(
      new Compound_Analytic_MultigroupOpacity(groups, models, rtt_cdi::ABSORPTION));

  // pack it
  vector<char> const packed = mg->pack();

  // now unpack it
  Compound_Analytic_MultigroupOpacity const opacity(packed);

  // now check it

  // check the interface to multigroup opacity
  {
    string const desc = "Compound Multigroup Absorption";

    FAIL_IF(opacity.data_in_tabular_form());
    FAIL_IF(opacity.getReactionType() != rtt_cdi::ABSORPTION);
    FAIL_IF(opacity.getModelType() != rtt_cdi::ANALYTIC);
    FAIL_IF(opacity.getNumTemperatures() != 0);
    FAIL_IF(opacity.getNumDensities() != 0);
    FAIL_IF(opacity.getTemperatureGrid() != vector<double>());
    FAIL_IF(opacity.getDensityGrid() != vector<double>());
    FAIL_IF(opacity.getNumGroups() != 3);
    FAIL_IF(opacity.getNumGroupBoundaries() != 4);
    FAIL_IF(opacity.getEnergyPolicyDescriptor() != "mg");
    FAIL_IF(opacity.getDataDescriptor() != desc);
    FAIL_IF(opacity.getDataFilename() != string());
  }

  // check the group structure
  vector<double> mg_groups = opacity.getGroupBoundaries();

  if (soft_equiv(mg_groups.begin(), mg_groups.end(), groups.begin(), groups.end())) {
    PASSMSG("Group boundaries for unpacked MG opacity match.");
  } else {
    FAILMSG("Group boundaries for unpacked MG do not match.");
  }

  // >>> get opacities

  // scalar density and temperature
  vector<double> sigma = opacity.getOpacity(2.0, 3.0);
  vector<double> ref(3, 0.0);
  {
    ref[0] = 100.0 / 8.0;
    ref[1] = 1.5;
    ref[2] = 3.0;
  }
  if (soft_equiv(sigma.begin(), sigma.end(), ref.begin(), ref.end())) {
    ostringstream message;
    message << "Analytic multigroup opacities for unpacked MG opacity are correct for scalar "
            << "temperature and scalar density.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Analytic multigroup opacities for unpacked MG opacity are NOT correct for scalar "
            << "temperature and scalar density.";
    FAILMSG(message.str());
  }

  // make sure we catch an assertion showing that we cannot unpack an
  // unregistered opacity

  // make a Marshak (user-defined) model for the first group
  models[0] = std::make_shared<rtt_cdi_analytic_test::Marshak_Model>(100.0);

  // make a Polynomial model for the second group
  models[1] =
      std::make_shared<rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model>(1.5, 0.0, 0.0, 0.0);

  // make a Constant model for the third group
  models[2] = std::make_shared<rtt_cdi_analytic::Constant_Analytic_Opacity_Model>(3.0);

  // make an analytic multigroup opacity object for absorption
  shared_ptr<const MultigroupOpacity> const mg2(
      new Compound_Analytic_MultigroupOpacity(groups, models, rtt_cdi::ABSORPTION));

  vector<char> const packed2 = mg2->pack();

  // we should catch an assertion when unpacking this because the Marshak_Model is not registered in
  // rtt_cdi::Opacity_Models
  bool caught = false;
  try {
    Compound_Analytic_MultigroupOpacity const nmg(packed2);
  } catch (const rtt_dsxx::assertion &ass) {
    caught = true;
    ostringstream message;
    message << "Caught the following assertion, " << ass.what();
    PASSMSG(message.str());
  }
  if (!caught) {
    FAILMSG("Failed to catch unregistered analyic model assertion");
  }
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    multigroup_test(ut);
    test_CDI(ut);
    packing_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstCompound_Analytic_MultigroupOpacity.cc
//------------------------------------------------------------------------------------------------//
