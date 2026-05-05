//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/test/tDummyOpacity.cc
 * \author Thomas M. Evans
 * \brief  GrayOpacity and Multigroup opacity test.
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "DummyGrayOpacity.hh"
#include "DummyMultigroupOpacity.hh"
#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "cdi/GrayOpacity.hh"
#include "cdi/MultigroupOpacity.hh"
#include "cdi/OpacityCommon.hh"
#include <cstddef>
#include <memory>
#include <sstream>
#include <vector>

namespace {

using namespace std;

using rtt_cdi::GrayOpacity;
using rtt_cdi::MultigroupOpacity;
using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void simple_tests(rtt_dsxx::UnitTest &ut) {
  // make shared_ptrs to gray and multigroup opacities
  std::shared_ptr<GrayOpacity> gray;
  std::shared_ptr<MultigroupOpacity> mg;

  // Assign and check gray opacity
  std::shared_ptr<rtt_cdi_test::DummyGrayOpacity> gray_total;
  std::shared_ptr<rtt_cdi_test::DummyGrayOpacity> gray_abs;
  gray_total = std::make_shared<rtt_cdi_test::DummyGrayOpacity>();
  gray_abs = std::make_shared<rtt_cdi_test::DummyGrayOpacity>(rtt_cdi::ABSORPTION);

  // check gray opacity for total opacities
  {
    gray = gray_total;

    FAIL_IF_NOT(gray->data_in_tabular_form());
    FAIL_IF_NOT(gray->getReactionType() == rtt_cdi::TOTAL);
    FAIL_IF_NOT(gray->getModelType() == rtt_cdi::ANALYTIC);
    FAIL_IF_NOT(gray->getOpacityModelType() == rtt_cdi::DUMMY_TYPE);

    vector<double> T_grid = {1.0, 2.0, 3.0};
    vector<double> rho_grid = {0.1, 0.2};

    std::vector<double> const T(gray->getTemperatureGrid());
    if (soft_equiv(T.begin(), T.end(), T_grid.begin(), T_grid.end()))
      PASSMSG("Gray temperature grid correct.");
    else
      FAILMSG("Gray temperature grid incorrect.");

    std::vector<double> const rho(gray->getDensityGrid());
    if (soft_equiv(rho.begin(), rho.end(), rho_grid.begin(), rho_grid.end()))
      PASSMSG("Gray density grid correct.");
    else
      FAILMSG("Gray density grid incorrect.");
  }

  // now reassign and check gray opacity for absorption opacities
  {
    gray = gray_abs;
    FAIL_IF(gray->getReactionType() != rtt_cdi::ABSORPTION);
  }

  // Assign and check multigroup opacity
  std::shared_ptr<rtt_cdi_test::DummyMultigroupOpacity> mg_total;
  std::shared_ptr<rtt_cdi_test::DummyMultigroupOpacity> mg_abs;
  mg_total = std::make_shared<rtt_cdi_test::DummyMultigroupOpacity>();
  mg_abs = std::make_shared<rtt_cdi_test::DummyMultigroupOpacity>(rtt_cdi::ABSORPTION);

  // check multigroup total opacities
  {
    mg = mg_total;

    FAIL_IF(!mg->data_in_tabular_form());
    FAIL_IF(mg->getReactionType() != rtt_cdi::TOTAL);
    FAIL_IF(mg->getModelType() != rtt_cdi::ANALYTIC);
    FAIL_IF(mg->getOpacityModelType() != rtt_cdi::DUMMY_TYPE);

    vector<double> T_grid = {1.0, 2.0, 3.0};
    vector<double> rho_grid = {0.1, 0.2};
    vector<double> e_groups = {0.05, 0.5, 5.0, 50.0};

    std::vector<double> const T(mg->getTemperatureGrid());
    if (soft_equiv(T.begin(), T.end(), T_grid.begin(), T_grid.end()))
      PASSMSG("Multigroup temperature grid correct.");
    else
      FAILMSG("Multigroup temperature grid incorrect.");

    std::vector<double> const rho(mg->getDensityGrid());
    if (soft_equiv(rho.begin(), rho.end(), rho_grid.begin(), rho_grid.end()))
      PASSMSG("Multigroup density grid correct.");
    else
      FAILMSG("Multigroup density grid incorrect.");

    std::vector<double> const bounds(mg->getGroupBoundaries());
    if (soft_equiv(bounds.begin(), bounds.end(), e_groups.begin(), e_groups.end()))
      PASSMSG("Multigroup energy boundaries correct.");
    else
      FAILMSG("Multigroup energy boundaries incorrect.");

    FAIL_IF_NOT(mg->getNumTemperatures() == 3);
    FAIL_IF_NOT(mg->getNumDensities() == 2);
    FAIL_IF_NOT(mg->getNumGroupBoundaries() == 4);
  }

  // now reassign and check multigroup opacities for absorption
  {
    mg = mg_abs;
    FAIL_IF_NOT(mg->getReactionType() == rtt_cdi::ABSORPTION);
  }
  return;
}

//------------------------------------------------------------------------------------------------//

void gray_opacity_test(rtt_dsxx::UnitTest &ut) {
  // ---------------------------- //
  // Create a GrayOpacity object. //
  // ---------------------------- //

  std::shared_ptr<GrayOpacity> spDGO;

  if ((spDGO = std::make_shared<rtt_cdi_test::DummyGrayOpacity>()), spDGO)
    PASSMSG("shared_ptr to new GrayOpacity object created.");
  else
    FAILMSG("Unable to create a shared_ptr to new GrayOpacity object.");

  // ------------------------ //
  // Dummy Gray Opacity Tests //
  // ------------------------ //

  double const temperature = 0.1;                               // keV
  double density = 27.0;                                        // g/cm^3
  double const tabulatedGrayOpacity = temperature + density / 1000.0; // cm^2/g

  double const opacity = spDGO->getOpacity(temperature, density);

  if (soft_equiv(opacity, tabulatedGrayOpacity)) {
    ostringstream message;
    message << spDGO->getDataDescriptor() << " getOpacity computation was good.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spDGO->getDataDescriptor() << " getOpacity value is out of spec.";
    FAILMSG(message.str());
  }

  // try using a vector of temps.

  std::vector<double> vtemperature(2);
  vtemperature[0] = 0.5; // keV
  vtemperature[1] = 0.7; // keV
  density = 0.35;        // g/cm^3
  std::vector<double> vRefOpacity(vtemperature.size());
  for (size_t i = 0; i < vtemperature.size(); ++i)
    vRefOpacity[i] = vtemperature[i] + density / 1000;

  std::vector<double> vOpacity = spDGO->getOpacity(vtemperature, density);

  if (soft_equiv(vOpacity.begin(), vOpacity.end(), vRefOpacity.begin(), vRefOpacity.end())) {
    ostringstream message;
    message << spDGO->getDataDescriptor()
            << " getOpacity computation was good for a vector of temps.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spDGO->getDataDescriptor()
            << " getOpacity value is out of spec. for a vector of temps.";
    FAILMSG(message.str());
  }

  // try using a vector of densities.

  std::vector<double> vdensity = {0.5, 1.0, 3.3, 5.0, 27.0};

  vRefOpacity.resize(vdensity.size());
  for (size_t i = 0; i < vdensity.size(); ++i)
    vRefOpacity[i] = temperature + vdensity[i] / 1000;

  vOpacity = spDGO->getOpacity(temperature, vdensity);

  if (soft_equiv(vOpacity.begin(), vOpacity.end(), vRefOpacity.begin(), vRefOpacity.end())) {
    ostringstream message;
    message << spDGO->getDataDescriptor() << " getOpacity computation was good for a vector of "
            << "densities.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spDGO->getDataDescriptor() << " getOpacity value is out of spec. for a vector of "
            << "densities.";
    FAILMSG(message.str());
  }
}

//------------------------------------------------------------------------------------------------//

void multigroup_opacity_test(rtt_dsxx::UnitTest &ut) {
  // ----------------------------------------- //
  // Create a Dummy Multigroup Opacity object. //
  // ----------------------------------------- //

  std::shared_ptr<MultigroupOpacity> spDmgO;

  if ((spDmgO = std::make_shared<rtt_cdi_test::DummyMultigroupOpacity>()), spDmgO) {
    PASSMSG("shared_ptr to new MultigroupOpacity object created.");
  }

  // --------------- //
  // MG Opacity test //
  // --------------- //

  // Setup the test point.
  double const temperature = 0.01; // keV
  double const density = 2.0;      // g/cm^3

  // declare vector temps and densities
  vector<double> vtemperature = {0.5, 0.7}; // keV
  vector<double> vdensity = {1.5, 2.0};     // g/cc

  // The dummy opacity object should have 3 groups.  Check it.
  size_t const ng = spDmgO->getNumGroupBoundaries() - 1;
  if (ng == 3) {
    PASSMSG("Correct number of groups found for MultigroupOpacity object.");
  } else {
    FAILMSG("Wrong number of groups found for MultigroupOpacity object.");
  }

  const std::vector<double> energyBoundaries = spDmgO->getGroupBoundaries();

  // Create a container that hold all the MG opacities for a specified temperature and density.
  // Fill this container with the values that DummyMultigroupOpacity should contain.
  std::vector<double> tabulatedMGOpacity(ng);
  for (size_t ig = 0; ig < ng; ++ig)
    tabulatedMGOpacity[ig] =
        2 * (temperature + density / 1000) / (energyBoundaries[ig] + energyBoundaries[ig + 1]);

  // Use the getOpacity accessor to obtain the MG opacities for a specified temperature and density.
  std::vector<double> mgOpacity = spDmgO->getOpacity(temperature, density);

  // Make sure the accessor values match the expected values.
  if (soft_equiv(mgOpacity.begin(), mgOpacity.end(), tabulatedMGOpacity.begin(),
                 tabulatedMGOpacity.end())) {
    PASSMSG(spDmgO->getDataDescriptor() + " getOpacity computation was good.");
  } else {
    FAILMSG(spDmgO->getDataDescriptor() + " getOpacity value is out of spec.");
  }

  // Repeat with a vector of temps.

  // Reference values.

  // The opacity container is a `vector<vector<double>>`.  Each nested vector contains all of the
  // group opacity values for a single temperature.

  // a MG opacity set for a single temperature, density combination can be extracted from this
  // container by using the following type of assignment.
  // * `std::vector< double > vec1 = vRefMgOpacity[0];`

  // the size of this vector is the number of temperatures, ***not*** the number of groups!
  std::vector<std::vector<double>> vRefMgOpacity(2);
  for (size_t it = 0; it < vtemperature.size(); ++it) {
    vRefMgOpacity[it].resize(ng);
    for (size_t ig = 0; ig < ng; ++ig)
      vRefMgOpacity[it][ig] = 2.0 * (vtemperature[it] + density / 1000.0) /
                              (energyBoundaries[ig] + energyBoundaries[ig + 1]);
  }

  // Retrieve the same set of opacity values via the getOpacity() accessor.
  std::vector<std::vector<double>> const vMgOpacity = spDmgO->getOpacity(vtemperature, density);

  // Compare the results.
  if (soft_equiv(vMgOpacity, vRefMgOpacity)) {
    ostringstream message;
    message << spDmgO->getDataDescriptor()
            << " getOpacity computation was good for a vector of  temperatures.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spDmgO->getDataDescriptor()
            << " getOpacity value is out of spec for a vector of temperatures.";
    FAILMSG(message.str());
  }

  // STL-like accessor (MG opacities)

  // We have added STL-like getOpacity functions to DummyMultigroupOpacity, these are not available
  // through the rtt_cdi::MultigroupOpacity base class so we test them as a DummyMultigroupOpacity.
  // This demonstrates that one could make an opacity class that contains extra functionality. Of
  // course this functionality is not available through CDI.

  std::shared_ptr<rtt_cdi_test::DummyMultigroupOpacity> spDumMgOp;
  if ((spDumMgOp = std::make_shared<rtt_cdi_test::DummyMultigroupOpacity>()), spDumMgOp) {
    PASSMSG("shared_ptr to new DummyMultigroupOpacity object created.");
  } else {
    FAILMSG("Unable to create a shared_ptr to a new DummyMultigroupOpacity object.");
  }

  // The STL-like accessors only work with 1-D containers.
  vector<double> vOpacity;
  vector<double> vRefOpacity;

  vOpacity.resize(vtemperature.size() * ng);
  vRefOpacity.resize(vtemperature.size() * ng);

  // Reference Values
  for (size_t it = 0; it < vtemperature.size(); ++it)
    for (size_t ig = 0; ig < ng; ++ig)
      vRefOpacity[it * ng + ig] = 2.0 * (vtemperature[it] + vdensity[it] / 1000.0) /
                                  (energyBoundaries[ig] + energyBoundaries[ig + 1]);

  // Obtain values using getOpacity() accessor.
  auto it = spDumMgOp->getOpacity(vtemperature.begin(), vtemperature.end(), vdensity.begin(),
                                  vdensity.end(), vOpacity.begin());
  FAIL_IF_NOT(it == vOpacity.end());

  // Compare the results:
  if (soft_equiv(vOpacity.begin(), vOpacity.end(), vRefOpacity.begin(), vRefOpacity.end())) {
    ostringstream message;
    message << spDumMgOp->getDataDescriptor() << " STL getOpacity() computation was good for a\n"
            << " vector of temps. and a vector of densities.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spDumMgOp->getDataDescriptor() << " STL getOpacity() value is out of spec. for a\n"
            << " vector of temps. and a vector of densities.";
    FAILMSG(message.str());
  }
  return;
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    simple_tests(ut);
    gray_opacity_test(ut);
    multigroup_opacity_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tDummyOpacity.cc
//------------------------------------------------------------------------------------------------//
