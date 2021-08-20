//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/test/tCDI.cc
 * \author Thomas M. Evans
 * \date   Tue Oct  9 15:52:01 2001
 * \brief  CDI test executable.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "DummyEoS.hh"
#include "DummyGrayOpacity.hh"
#include "DummyMultigroupOpacity.hh"
#include "cdi/CDI.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>

using namespace std;

using rtt_cdi::CDI;
using rtt_cdi::EoS;
using rtt_cdi::GrayOpacity;
using rtt_cdi::MultigroupOpacity;
using rtt_cdi_test::DummyEoS;
using rtt_cdi_test::DummyGrayOpacity;
using rtt_cdi_test::DummyMultigroupOpacity;
using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void check_CDI(rtt_dsxx::UnitTest &ut, CDI const &cdi) {
  // check cdi, note that the different combinations of rtt_cdi::Model and rtt_cdi::Reaction will
  // yield the same results because DummyOpacity, DummyMultigroupOpacity, and DummyEoS all yield the
  // same stuff.  These have all been tested in tDummyOpacity and tDummyEoS.  Here we just check the
  // types

  // Work around an MSVC ICE in MSC_VER 1925 (VS 2019, 16.5.0 preview 2)
#if !defined(_MSC_VER) || _MSC_VER < 1900 || _MSC_VER > 1925
  // check for gray
  if (typeid(*cdi.gray(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION)) == typeid(DummyGrayOpacity))
    PASSMSG("CDI gray() returned the correct type!");
  else
    FAILMSG("CDI gray() did not return the correct type!");

  if (typeid(*cdi.gray(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)) == typeid(DummyGrayOpacity))
    PASSMSG("CDI gray() returned the correct type!");
  else
    FAILMSG("CDI gray() did not return the correct type!");

  // check for multigroup
  if (typeid(*cdi.mg(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION)) == typeid(DummyMultigroupOpacity))
    PASSMSG("CDI mg() returned the correct type!");
  else
    FAILMSG("CDI mg() did not return the correct type!");

  if (typeid(*cdi.mg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)) == typeid(DummyMultigroupOpacity))
    PASSMSG("CDI mg() returned the correct type!");
  else
    FAILMSG("CDI mg() did not return the correct type!");

  if (typeid(*cdi.eos()) == typeid(DummyEoS))
    PASSMSG("CDI eos() returned the correct type!");
  else
    FAILMSG("CDI eos() did not return the correct type!");
#endif

  // gray test case: Find the value of opacity at T=0.35 keV and rho = 27.2 g/cm^3.  For
  // DummyGrayOpacity the value should be .35272 cm^2/g.

  double temp = 0.35;                       // keV
  double dens = 27.2;                       // g/cm^3
  double refOpacity = temp + dens / 1000.0; // cm^2/g

  double opacity = cdi.gray(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION)->getOpacity(temp, dens);

  if (soft_equiv(opacity, refOpacity))
    PASSMSG("CDI.gray()->getOpacity is ok.");
  else
    FAILMSG("CDI.gray()->getOpacity is not ok.");

  // mg test case: Find the mg opacities at T=0.35 keV and rho = 27.2 g/cm^3.  For
  // DummyMultigroupOpacity the values should be { }.  Three groups
  size_t ng = 3;

  // The energy groups in DummyMultigroupOpacity are hardwired to be { 0.05, 0.5, 5.0, 50.0 } keV.
  std::vector<double> energyBoundary(ng + 1);
  energyBoundary[0] = 0.05;
  energyBoundary[1] = 0.5;
  energyBoundary[2] = 5.0;
  energyBoundary[3] = 50.0;

  if (cdi.mg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)->getNumGroups() == ng)
    PASSMSG("CDI.mg()->getNumGroups() is ok.");
  else
    FAILMSG("CDI.mg()->getNumGroups() is not ok.");

  std::vector<double> vRefOpacity(3);
  for (size_t group = 0; group < ng; ++group)
    vRefOpacity[group] =
        2.0 * (temp + dens / 1000.0) / (energyBoundary[group] + energyBoundary[group + 1]);

  std::vector<double> vOpacity(ng);
  vOpacity = cdi.mg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)->getOpacity(temp, dens);

  if (soft_equiv(vOpacity.begin(), vOpacity.end(), vRefOpacity.begin(), vRefOpacity.end()))
    PASSMSG("CDI.mg()->getOpacity(T,rho) is ok.");
  else
    FAILMSG("CDI.mg()->getOpacity(T,rho) is not ok.");

  // Test the EoS plug-in component of this CDI.

  double refCve = temp + dens / 1000.0;

  double Cve = cdi.eos()->getElectronHeatCapacity(temp, dens);

  if (soft_equiv(Cve, refCve)) {
    ostringstream message;
    message << "CDI.eos()->getElectronHeatCapacity( dbl, dbl )\n\t"
            << "returned a value that matched the reference value.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "CDI.eos()->getElectronHeatCapacity( dbl, dbl )\n\t"
            << "returned a value that was out of spec.";
    FAILMSG(message.str());
  }

  double U = 3.1415;
  double Ti = cdi.eos()->getIonTemperature(dens, U, 1.0);
  double Te = cdi.eos()->getElectronTemperature(dens, U, 1.0);

  if (soft_equiv(Ti, 1.772427713617681)) {
    ostringstream message;
    message << "CDI.eos()->getIonTemperature( dbl, dbl )\n\t"
            << "returned a value that matched the reference value.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "CDI.eos()->getIonTemperature( dbl, dbl )\n\t"
            << "returned a value that was out of spec.";
    FAILMSG(message.str());
  }
  if (soft_equiv(Te, 31.415)) {
    ostringstream message;
    message << "CDI.eos()->getElectronTemperature( dbl, dbl )\n\t"
            << "returned a value that matched the reference value.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "CDI.eos()->getElectronTemperature( dbl, dbl )\n\t"
            << "returned a value that was out of spec.";
    FAILMSG(message.str());
  }

  // Print the material ID string

  // Test Rosseland integration with MG opacities.
  {
    // integrate on the interval [0.5, 5.0] keV for T=1.0 keV.
    double int_total1 = rtt_cdi::integrateRosselandSpectrum(0.5, 5.0, 1.0);
    // group 1 has the same energy range.
    double int_total2 = cdi.integrateRosselandSpectrum(2, 1.0);
    if (soft_equiv(int_total1, int_total2, 1.0e-7)) {
      ostringstream message;
      message.precision(10);
      message << "Calculated a total normalized Rosseland integral of " << setw(12)
              << setiosflags(ios::fixed) << int_total2;
      PASSMSG(message.str());
    } else {
      ostringstream message;
      message.precision(10);
      message << "Calculated a total normalized Rosseland integral of " << setw(12)
              << setiosflags(ios::fixed) << int_total2 << " instead of " << setw(12)
              << setiosflags(ios::fixed) << int_total1 << "." << endl;
      FAILMSG(message.str());
    }
  }
  return;
}

//------------------------------------------------------------------------------------------------//
void test_CDI(rtt_dsxx::UnitTest &ut) {
  // make shared_ptrs to opacity and EoS objects
  std::shared_ptr<const GrayOpacity> gray_planck_abs;
  std::shared_ptr<const GrayOpacity> gray_iso_scatter;
  std::shared_ptr<const MultigroupOpacity> mg_planck_abs;
  std::shared_ptr<const MultigroupOpacity> mg_iso_scatter;
  std::shared_ptr<const MultigroupOpacity> mg_diff_bound;
  std::shared_ptr<const EoS> eos;

  // assign to dummy state objects
  gray_planck_abs = std::make_shared<DummyGrayOpacity>(rtt_cdi::ABSORPTION, rtt_cdi::PLANCK);
  gray_iso_scatter = std::make_shared<DummyGrayOpacity>(rtt_cdi::SCATTERING, rtt_cdi::ISOTROPIC);
  mg_planck_abs = std::make_shared<DummyMultigroupOpacity>(rtt_cdi::ABSORPTION, rtt_cdi::PLANCK);
  mg_iso_scatter =
      std::make_shared<DummyMultigroupOpacity>(rtt_cdi::SCATTERING, rtt_cdi::ISOTROPIC);
  //multigroup with different boundaries
  mg_diff_bound =
      std::make_shared<DummyMultigroupOpacity>(rtt_cdi::SCATTERING, rtt_cdi::THOMSON, 6);

  // EOS
  eos = std::make_shared<DummyEoS>();

  // make a CDI, it should be empty
  std::string const matName("dummyMaterial");
  CDI cdi(matName);
  for (size_t i = 0; i < rtt_cdi::constants::num_Models; i++)
    for (size_t j = 0; j < rtt_cdi::constants::num_Reactions; j++) {
      // casting these is inherently dangerous, but because this is a test I won't go nuts
      auto const m = static_cast<rtt_cdi::Model>(i);
      auto const r = static_cast<rtt_cdi::Reaction>(j);

      if (cdi.isGrayOpacitySet(m, r))
        ITFAILS;
      if (cdi.isMultigroupOpacitySet(m, r))
        ITFAILS;
    }
  if (cdi.isEoSSet())
    ITFAILS;

  if (cdi.isEICouplingSet())
    ITFAILS;

  if (cdi.getMatID() == matName)
    PASSMSG("Good, the material identifier was set and fetched correctly.");
  else
    FAILMSG("Oh-ho, the material identifier was not set and fetched correctly.");

  // there should be no energy group boundaries set yet
  if (CDI::getFrequencyGroupBoundaries().empty())
    PASSMSG("Good, no frequency group boundaries defined yet.");
  else
    FAILMSG("Oh-oh, frequency boundaries are defined.");

  // now assign stuff to it
  cdi.setGrayOpacity(gray_planck_abs);
  cdi.setGrayOpacity(gray_iso_scatter);
  cdi.setMultigroupOpacity(mg_planck_abs);
  cdi.setMultigroupOpacity(mg_iso_scatter);
  cdi.setEoS(eos);

  // check the energy group boundaries
  {
    vector<double> b1 = CDI::getFrequencyGroupBoundaries();
    vector<double> b2 = cdi.mg(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION)->getGroupBoundaries();
    vector<double> b3 = cdi.mg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)->getGroupBoundaries();

    // these should all be equal
    bool test = true;

    if (!soft_equiv(b1.begin(), b1.end(), b2.begin(), b2.end()))
      test = false;

    if (!soft_equiv(b1.begin(), b1.end(), b3.begin(), b3.end()))
      test = false;

    if (test)
      PASSMSG("All multigroup data has consistent energy groups.");
    else
      FAILMSG("Multigroup data has inconsistent energy groups.");
  }

  // catch an exception when we try to assign a multigroup opacity to CDI that has a different
  // frequency group structure
  bool caught = false;
  try {
    cdi.setMultigroupOpacity(mg_diff_bound);
  } catch (const rtt_dsxx::assertion &error) {
    ostringstream message;
    message << "Good, we caught the following exception: \n" << error.what();
    PASSMSG(message.str());
    caught = true;
  }
  if (!caught) {
    ostringstream message;
    message << "Failed to catch an exception for setting a different "
            << "frequency group structure.";
    FAILMSG(message.str());
  }

  // make sure these are assigned
  if (cdi.isGrayOpacitySet(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION))
    PASSMSG("Gray planck absorption set!");
  else
    FAILMSG("Gray planck absorption not set!");

  if (cdi.isGrayOpacitySet(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING))
    PASSMSG("Gray isotropic scattering set!");
  else
    FAILMSG("Gray isotropic scattering not set!");

  if (cdi.isMultigroupOpacitySet(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION))
    PASSMSG("Multigroup planck (in-group) absorption set!");
  else
    FAILMSG("Multigroup planck (in-group) absorption not set!");

  if (cdi.isMultigroupOpacitySet(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING))
    PASSMSG("Multigroup isotropic scattering set!");
  else
    FAILMSG("Multigroup isotropic scattering not set!");

  if (cdi.isEoSSet())
    PASSMSG("EoS set!");
  else
    FAILMSG("EoS not set!");

  // catch some exceptions
  caught = false;
  try {
    cdi.setGrayOpacity(gray_planck_abs);
  } catch (const rtt_dsxx::assertion &error) {
    ostringstream message;
    message << "Good, we caught the following exception: \n" << error.what();
    PASSMSG(message.str());
    caught = true;
  }
  if (!caught) {
    FAILMSG("Failed to catch overset exception!");
  }

  caught = false;
  try {
    cdi.mg(rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION);
  } catch (const rtt_dsxx::assertion &error) {
    ostringstream message;
    message << "Good, we caught the following exception: \n" << error.what();
    PASSMSG(message.str());
    caught = true;
  }
  if (!caught)
    FAILMSG("Failed to catch an illegal access exception!");

  // check the cdi through a function call
  check_CDI(ut, cdi);

  // reset and make sure we are empty
  cdi.reset();
  if (!cdi.isGrayOpacitySet(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION))
    PASSMSG("Gray planck absorption unset!");
  else
    FAILMSG("Gray planck absorption is still set!");

  if (!cdi.isGrayOpacitySet(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING))
    PASSMSG("Gray isotropic scattering unset!");
  else
    FAILMSG("Gray isotropic scattering is still set!");

  if (!cdi.isMultigroupOpacitySet(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION))
    PASSMSG("Multigroup planck (in-group) absorption unset!");
  else
    FAILMSG("Multigroup planck (in-group) absorption is still set!");

  if (!cdi.isMultigroupOpacitySet(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING))
    PASSMSG("Multigroup isotropic scattering unset!");
  else
    FAILMSG("Multigroup isotropic scattering is still set!");

  if (!cdi.isEoSSet())
    PASSMSG("EoS unset!");
  else
    FAILMSG("EoS is still set!");

  caught = false;
  try {
    cdi.mg(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION);
  } catch (const rtt_dsxx::assertion &error) {
    ostringstream message;
    message << "Good, we caught the following exception: \n" << error.what();
    PASSMSG(message.str());
    caught = true;
  }
  if (!caught) {
    FAILMSG("Failed to catch an illegal access exception!");
  }

  // there should be no energy group boundaries set yet
  if (CDI::getFrequencyGroupBoundaries().empty())
    PASSMSG("Good, no frequency group boundaries defined after reset.");
  else
    FAILMSG("Oh-oh, frequency boundaries are defined after reset.");

  if (ut.numFails == 0) {
    PASSMSG("Fundamental CDI Operations are ok.");
    cout << endl;
  }
}

//------------------------------------------------------------------------------------------------//
void test_planck_integration(rtt_dsxx::UnitTest &ut) {
  // We have not defined any group structure yet; thus, the Insist will always fire if an
  // integration is requested over a non-existent group
  bool caught = false;
  try {
    rtt_cdi::CDI::integratePlanckSpectrum(1, 1.0);
  } catch (const rtt_dsxx::assertion &error) {
    ostringstream message;
    message << "Caught illegal Planck calculation exception: \n"
            << "\t" << error.what();
    PASSMSG(message.str());
    caught = true;
  }
  if (!caught) {
    FAILMSG("Did not catch an exception for calculating Planck integral.");
  }

  // check some planck integrals
  double int_total = rtt_cdi::integratePlanckSpectrum(0.0, 100.0, 1.0);
  if (soft_equiv(int_total, 1.0, 3.5e-10)) {
    ostringstream message;
    message.precision(10);
    message << "Calculated a total normalized Planck integral of " << setw(12)
            << setiosflags(ios::fixed) << int_total;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message.precision(10);
    message << "Calculated a total normalized Planck integral of " << setw(12)
            << setiosflags(ios::fixed) << int_total << " instead of 1.0.";
    FAILMSG(message.str());
  }

  double int_1 = rtt_cdi::integratePlanckSpectrum(0.0, 5.0, 10.0);
  if (soft_equiv(int_1, .00529316, 1.0e-6)) {
    ostringstream message;
    message.precision(10);
    message << "Calculated a normalized Planck integral of " << setw(12) << setiosflags(ios::fixed)
            << int_1;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message.precision(10);
    message << "Calculated a normalized Planck integral of " << setw(12) << setiosflags(ios::fixed)
            << int_1 << " instead of .00529316.";
    FAILMSG(message.str());
  }

  double int_2 = rtt_cdi::integratePlanckSpectrum(0.0, .50, 10.0);
  if (soft_equiv(int_2, 6.29674e-6, 1.0e-6)) {
    ostringstream message;
    message.precision(10);
    message << "Calculated a normalized Planck integral of " << setw(12) << setiosflags(ios::fixed)
            << int_2;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message.precision(10);
    message << "Calculated a normalized Planck integral of " << setw(12) << setiosflags(ios::fixed)
            << int_2 << " instead of 6.29674e-6.";
    FAILMSG(message.str());
  }

  double int_0 = rtt_cdi::integratePlanckSpectrum(0.0, 0.0, 10.0);
  if (!soft_equiv(int_0, 0.0, 1.e-6))
    ITFAILS;

  double int_range = rtt_cdi::integratePlanckSpectrum(0.1, 1.0, 1.0);
  if (!soft_equiv(int_range, 0.0345683, 3.e-5))
    ITFAILS;

  //
  // Extreme cases: test early return logic for integration routines.
  //
  double const tol = std::numeric_limits<double>::epsilon();

  // Extreme case 1. This should do the normal computation, but the result is zero.
  {
    double const integrate_range_cold = rtt_cdi::integratePlanckSpectrum(0.001, 1.0, 1.0e-30);
    if (!soft_equiv(integrate_range_cold, 0.0))
      ITFAILS;
  }
  // Extreme case 2. T < numeric_limits<double>::min() --> follow special logic and return zero.
  {
    double const integrate_range_cold = rtt_cdi::integratePlanckSpectrum(0.1, 1.0, 1.0e-308);
    if (!soft_equiv(integrate_range_cold, 0.0, tol))
      ITFAILS;
  }
  // Extreme case 3. This should do the normal computation, but the result is zero.
  {
    std::vector<double> const bounds = {1.0, 3.0, 30.0};
    // If T < sqrt(numeric_limits<double>::min()), the integration will fail.
    double const T_eval = 1.0e-300;
    std::vector<double> planck;
    CDI::integrate_Planckian_Spectrum(bounds, T_eval, planck);

    // make sure the alternate calling method returns the same values
    auto planck_alt = CDI::integrate_Planckian_Spectrum(bounds, T_eval);

    if (planck_alt != planck)
      ITFAILS;
    if (!soft_equiv(planck[0], 0.0, tol))
      ITFAILS;
  }
  // Extreme case 4a. T < numeric_limits<double>::min() --> catch the thrown exception
  caught = true;
  Remember(caught = false;);
  try {
    std::vector<double> const bounds = {1.0, 3.0, 30.0};
    double const T_eval = 1.0e-308;
    std::vector<double> planck;
    CDI::integrate_Planckian_Spectrum(bounds, T_eval, planck);
  } catch (const rtt_dsxx::assertion &error) {
    ostringstream message;
    message << "Good, we caught the following exception: \n" << error.what();
    PASSMSG(message.str());
    caught = true;
  }
  if (!caught) {
    ostringstream message;
    message << "Failed to catch an exception when passing a denorm temperature.";
  }

  // Extreme case 4b. bounds < numeric_limits<double>::min() --> return 1.0, 0.0
  {
    std::vector<double> const bounds = {0.0, 3.0, 30.0};
    double const T_eval = 1.0e-300;
    std::vector<double> planck;
    CDI::integrate_Planckian_Spectrum(bounds, T_eval, planck);
    // make sure the alternate calling method returns the same values
    auto planck_alt = CDI::integrate_Planckian_Spectrum(bounds, T_eval);
    if (planck_alt != planck)
      ITFAILS;
    if (!soft_equiv(planck[0], 1.0, tol))
      ITFAILS;
    if (!soft_equiv(planck[1], 0.0, tol))
      ITFAILS;
  }

  // Catch an illegal group assertion.
  CDI cdi;
  std::shared_ptr<const MultigroupOpacity> mg(
      new DummyMultigroupOpacity(rtt_cdi::SCATTERING, rtt_cdi::THOMSON));
  cdi.setMultigroupOpacity(mg);

  // Check the normalized planck integrals.
  if (CDI::getNumberFrequencyGroups() != 3)
    ITFAILS;

  double g1_integral = rtt_cdi::CDI::integratePlanckSpectrum(1, 1.0);
  double g2_integral = rtt_cdi::CDI::integratePlanckSpectrum(2, 1.0);
  double g3_integral = rtt_cdi::CDI::integratePlanckSpectrum(3, 1.0);
  double g_total = rtt_cdi::CDI::integratePlanckSpectrum(1.0);

  if (soft_equiv(g1_integral, 0.00528686, 1.e-6))
    PASSMSG("Group 1 integral within tolerance.");
  else
    FAILMSG("Group 1 integral fails tolerance.");

  if (soft_equiv(g2_integral, 0.74924, 1.e-6))
    PASSMSG("Group 2 integral within tolerance.");
  else
    FAILMSG("Group 2 integral fails tolerance.");

  if (soft_equiv(g3_integral, 0.245467, 1.e-6))
    PASSMSG("Group 3 integral within tolerance.");
  else
    FAILMSG("Group 3 integral fails tolerance.");

  if (soft_equiv(g_total, 0.999994, 1.e-6))
    PASSMSG("Total integral over groups within tolerance.");
  else
    FAILMSG("Total integral over groups fails tolerance.");

  // Test that a zero temperature returns a zero.
  if (soft_equiv(rtt_cdi::integratePlanckSpectrum(0.0, 100.0, 0.0), 0.0))
    PASSMSG("Planck integral from hnu=0 to 100 at T=0 is zero: good!");
  else
    FAILMSG("Planck integral from hnu=0 to 100 at T=0 is not zero: BAD!");

  if (soft_equiv(rtt_cdi::CDI::integratePlanckSpectrum(1, 0.0), 0.0))
    PASSMSG("Planck integral of group 1 at T=0 is zero: good!");
  else
    FAILMSG("Planck integral of group 1 at T=0 is not zero: BAD!");

  if (soft_equiv(rtt_cdi::CDI::integratePlanckSpectrum(2, 0.0), 0.0))
    PASSMSG("Planck integral of group 2 at T=0 is zero: good!");
  else
    FAILMSG("Planck integral of group 2 at T=0 is not zero: BAD!");

  if (soft_equiv(rtt_cdi::CDI::integratePlanckSpectrum(3, 0.0), 0.0))
    PASSMSG("Planck integral of group 3 at T=0 is zero: good!");
  else
    FAILMSG("Planck integral of group 3 at T=0 is not zero: BAD!");

  if (soft_equiv(rtt_cdi::CDI::integratePlanckSpectrum(0.0), 0.0))
    PASSMSG("Planck integral over all groups at T=0 is zero: good!");
  else
    FAILMSG("Planck integral over all groups at T=0 is not zero: BAD!");

  if (ut.numFails == 0) {
    PASSMSG("All Planckian integral tests ok.");
    cout << endl;
  }

  // Compare the integration over all groups to the individual groups:

  std::vector<double> group_bounds(CDI::getFrequencyGroupBoundaries());
  if (group_bounds.size() != 4)
    ITFAILS;

  std::vector<double> planck;

  CDI::integrate_Planckian_Spectrum(group_bounds, 1.0, planck);
  // make sure the alternate calling method returns the same values
  auto planck_alt = CDI::integrate_Planckian_Spectrum(group_bounds, 1.0);

  for (int group_index = 1; group_index <= 3; ++group_index) {

    double planck_g = rtt_cdi::CDI::integratePlanckSpectrum(group_index, 1.0);

    if (!soft_equiv(planck[group_index - 1], planck_g))
      ITFAILS;
  }
  if (planck_alt != planck)
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("Group-wise and Full spectrum Planckian and Rosseland integrals "
            "match.");
  else
    FAILMSG("Group-wise and Full spectrum Planckian and Rosseland integrals do "
            "not match.");
  return;
}

//------------------------------------------------------------------------------------------------//
void test_rosseland_integration(rtt_dsxx::UnitTest &ut) {
  // Only report this as a failure if 1) the error was not caught AND 2) the Require macro is
  // available.
  if (ut.dbcRequire() && !ut.dbcNothrow()) {
    bool caught = false;
    try {
      CDI::integrateRosselandSpectrum(0, 1.0);
    } catch (const rtt_dsxx::assertion & /*error*/) {
      PASSMSG("Caught illegal Rosseland calculation exception:");
      caught = true;
    }
    if (!caught)
      FAILMSG(string("Did not catch an exception for calculating ") + "Rosseland integral.");

    // catch our assertion
    caught = false;
    double P, R;
    try {
      CDI::integrate_Rosseland_Planckian_Spectrum(0, 1.0, P, R);
    } catch (const rtt_dsxx::assertion & /*error*/) {
      PASSMSG(string("Caught illegal Rosseland and Planckian ") + "calculation exception:");
      caught = true;
    }
    if (!caught)
      FAILMSG(string("Did not catch an exception for calculating ") +
              "Rosseland and Planckian integral.");
  }

  // check some planck integrals
  double int_total = rtt_cdi::integrateRosselandSpectrum(0, 100.0, 1.0);
  if (soft_equiv(int_total, 1.0, 1.0e-7)) {
    ostringstream message;
    message.precision(10);
    message << "Calculated a total normalized Rosseland integral of " << setw(12)
            << setiosflags(ios::fixed) << int_total;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message.precision(10);
    message << "Calculated a total normalized Rosseland integral of " << setw(12)
            << setiosflags(ios::fixed) << int_total << " instead of 1.0.";
    FAILMSG(message.str());
  }

  double int_1 = rtt_cdi::integratePlanckSpectrum(0.1, 1.0, 1.0);
  if (soft_equiv(int_1, .0345683, 1.0e-5)) {
    ostringstream message;
    message.precision(10);
    message << "Calculated a normalized Planck integral of " << setw(12) << setiosflags(ios::fixed)
            << int_1;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message.precision(10);
    message << "Calculated a normalized Planck integral of " << setw(12) << setiosflags(ios::fixed)
            << int_1 << " instead of .0345683";
    FAILMSG(message.str());
  }

  double int_2 = rtt_cdi::integrateRosselandSpectrum(0.1, 1.0, 1.0);
  if (soft_equiv(int_2, 0.01220025, 1.0e-5)) {
    ostringstream message;
    message.precision(10);
    message << "Calculated a normalized Rosseland integral of " << setw(12)
            << setiosflags(ios::fixed) << int_2;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message.precision(10);
    message << "Calculated a normalized Rosseland integral of " << setw(12)
            << setiosflags(ios::fixed) << int_2 << " instead of 0.01220025";
    FAILMSG(message.str());
  }

  double PL, ROS;
  rtt_cdi::integrate_Rosseland_Planckian_Spectrum(0.1, 1.0, 1.0, PL, ROS);
  if (soft_equiv(PL, .0345683, 1.e-5)) {
    ostringstream message;
    message.precision(10);
    message << "Calculated a normalized Planck integral for RosselandSpectrum " << setw(12)
            << setiosflags(ios::fixed) << int_2;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message.precision(10);
    message << "Calculated a normalized Planck integral for RosselandSpectrum " << setw(12)
            << setiosflags(ios::fixed) << int_2 << " instead of .0345683.";
    FAILMSG(message.str());
  }

  if (soft_equiv(ROS, 0.01220025, 1.e-5)) {
    ostringstream message;
    message.precision(10);
    message << "Calculated a normalized Rosseland integral for RosselandSpectrum " << setw(12)
            << setiosflags(ios::fixed) << int_2;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message.precision(10);
    message << "Calculated a normalized Rosseland integral for RosselandSpectrum " << setw(12)
            << setiosflags(ios::fixed) << int_2 << " instead of 0.01220025.";
    FAILMSG(message.str());
  }

  rtt_cdi::integrate_Rosseland_Planckian_Spectrum(0.1, 1.0, 1.0, PL, ROS);
  if (soft_equiv(PL, .0345683, 1.e-5)) {
    ostringstream message;
    message.precision(10);
    message << "Calculated a normalized Planck integral for "
               "RosselandPlanckianSpectrum "
            << setw(12) << setiosflags(ios::fixed) << int_2;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message.precision(10);
    message << "Calculated a normalized Planck integral for "
               "RosselandPlanckianSpectrum "
            << setw(12) << setiosflags(ios::fixed) << int_2 << " instead of .0345683.";
    FAILMSG(message.str());
  }

  if (soft_equiv(ROS, 0.01220025, 1.e-5)) {
    ostringstream message;
    message.precision(10);
    message << "Calculated a  normalized Rosseland integral for "
               "RosselandPlanckianSpectrum"
            << setw(12) << setiosflags(ios::fixed) << int_2;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message.precision(10);
    message << "Calculated a normalized Rosseland integral for "
               "RosselandPlanckianSpectrum"
            << setw(12) << setiosflags(ios::fixed) << int_2 << " instead of 0.01220025.";
    FAILMSG(message.str());
  }

  // Test that a zero temperature returns a zero
  if (soft_equiv(rtt_cdi::integrateRosselandSpectrum(0.0, 100.0, 0.0), 0.0))
    PASSMSG("Rosseland integral from hnu=0 to 100 at T=0 is zero: good!");
  else
    FAILMSG("Rosseland integral from hnu=0 to 100 at T=0 is not zero: BAD!");
  rtt_cdi::integrate_Rosseland_Planckian_Spectrum(0.0, 100.0, 0.0, PL, ROS);
  if (soft_equiv(PL, 0.0))
    PASSMSG("Rosseland call for Planck integral  at T=0 is zero: good!");
  else
    FAILMSG("Rosseland call for Planck integral at T=0 is not zero: BAD!");
  if (soft_equiv(ROS, 0.0))
    PASSMSG("Rosseland integral  at T=0 is zero: good!");
  else
    FAILMSG("Rosseland integral  at T=0 is not zero: BAD!");

  // check the normalized planck integrals
  if (CDI::getNumberFrequencyGroups() != 3)
    ITFAILS;

  // ----- Rosseland+Planckian ----- //

  // First group
  CDI::integrate_Rosseland_Planckian_Spectrum(1, 1.0, PL, ROS);
  if (!soft_equiv(PL, 0.005286862763740451, 1.e-6))
    ITFAILS;
  if (!soft_equiv(ROS, 0.00158258277444842, 1.e-5))
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("Group 1 Rosseland and Planck integrals ok.");
  else
    FAILMSG("Group 1 Rosseland and Planck integrals failed.");

  // Second group
  CDI::integrate_Rosseland_Planckian_Spectrum(2, 1.0, PL, ROS);
  if (!soft_equiv(PL, 0.7492399297, 1.e-6))
    ITFAILS;
  if (!soft_equiv(ROS, 0.5897280880, 1.e-6))
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("Group 2 Rosseland and Planck integrals ok.");
  else
    FAILMSG("Group 2 Rosseland and Planck integrals failed.");

  // Third group
  CDI::integrate_Rosseland_Planckian_Spectrum(3, 1.0, PL, ROS);
  if (!soft_equiv(PL, 0.2454669108, 1.e-6))
    ITFAILS;
  if (!soft_equiv(ROS, 0.4086877254, 1.e-6))
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("Group 3 Rosseland and Planck integrals ok.");
  else
    FAILMSG("Group 3 Rosseland and Planck integrals failed.");

  // All groups:
  std::vector<double> group_bounds(CDI::getFrequencyGroupBoundaries());
  if (group_bounds.size() != 4)
    ITFAILS;

  // Vectors for testing the combined planck/rosseland integration
  std::vector<double> planck;
  std::vector<double> rosseland;
  // Vector for testing the rosseland-only integration
  std::vector<double> rosseland_only;

  // Test using both combined Rosseland/Planckian integration and stand-alone Rosseland integration
  CDI::integrate_Rosseland_Planckian_Spectrum(group_bounds, 1.0, planck, rosseland);
  CDI::integrate_Rosseland_Spectrum(group_bounds, 1.0, rosseland_only);

  for (int group_index = 1; group_index <= 3; ++group_index) {
    CDI::integrate_Rosseland_Planckian_Spectrum(group_index, 1.0, PL, ROS);

    if (!soft_equiv(planck[group_index - 1], PL))
      ITFAILS;
    if (!soft_equiv(rosseland[group_index - 1], ROS))
      ITFAILS;
    if (!soft_equiv(rosseland_only[group_index - 1], ROS))
      ITFAILS;
  }

  // Special case of zero temperature
  CDI::integrate_Rosseland_Planckian_Spectrum(3, 0.0, PL, ROS);
  if (!soft_equiv(PL, 0.0, 1.e-6))
    ITFAILS;
  if (!soft_equiv(ROS, 0.0, 1.e-6))
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("Zero T Rosseland and Planck integrals ok.");
  else
    FAILMSG("Zero T Rosseland and Planck integrals failed.");

  if (ut.numFails == 0)
    PASSMSG(string("Group-wize and Full spectrum Planckian and ") + "Rosseland integrals match.");
  else
    FAILMSG(string("Group-wize and Full spectrum Planckian and ") +
            "Rosseland integrals do not match.");

  // More checks for extreme inputs

  // Extreme case 1. This should do the normal computation, but the result is zero.
  {
    std::vector<double> const bounds = {0.1, 0.3, 1.0, 3.0, 30.0};
    double const T_eval = 1.0e-30;
    std::vector<double> lrosseland;
    CDI::integrate_Rosseland_Spectrum(bounds, T_eval, lrosseland);
    if (!soft_equiv(lrosseland, std::vector<double>(lrosseland.size(), 0.0)))
      ITFAILS;
  }
  // Extreme case 2. T < numeric_limits<double>::min() --> follow special logic and return zero.
  bool caught = true;
  Remember(caught = false;);
  try {
    std::vector<double> const bounds = {0.1, 0.3, 1.0, 3.0, 30.0};
    double const T_eval = 1.0e-308;
    std::vector<double> lrosseland;
    CDI::integrate_Rosseland_Spectrum(bounds, T_eval, lrosseland);
  } catch (const rtt_dsxx::assertion &error) {
    ostringstream message;
    message << "Good, we caught the following exception: \n" << error.what();
    PASSMSG(message.str());
    caught = true;
  }
  if (!caught) {
    ostringstream message;
    message << "Failed to catch an exception when passing a denorm temperature.";
  }

  // Extreme case 3. This should do the normal computation, but the result is
  // zero.
  {
    std::vector<double> const bounds = {0.1, 0.3, 1.0, 3.0, 30.0};
    double const T_eval = 1.0e-30;
    std::vector<double> lplanck;
    std::vector<double> lrosseland;
    CDI::integrate_Rosseland_Planckian_Spectrum(bounds, T_eval, lplanck, lrosseland);
    if (!soft_equiv(lrosseland, std::vector<double>(lrosseland.size(), 0.0)))
      ITFAILS;
    if (!soft_equiv(lplanck, std::vector<double>(lplanck.size(), 0.0)))
      ITFAILS;
  }
  // Extreme case 4. T < numeric_limits<double>::min() --> follow special logic and return zero.
  caught = true;
  Remember(caught = false;);
  try {
    std::vector<double> const bounds = {0.1, 0.3, 1.0, 3.0, 30.0};
    double const T_eval = 1.0e-308;
    std::vector<double> lplanck;
    std::vector<double> lrosseland;
    CDI::integrate_Rosseland_Planckian_Spectrum(bounds, T_eval, lplanck, lrosseland);
  } catch (const rtt_dsxx::assertion &error) {
    ostringstream message;
    message << "Good, we caught the following exception: \n" << error.what();
    PASSMSG(message.str());
    caught = true;
  }
  if (!caught) {
    ostringstream message;
    message << "Failed to catch an exception when passing a denorm temperature.";
  }

  if (ut.numFails == 0) {
    PASSMSG("All Rosseland and Rosseland/Planckian integral tests ok.");
    cout << endl;
  }
  return;
}

//------------------------------------------------------------------------------------------------//
void test_mgopacity_collapse(rtt_dsxx::UnitTest &ut) {
  std::cout << "Running test test_mgopacity_collapse(ut)..." << std::endl;
  size_t const numFailCheckpoint(ut.numFails);

  // Test functions that collapse MG opacity data into one-group data using either Planckian or
  // Rosseland weight functions:

  std::shared_ptr<const MultigroupOpacity> mg_planck_abs(
      new DummyMultigroupOpacity(rtt_cdi::ABSORPTION, rtt_cdi::PLANCK));

  // bounds = { 0.05, 0.5, 5, 50 }
  std::vector<double> const bounds(mg_planck_abs->getGroupBoundaries());
  double const matTemp(1.0);
  double const matDens(1.0);
  // mgOpacities = { 3.64, 0.364, 0.0364 }
  std::vector<double> mgOpacities(mg_planck_abs->getOpacity(matTemp, matDens));

  // Vectors to hold the plankian and rosseland integrals for each cell.
  std::vector<double> planck_spectrum(bounds.size() - 1);
  std::vector<double> rosseland_spectrum(bounds.size() - 1);
  std::vector<double> rosseland_only_spectrum(bounds.size() - 1);
  std::vector<double> emission_group_cdf(bounds.size() - 1);

  // Simple test:
  {
    // Force the spectrum to be flat.
    planck_spectrum[0] = 1.0 / 3.0;
    planck_spectrum[1] = 1.0 / 3.0;
    planck_spectrum[2] = 1.0 / 3.0;
    rosseland_spectrum[0] = 1.0 / 3.0;
    rosseland_spectrum[1] = 1.0 / 3.0;
    rosseland_spectrum[2] = 1.0 / 3.0;

    // Generate reference solutions
    std::vector<double> emission_group_cdf_ref(bounds.size() - 1);
    double opacity_pl_ref(0.0);
    double opacity_pl_recip_ref(0.0);
    double opacity_ross_ref(0.0);
    for (size_t ig = 0; ig < bounds.size() - 1; ++ig) {
      opacity_pl_ref += planck_spectrum[ig] * mgOpacities[ig];
      opacity_pl_recip_ref += planck_spectrum[ig] / mgOpacities[ig];
      emission_group_cdf_ref[ig] = opacity_pl_ref;
      opacity_ross_ref += rosseland_spectrum[ig] / mgOpacities[ig];
    }
    opacity_ross_ref = 1.0 / opacity_ross_ref;

    // Collapse the opacities:
    double const opacity_pl = CDI::collapseMultigroupOpacitiesPlanck(
        bounds, mgOpacities, planck_spectrum, emission_group_cdf);
    // Make sure you get the same answer with the version that doesn't calculate emission CDF
    double const opacity_pl_alt =
        CDI::collapseMultigroupOpacitiesPlanck(bounds, mgOpacities, planck_spectrum);
    double const opacity_pl_recip =
        CDI::collapseMultigroupReciprocalOpacitiesPlanck(bounds, mgOpacities, planck_spectrum);
    double const opacity_ross =
        CDI::collapseMultigroupOpacitiesRosseland(bounds, mgOpacities, rosseland_spectrum);

    if (!soft_equiv(opacity_pl, opacity_pl_alt, 1.0e-12))
      ITFAILS;
    if (!soft_equiv(opacity_pl, opacity_pl_ref))
      ITFAILS;
    if (!soft_equiv(opacity_pl_recip, opacity_pl_recip_ref))
      ITFAILS;
    if (!soft_equiv(opacity_ross, opacity_ross_ref))
      ITFAILS;
    if (!soft_equiv(emission_group_cdf.begin(), emission_group_cdf.end(),
                    emission_group_cdf_ref.begin(), emission_group_cdf_ref.end()))
      ITFAILS;
  }

  // Standard use:
  {
    // Compute the planck and rosseland integrals for all groups in the spectrum with the
    // temperature of this cell.  This will return the following data: planck_spectrum
    //   = { 0.00528686276374045, 0.749239929709154, 0.24546691079067 }
    // rosseland_spectrum
    //   = { 0.00158258277444842, 0.589728087989656, 0.408687725376658 }
    CDI::integrate_Rosseland_Planckian_Spectrum(bounds, matTemp,                      // <- input
                                                planck_spectrum, rosseland_spectrum); // <- output
    // We also generate the stand-alone Rosseland spectrum for comparison (the rosseland and
    // rosseland_only vectors should be the same).
    CDI::integrate_Rosseland_Spectrum(bounds, matTemp,          // <- input
                                      rosseland_only_spectrum); // <- output

    // Collapse the opacities:
    double const opacity_pl = CDI::collapseMultigroupOpacitiesPlanck(
        bounds, mgOpacities, planck_spectrum, emission_group_cdf);
    // Make sure you get the same answer with the version that doesn't calculate emission CDF
    double const opacity_pl_alt =
        CDI::collapseMultigroupOpacitiesPlanck(bounds, mgOpacities, planck_spectrum);

    double const opacity_pl_recip =
        CDI::collapseMultigroupReciprocalOpacitiesPlanck(bounds, mgOpacities, planck_spectrum);
    double const opacity_ross =
        CDI::collapseMultigroupOpacitiesRosseland(bounds, mgOpacities, rosseland_spectrum);
    double const opacity_ross_only =
        CDI::collapseMultigroupOpacitiesRosseland(bounds, mgOpacities, rosseland_only_spectrum);

    std::vector<double> emission_group_cdf_ref(bounds.size() - 1);
    emission_group_cdf_ref[0] = 0.019244301636310527;
    emission_group_cdf_ref[1] = 0.29196935332821244;
    emission_group_cdf_ref[2] = 0.30090440514265909;
    double const opacity_pl_ref(0.30090440514265909);
    double const opacity_pl_recip_ref(8.80345577340399);
    double const opacity_ross_ref(0.0778314764921229);

    if (!soft_equiv(opacity_pl, opacity_pl_alt, 1.0e-12))
      ITFAILS;
    if (!soft_equiv(opacity_pl, opacity_pl_ref))
      ITFAILS;
    if (!soft_equiv(opacity_pl_recip, opacity_pl_recip_ref))
      ITFAILS;
    if (!soft_equiv(opacity_ross, opacity_ross_ref))
      ITFAILS;
    if (!soft_equiv(opacity_ross_only, opacity_ross_ref))
      ITFAILS;
    if (!soft_equiv(emission_group_cdf.begin(), emission_group_cdf.end(),
                    emission_group_cdf_ref.begin(), emission_group_cdf_ref.end()))
      ITFAILS;
  }

  // Special case 1 (opacity==0)
  {
    std::vector<double> mgOpac(3, 0.0);

    // Force the spectrum to be flat.
    planck_spectrum[0] = 1.0 / 3.0;
    planck_spectrum[1] = 1.0 / 3.0;
    planck_spectrum[2] = 1.0 / 3.0;
    rosseland_spectrum[0] = 1.0 / 3.0;
    rosseland_spectrum[1] = 1.0 / 3.0;
    rosseland_spectrum[2] = 1.0 / 3.0;

    // Generate reference solutions
    std::vector<double> emission_group_cdf_ref(bounds.size() - 1);
    double opacity_pl_ref(0.0);
    auto opacity_pl_recip_ref(static_cast<double>(std::numeric_limits<float>::max()));
    double opacity_ross_ref(0.0);
    for (size_t ig = 0; ig < bounds.size() - 1; ++ig) {
      emission_group_cdf_ref[ig] = opacity_pl_ref;
    }

    // Collapse the opacities:
    double const opacity_pl =
        CDI::collapseMultigroupOpacitiesPlanck(bounds, mgOpac, planck_spectrum, emission_group_cdf);
    // Make sure you get the same answer with the version that doesn't calculate emission CDF
    double const opacity_pl_alt =
        CDI::collapseMultigroupOpacitiesPlanck(bounds, mgOpac, planck_spectrum);

    double const opacity_pl_recip =
        CDI::collapseMultigroupReciprocalOpacitiesPlanck(bounds, mgOpac, planck_spectrum);
    double const opacity_ross =
        CDI::collapseMultigroupOpacitiesRosseland(bounds, mgOpac, rosseland_spectrum);

    if (!soft_equiv(opacity_pl, opacity_pl_alt, 1.0e-12))
      ITFAILS;
    if (!soft_equiv(opacity_pl, opacity_pl_ref))
      ITFAILS;
    if (!soft_equiv(opacity_pl_recip, opacity_pl_recip_ref))
      ITFAILS;
    if (!soft_equiv(opacity_ross, opacity_ross_ref))
      ITFAILS;
    if (!soft_equiv(emission_group_cdf.begin(), emission_group_cdf.end(),
                    emission_group_cdf_ref.begin(), emission_group_cdf_ref.end()))
      ITFAILS;
  }

  // Special case 2 (rosseland_spectrum == 0)
  {
    // Force the spectrum to be flat.
    planck_spectrum[0] = 1.0 / 3.0;
    planck_spectrum[1] = 1.0 / 3.0;
    planck_spectrum[2] = 1.0 / 3.0;
    rosseland_spectrum[0] = 0.0;
    rosseland_spectrum[1] = 0.0;
    rosseland_spectrum[2] = 0.0;

    // Generate reference solutions
    std::vector<double> emission_group_cdf_ref(bounds.size() - 1);
    double opacity_pl_ref(0.0);
    double opacity_ross_ref(mgOpacities[0]);
    for (size_t ig = 0; ig < bounds.size() - 1; ++ig) {
      opacity_pl_ref += planck_spectrum[ig] * mgOpacities[ig];
      emission_group_cdf_ref[ig] = opacity_pl_ref;
    }

    // Collapse the opacities:
    double const opacity_pl = CDI::collapseMultigroupOpacitiesPlanck(
        bounds, mgOpacities, planck_spectrum, emission_group_cdf);
    // Make sure you get the same answer with the version that doesn't calculate emission CDF
    double const opacity_pl_alt =
        CDI::collapseMultigroupOpacitiesPlanck(bounds, mgOpacities, planck_spectrum);

    double const opacity_ross =
        CDI::collapseMultigroupOpacitiesRosseland(bounds, mgOpacities, rosseland_spectrum);

    if (!soft_equiv(opacity_pl, opacity_pl_alt, 1.0e-12))
      ITFAILS;
    if (!soft_equiv(opacity_pl, opacity_pl_ref))
      ITFAILS;
    if (!soft_equiv(opacity_ross, opacity_ross_ref))
      ITFAILS;
    if (!soft_equiv(emission_group_cdf.begin(), emission_group_cdf.end(),
                    emission_group_cdf_ref.begin(), emission_group_cdf_ref.end()))
      ITFAILS;
  }

  // Report
  if (ut.numFails == numFailCheckpoint)
    PASSMSG("test_mgopacity_collapse completed successfully.");
  else
    FAILMSG("test_mgopacity_collapse reported at least one failure.");

  return;
}

//------------------------------------------------------------------------------------------------//
void test_cold(rtt_dsxx::UnitTest &ut) {
  // Test common corner case of cold cell

  std::vector<double> planck, rosseland;

  unsigned const ng = 3;
  std::vector<double> energyBoundary(ng + 1);
  energyBoundary[0] = 0.05;
  energyBoundary[1] = 0.5;
  energyBoundary[2] = 5.0;
  energyBoundary[3] = 50.0;

  cout << "cold case tests:" << endl;

  CDI::integrate_Planckian_Spectrum(energyBoundary, 0.0, planck);
  double norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 0.0), "planck extended integral, first form");

  planck = CDI::integrate_Planckian_Spectrum(energyBoundary, 0.0);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 0.0), "planck extended integral, second form");

  CDI::integrate_Rosseland_Spectrum(energyBoundary, 0.0, rosseland);
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 0.0), "rosseland extended integral");

  CDI::integrate_Rosseland_Planckian_Spectrum(energyBoundary, 0.0, planck, rosseland);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 0.0), "planck extended integral, paired");
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 0.0), "rosseland extended integral, paired");

  energyBoundary[0] = 0.0;

  CDI::integrate_Planckian_Spectrum(energyBoundary, 0.0, planck);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "planck extended integral, first form");

  planck = CDI::integrate_Planckian_Spectrum(energyBoundary, 0.0);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "planck extended integral, second form");

  CDI::integrate_Rosseland_Spectrum(energyBoundary, 0.0, rosseland);
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "rosseland extended integral");

  CDI::integrate_Rosseland_Planckian_Spectrum(energyBoundary, 0.0, planck, rosseland);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "planck extended integral, paired");
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "rosseland extended integral, paired");
}

//------------------------------------------------------------------------------------------------//
void test_hot(rtt_dsxx::UnitTest &ut) {
  // Test common corner case of hot cell

  std::vector<double> planck, rosseland;

  unsigned const ng = 3;
  std::vector<double> energyBoundary(ng + 1);
  energyBoundary[0] = 0.05;
  energyBoundary[1] = 0.5;
  energyBoundary[2] = 5.0;
  energyBoundary[3] = 50.0;

  CDI::integrate_Planckian_Spectrum(energyBoundary, 1e100, planck);
  double norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 0.0), "hot planck integral, first form");

  planck = CDI::integrate_Planckian_Spectrum(energyBoundary, 1e100);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 0.0), "hot planck integral, second form");

  CDI::integrate_Rosseland_Spectrum(energyBoundary, 1e100, rosseland);
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 0.0), "hot rosseland integral");

  CDI::integrate_Rosseland_Planckian_Spectrum(energyBoundary, 1e100, planck, rosseland);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 0.0), "hot planck integral, paired");
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 0.0), "hot rosseland integral, paired");

  energyBoundary[0] = 0.0;

  CDI::integrate_Planckian_Spectrum(energyBoundary, 1e100, planck);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 0.0), "hot planck integral, first form");

  planck = CDI::integrate_Planckian_Spectrum(energyBoundary, 1e100);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 0.0), "hot planck integral, second form");

  CDI::integrate_Rosseland_Spectrum(energyBoundary, 1e100, rosseland);
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 0.0), "hot rosseland integral");

  CDI::integrate_Rosseland_Planckian_Spectrum(energyBoundary, 1e100, planck, rosseland);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 0.0), "hot planck integral, paired");
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 0.0), "hot rosseland integral, paired");
}

//------------------------------------------------------------------------------------------------//
void test_extend(rtt_dsxx::UnitTest &ut) {
  // Test extended integration option

  CDI::setExtend();
  ut.check(CDI::getExtend(), "extend boundaries flag");

  std::vector<double> planck, rosseland;

  unsigned const ng = 3;
  std::vector<double> energyBoundary(ng + 1);
  energyBoundary[0] = 0.05;
  energyBoundary[1] = 0.5;
  energyBoundary[2] = 5.0;
  energyBoundary[3] = 50.0;

  CDI::integrate_Planckian_Spectrum(energyBoundary, 10.0, planck);
  double norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "planck extended integral, first form");

  planck = CDI::integrate_Planckian_Spectrum(energyBoundary, 5.0);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "planck extended integral, second form");

  CDI::integrate_Rosseland_Spectrum(energyBoundary, 10.0, rosseland);
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "rosseland extended integral");

  CDI::integrate_Rosseland_Planckian_Spectrum(energyBoundary, 5.0, planck, rosseland);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "planck extended integral, paired");
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "rosseland extended integral, paired");
}

//------------------------------------------------------------------------------------------------//
void test_extend_cold(rtt_dsxx::UnitTest &ut) {
  // Test extended integration option

  std::vector<double> planck, rosseland;

  unsigned const ng = 3;
  std::vector<double> energyBoundary(ng + 1);
  energyBoundary[0] = 0.05;
  energyBoundary[1] = 0.5;
  energyBoundary[2] = 5.0;
  energyBoundary[3] = 50.0;

  CDI::integrate_Planckian_Spectrum(energyBoundary, 0.0, planck);
  double norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "cold planck extended integral, first form");

  planck = CDI::integrate_Planckian_Spectrum(energyBoundary, 5.0);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "cold planck extended integral, second form");

  CDI::integrate_Rosseland_Spectrum(energyBoundary, 0.0, rosseland);
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "cold rosseland extended integral");

  CDI::integrate_Rosseland_Planckian_Spectrum(energyBoundary, 5.0, planck, rosseland);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "cold planck extended integral, paired");
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "cold rosseland extended integral, paired");

  energyBoundary[0] = 0.0;
  CDI::integrate_Planckian_Spectrum(energyBoundary, 0.0, planck);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "cold planck extended integral, first form");

  planck = CDI::integrate_Planckian_Spectrum(energyBoundary, 5.0);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "cold planck extended integral, second form");

  CDI::integrate_Rosseland_Spectrum(energyBoundary, 0.0, rosseland);
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "cold rosseland extended integral");

  CDI::integrate_Rosseland_Planckian_Spectrum(energyBoundary, 5.0, planck, rosseland);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "cold planck extended integral, paired");
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "cold rosseland extended integral, paired");
}

//------------------------------------------------------------------------------------------------//
void test_extend_hot(rtt_dsxx::UnitTest &ut) {
  // Test extended integration option

  std::vector<double> planck, rosseland;

  unsigned const ng = 3;
  std::vector<double> energyBoundary(ng + 1);
  energyBoundary[0] = 0.05;
  energyBoundary[1] = 0.5;
  energyBoundary[2] = 5.0;
  energyBoundary[3] = 50.0;

  CDI::integrate_Planckian_Spectrum(energyBoundary, 1e100, planck);
  double norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "hot planck extended integral, first form");

  planck = CDI::integrate_Planckian_Spectrum(energyBoundary, 1e100);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "hot planck extended integral, second form");

  CDI::integrate_Rosseland_Spectrum(energyBoundary, 1e100, rosseland);
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "hot rosseland extended integral");

  CDI::integrate_Rosseland_Planckian_Spectrum(energyBoundary, 1e100, planck, rosseland);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "hot planck extended integral, paired");
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "hot rosseland extended integral, paired");

  energyBoundary[0] = 0.0;
  CDI::integrate_Planckian_Spectrum(energyBoundary, 1e100, planck);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "hot planck extended integral, first form");

  planck = CDI::integrate_Planckian_Spectrum(energyBoundary, 1e100);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "hot planck extended integral, second form");

  CDI::integrate_Rosseland_Spectrum(energyBoundary, 1e100, rosseland);
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "hot rosseland extended integral");

  CDI::integrate_Rosseland_Planckian_Spectrum(energyBoundary, 1e100, planck, rosseland);
  norm = std::accumulate(planck.begin(), planck.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "hot planck extended integral, paired");
  norm = std::accumulate(rosseland.begin(), rosseland.end(), 0.0);
  ut.check(soft_equiv(norm, 1.0), "hot rosseland extended integral, paired");
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_CDI(ut);
    test_planck_integration(ut);
    test_rosseland_integration(ut);
    test_mgopacity_collapse(ut);
    test_cold(ut);
    test_hot(ut);
    test_extend(ut);
    test_extend_cold(ut);
    test_extend_hot(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tCDI.cc
//------------------------------------------------------------------------------------------------//
