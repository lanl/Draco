//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/test/tstAnalytic_CP_Eloss.cc
 * \author Kendra P. Long
 * \date   Fri Aug  2 14:28:08 2019
 * \brief  Analytic_CP_Eloss test.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi/CDI.hh"
#include "cdi_CPEloss/Analytic_CP_Eloss.hh"
#include "cdi_CPEloss/Analytic_Constant_Eloss_Model.hh"
#include "cdi_CPEloss/Analytic_KP_Alpha_Eloss_Model.hh"
#include "cdi_CPEloss/Analytic_Spitzer_Eloss_Model.hh"
#include "cdi_CPEloss/Analytic_TR_Eloss_Model.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/dbc.hh"
#include <sstream>

using namespace std;

using rtt_cdi_cpeloss::Analytic_CP_Eloss;
using rtt_cdi_cpeloss::Analytic_Eloss_Model;
//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void KP_alpha_test(rtt_dsxx::UnitTest &ut) {

  // Deuterium:
  int32_t deuterium_zaid = 1002;
  double deuterium_mass = 3.34358e-24;
  rtt_cdi::CParticle target_in(deuterium_zaid, deuterium_mass);
  // Alpha particle:
  int32_t alpha_zaid = 2004;
  double alpha_mass = 6.64424e-24;
  rtt_cdi::CParticle projectile_in(alpha_zaid, alpha_mass);

  std::shared_ptr<Analytic_Eloss_Model> model_in(
      new rtt_cdi_cpeloss::Analytic_KP_Alpha_Eloss_Model(target_in, projectile_in));

  Analytic_CP_Eloss eloss_mod(model_in, target_in, projectile_in,
                              rtt_cdi::CPModelAngleCutoff::NONE);

  // Check that basic accessors return correct result:
  // Analytic model should match that passed to constructor
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_mod.getEloss(1., 10., 1.),
                                   eloss_mod.get_Analytic_Model()->calculate_eloss(1., 10., 1.),
                                   1.0e-3));

  // Model type better be analytic:
  FAIL_IF_NOT(eloss_mod.getModelType() == rtt_cdi::CPModelType::ANALYTIC_ETYPE);

  // NOT tabular data
  FAIL_IF(eloss_mod.is_data_in_tabular_form());

  // All sizes should be 0 (again, not tabular data)
  FAIL_IF_NOT(eloss_mod.getTemperatureGrid().size() == eloss_mod.getNumTemperatures());
  FAIL_IF_NOT(eloss_mod.getDensityGrid().size() == eloss_mod.getNumDensities());
  FAIL_IF_NOT(eloss_mod.getEnergyGrid().size() == eloss_mod.getNumEnergies());

  // Check that accessors return the correct target and projectile:
  FAIL_IF_NOT(target_in.get_zaid() == eloss_mod.getTarget().get_zaid());
  FAIL_IF_NOT(projectile_in.get_zaid() == eloss_mod.getProjectile().get_zaid());

  // Check that accessor returns the correct model angle cutoff
  FAIL_IF_NOT(eloss_mod.getModelAngleCutoff() == rtt_cdi::CPModelAngleCutoff::NONE);

  // Get eloss values for some sample data:
  {
    double T = 1.0;    // keV
    double rho = 10.0; // g / cc
    double vel0 = 1.0; // cm / shk

    double eloss_coeff = eloss_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 212.287, 1.0e-3));
  }

  // Point near the maximum on rho in [0, 20] g/cc, T in [0, 100] keV:
  {
    double T = 15.0;   // keV
    double rho = 20.0; // g / cc
    double vel0 = 1.0; // cm / shk

    double eloss_coeff = eloss_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 78.5845, 1.0e-4));
  }

  // Point where the range fit goes negative -- eloss should return 1.0e25:
  // (equivalent to the particle range limiting to a very small number)
  {
    double T = 90.0;   // keV
    double rho = 17.5; // g / cc
    double vel0 = 1.0; // cm / shk

    double eloss_coeff = eloss_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 1.0e25));
  }

  if (ut.numFails == 0)
    PASSMSG("KP_Alpha CPEloss test passes.");
  else
    FAILMSG("KP_Alpha CPEloss test fails.");

  return;
}

void Spitzer_test(rtt_dsxx::UnitTest &ut) {

  // Deuterium:
  int32_t deuterium_zaid = 1002;
  double deuterium_mass = 3.34358e-24;
  rtt_cdi::CParticle target_in(deuterium_zaid, deuterium_mass);
  // Alpha particle:
  int32_t alpha_zaid = 2004;
  double alpha_mass = 6.64424e-24;
  rtt_cdi::CParticle projectile_in(alpha_zaid, alpha_mass);

  std::shared_ptr<Analytic_Eloss_Model> model_in(
      new rtt_cdi_cpeloss::Analytic_Spitzer_Eloss_Model(target_in, projectile_in));

  Analytic_CP_Eloss eloss_mod(model_in, target_in, projectile_in,
                              rtt_cdi::CPModelAngleCutoff::NONE);

  // Check that basic accessors return correct result:
  // Analytic model should match that passed to constructor
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_mod.getEloss(1., 10., 1.),
                                   eloss_mod.get_Analytic_Model()->calculate_eloss(1., 10., 1.),
                                   1.0e-3));

  // Model type better be analytic:
  FAIL_IF_NOT(eloss_mod.getModelType() == rtt_cdi::CPModelType::ANALYTIC_ETYPE);

  // NOT tabular data
  FAIL_IF(eloss_mod.is_data_in_tabular_form());

  // All sizes should be 0 (again, not tabular data)
  FAIL_IF_NOT(eloss_mod.getTemperatureGrid().size() == eloss_mod.getNumTemperatures());
  FAIL_IF_NOT(eloss_mod.getDensityGrid().size() == eloss_mod.getNumDensities());
  FAIL_IF_NOT(eloss_mod.getEnergyGrid().size() == eloss_mod.getNumEnergies());

  // Check that accessors return the correct target and projectile:
  FAIL_IF_NOT(target_in.get_zaid() == eloss_mod.getTarget().get_zaid());
  FAIL_IF_NOT(projectile_in.get_zaid() == eloss_mod.getProjectile().get_zaid());

  // Check that accessor returns the correct model angle cutoff
  FAIL_IF_NOT(eloss_mod.getModelAngleCutoff() == rtt_cdi::CPModelAngleCutoff::NONE);

  // Get eloss values for some sample data:
  {
    double T = 1.0;    // keV
    double rho = 10.0; // g / cc
    double vel0 = 1.0; // cm / shk

    double eloss_coeff = eloss_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 9.823342835413329303e+05, 1.0e-8));
  }

  // Check for large value when 2kT > E
  {
    double T = 1.e2;   // keV
    double rho = 10.0; // g / cc
    double vel0 = 1.0; // cm / shk

    double eloss_coeff = eloss_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 1.e15, 1.0e-8));
  }

  if (ut.numFails == 0)
    PASSMSG("Spitzer CPEloss test passes.");
  else
    FAILMSG("Spitzer CPEloss test fails.");

  return;
}

void TR_test(rtt_dsxx::UnitTest &ut) {

  // Deuterium:
  int32_t deuterium_zaid = 1002;
  double deuterium_mass = 3.34358e-24;
  rtt_cdi::CParticle target_ion(deuterium_zaid, deuterium_mass);
  // Electron:
  int32_t electron_zaid = -1;
  double electron_mass = 9.10938291e-28;
  rtt_cdi::CParticle target_electron(electron_zaid, electron_mass);
  // Alpha particle:
  int32_t alpha_zaid = 2004;
  double alpha_mass = 6.64424e-24;
  rtt_cdi::CParticle projectile_in(alpha_zaid, alpha_mass);

  std::shared_ptr<Analytic_Eloss_Model> ion_model_in(
      new rtt_cdi_cpeloss::Analytic_TR_Eloss_Model(target_ion, projectile_in));
  std::shared_ptr<Analytic_Eloss_Model> electron_model_in(
      new rtt_cdi_cpeloss::Analytic_TR_Eloss_Model(target_electron, projectile_in));

  Analytic_CP_Eloss eloss_ion_mod(ion_model_in, target_ion, projectile_in,
                                  rtt_cdi::CPModelAngleCutoff::NONE);
  Analytic_CP_Eloss eloss_ele_mod(electron_model_in, target_electron, projectile_in,
                                  rtt_cdi::CPModelAngleCutoff::NONE);

  // Check that basic accessors return correct result:
  // Analytic model should match that passed to constructor
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_ion_mod.getEloss(1., 10., 1.),
                                   eloss_ion_mod.get_Analytic_Model()->calculate_eloss(1., 10., 1.),
                                   1.0e-3));

  // Model type better be analytic:
  FAIL_IF_NOT(eloss_ion_mod.getModelType() == rtt_cdi::CPModelType::ANALYTIC_ETYPE);

  // NOT tabular data
  FAIL_IF(eloss_ion_mod.is_data_in_tabular_form());

  // All sizes should be 0 (again, not tabular data)
  FAIL_IF_NOT(eloss_ion_mod.getTemperatureGrid().size() == eloss_ion_mod.getNumTemperatures());
  FAIL_IF_NOT(eloss_ion_mod.getDensityGrid().size() == eloss_ion_mod.getNumDensities());
  FAIL_IF_NOT(eloss_ion_mod.getEnergyGrid().size() == eloss_ion_mod.getNumEnergies());

  // Check that accessors return the correct target and projectile:
  FAIL_IF_NOT(target_ion.get_zaid() == eloss_ion_mod.getTarget().get_zaid());
  FAIL_IF_NOT(target_electron.get_zaid() == eloss_ele_mod.getTarget().get_zaid());
  FAIL_IF_NOT(projectile_in.get_zaid() == eloss_ion_mod.getProjectile().get_zaid());

  // Check that accessor returns the correct model angle cutoff
  FAIL_IF_NOT(eloss_ion_mod.getModelAngleCutoff() == rtt_cdi::CPModelAngleCutoff::NONE);

  // Get eloss values for some sample ion data:
  {
    double T = 1.0;    // keV
    double rho = 10.0; // g / cc
    double vel0 = 1.0; // cm / shk

    double eloss_coeff = eloss_ion_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 9.927570982586801983e+05, 1.0e-8));
  }

  // Get eloss values for some sample electron data:
  {
    double T = 1.0;     // keV
    double rho = 1.e-2; // g / cc
    double vel0 = 1.0;  // cm / shk

    double eloss_coeff = eloss_ele_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 4.993807286690149340e+05, 1.0e-8));
  }

  // Check for large value when 2kT > E
  {
    double T = 1.e2;   // keV
    double rho = 10.0; // g / cc
    double vel0 = 1.0; // cm / shk

    double eloss_coeff = eloss_ion_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 1.e15, 1.0e-8));
  }

  if (ut.numFails == 0)
    PASSMSG("TR CPEloss test passes.");
  else
    FAILMSG("TR CPEloss test fails.");

  return;
}

void Constant_test(rtt_dsxx::UnitTest &ut) {

  // Deuterium:
  int32_t deuterium_zaid = 1002;
  double deuterium_mass = 3.34358e-24;
  rtt_cdi::CParticle target_in(deuterium_zaid, deuterium_mass);
  // Alpha particle:
  int32_t alpha_zaid = 2004;
  double alpha_mass = 6.64424e-24;
  rtt_cdi::CParticle projectile_in(alpha_zaid, alpha_mass);

  std::shared_ptr<Analytic_Eloss_Model> zero_model_in(
      new rtt_cdi_cpeloss::Analytic_Constant_Eloss_Model(target_in, projectile_in, 0.));
  std::shared_ptr<Analytic_Eloss_Model> one_model_in(
      new rtt_cdi_cpeloss::Analytic_Constant_Eloss_Model(target_in, projectile_in, 1.));

  Analytic_CP_Eloss eloss_zero_mod(zero_model_in, target_in, projectile_in,
                                   rtt_cdi::CPModelAngleCutoff::NONE);
  Analytic_CP_Eloss eloss_one_mod(one_model_in, target_in, projectile_in,
                                  rtt_cdi::CPModelAngleCutoff::NONE);

  // Check that basic accessors return correct result:
  // Analytic model should match that passed to constructor
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_zero_mod.getEloss(1., 2., 3.), 0., 1.0e-3));

  FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_one_mod.getEloss(1., 2., 3.), 1., 1.0e-3));

  if (ut.numFails == 0)
    PASSMSG("Constant CPEloss test passes.");
  else
    FAILMSG("Constant CPEloss test fails.");

  return;
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    KP_alpha_test(ut);
    Spitzer_test(ut);
    TR_test(ut);
    Constant_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstAnalytic_CP_Eloss.cc
//------------------------------------------------------------------------------------------------//
