//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/test/tstAnalytic_CP_Eloss.cc
 * \author Kendra P. Long
 * \date   Fri Aug  2 14:28:08 2019
 * \brief  Analytic_CP_Eloss test.
 * \note   Copyright (C) 2020-2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "cdi/CPCommon.hh"
#include "cdi_CPEloss/Analytic_CP_Eloss.hh"
#include "cdi_CPEloss/Analytic_Constant_Eloss_Model.hh"
#include "cdi_CPEloss/Analytic_Eloss_Model.hh"
#include "cdi_CPEloss/Analytic_KP_Alpha_Eloss_Model.hh"
#include "cdi_CPEloss/Analytic_Spitzer_Eloss_Model.hh"
#include "cdi_CPEloss/Analytic_TR_Eloss_Model.hh"
#include <cstdint>
#include <memory>

using namespace std;
using rtt_cdi_cpeloss::Analytic_CP_Eloss;
using rtt_cdi_cpeloss::Analytic_Eloss_Model;

// Deuterium:
int32_t constexpr deuterium_zaid{1002};
double constexpr deuterium_mass{3.34358e-24};
// Alpha particle:
int32_t constexpr alpha_zaid{2004};
double constexpr alpha_mass{6.64424e-24};
// Electron:
int32_t constexpr electron_zaid{-1};
double constexpr electron_mass{9.10938291e-28};

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void KP_alpha_test(rtt_dsxx::UnitTest &ut) {

  rtt_cdi::CParticle const target_in(deuterium_zaid, deuterium_mass);
  rtt_cdi::CParticle const projectile_in(alpha_zaid, alpha_mass);

  std::shared_ptr<Analytic_Eloss_Model> const model_in(
      new rtt_cdi_cpeloss::Analytic_KP_Alpha_Eloss_Model(target_in, projectile_in));

  Analytic_CP_Eloss const eloss_mod(model_in, target_in, projectile_in,
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
    double const T = 1.0;                             // keV
    double const ndens = 10.0 / target_in.get_mass(); // a / cc
    double const vel0 = 1.0;                          // cm / shk

    double const eloss_coeff = eloss_mod.getEloss(T, ndens, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 212.287, 1.0e-3));
  }

  // Point near the maximum on rho in [0, 20] g/cc, T in [0, 100] keV:
  {
    double const T = 15.0;                            // keV
    double const ndens = 20.0 / target_in.get_mass(); // a / cc
    double const vel0 = 1.0;                          // cm / shk

    double const eloss_coeff = eloss_mod.getEloss(T, ndens, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 78.5845, 1.0e-4));
  }

  // Point where the range fit goes negative -- eloss should return 1.0e25:
  // (equivalent to the particle range limiting to a very small number)
  {
    double const T = 90.0;                            // keV
    double const ndens = 17.5 / target_in.get_mass(); // a / cc
    double const vel0 = 1.0;                          // cm / shk

    double const eloss_coeff = eloss_mod.getEloss(T, ndens, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 1.0e25));
  }

  if (ut.numFails == 0)
    PASSMSG("KP_Alpha CPEloss test passes.");
  else
    FAILMSG("KP_Alpha CPEloss test fails.");

  return;
}

//------------------------------------------------------------------------------------------------//
void Spitzer_test(rtt_dsxx::UnitTest &ut) {

  rtt_cdi::CParticle const target_in(deuterium_zaid, deuterium_mass);
  rtt_cdi::CParticle const projectile_in(alpha_zaid, alpha_mass);

  std::shared_ptr<Analytic_Eloss_Model> const model_in(
      new rtt_cdi_cpeloss::Analytic_Spitzer_Eloss_Model(target_in, projectile_in));

  Analytic_CP_Eloss const eloss_mod(model_in, target_in, projectile_in,
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
    double const T = 1.0;                             // keV
    double const ndens = 10.0 / target_in.get_mass(); // a / cc
    double const vel0 = 1.0;                          // cm / shk

    double const eloss_coeff = eloss_mod.getEloss(T, ndens, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 9.823342835413329303e+05, 1.0e-8));
  }

  // Check for large value when 2kT > E
  {
    double const T = 1.e2;     // keV
    double const ndens = 10.0; // g / cc
    double const vel0 = 1.0;   // cm / shk

    double const eloss_coeff = eloss_mod.getEloss(T, ndens, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 1.e15, 1.0e-8));
  }

  if (ut.numFails == 0)
    PASSMSG("Spitzer CPEloss test passes.");
  else
    FAILMSG("Spitzer CPEloss test fails.");

  return;
}

//------------------------------------------------------------------------------------------------//
void TR_test(rtt_dsxx::UnitTest &ut) {

  rtt_cdi::CParticle const target_ion(deuterium_zaid, deuterium_mass);
  rtt_cdi::CParticle const target_electron(electron_zaid, electron_mass);
  rtt_cdi::CParticle const projectile_in(alpha_zaid, alpha_mass);

  std::shared_ptr<Analytic_Eloss_Model> const ion_model_in(
      new rtt_cdi_cpeloss::Analytic_TR_Eloss_Model(target_ion, projectile_in));
  std::shared_ptr<Analytic_Eloss_Model> const electron_model_in(
      new rtt_cdi_cpeloss::Analytic_TR_Eloss_Model(target_electron, projectile_in));

  Analytic_CP_Eloss const eloss_ion_mod(ion_model_in, target_ion, projectile_in,
                                        rtt_cdi::CPModelAngleCutoff::NONE);
  Analytic_CP_Eloss const eloss_ele_mod(electron_model_in, target_electron, projectile_in,
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
    double const T = 1.0;                              // keV
    double const ndens = 10.0 / target_ion.get_mass(); // a / cc
    double const vel0 = 1.0;                           // cm / shk

    double const eloss_coeff = eloss_ion_mod.getEloss(T, ndens, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 9.927570982586801983e+05, 1.0e-8));
  }

  // Get eloss values for some sample electron data:
  {
    double const T = 1.0;                                    // keV
    double const ndens = 1.e-2 / target_electron.get_mass(); // el / cc
    double const vel0 = 1.0;                                 // cm / shk

    double const eloss_coeff = eloss_ele_mod.getEloss(T, ndens, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 4.993807286690149340e+05, 1.0e-8));
  }

  // Check for large value when 2kT > E
  {
    double const T = 1.e2;                             // keV
    double const ndens = 10.0 / target_ion.get_mass(); // g / cc
    double const vel0 = 1.0;                           // cm / shk

    double const eloss_coeff = eloss_ion_mod.getEloss(T, ndens, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 1.e15, 1.0e-8));
  }

  if (ut.numFails == 0)
    PASSMSG("TR CPEloss test passes.");
  else
    FAILMSG("TR CPEloss test fails.");

  return;
}

//------------------------------------------------------------------------------------------------//
void Constant_test(rtt_dsxx::UnitTest &ut) {

  rtt_cdi::CParticle const target_in(deuterium_zaid, deuterium_mass);
  rtt_cdi::CParticle const projectile_in(alpha_zaid, alpha_mass);

  std::shared_ptr<Analytic_Eloss_Model> const zero_model_in(
      new rtt_cdi_cpeloss::Analytic_Constant_Eloss_Model(target_in, projectile_in, 0.));
  std::shared_ptr<Analytic_Eloss_Model> const one_model_in(
      new rtt_cdi_cpeloss::Analytic_Constant_Eloss_Model(target_in, projectile_in, 1.));

  Analytic_CP_Eloss const eloss_zero_mod(zero_model_in, target_in, projectile_in,
                                         rtt_cdi::CPModelAngleCutoff::NONE);
  Analytic_CP_Eloss const eloss_one_mod(one_model_in, target_in, projectile_in,
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
