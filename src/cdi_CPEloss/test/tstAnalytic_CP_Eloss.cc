//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/test/tstAnalytic_CP_Eloss.cc
 * \author Kendra P. Long
 * \brief  Analytic_CP_Eloss test.
 * \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "cdi/CPCommon.hh"
#include "cdi_CPEloss/Analytic_CP_Eloss.hh"
#include "cdi_CPEloss/Analytic_Constant_Eloss_Model.hh"
#include "cdi_CPEloss/Analytic_KP_Alpha_Eloss_Model.hh"
#include "cdi_CPEloss/Analytic_LiPet_Eloss_Model.hh"
#include "cdi_CPEloss/Analytic_Spitzer_Eloss_Model.hh"
#include "cdi_CPEloss/Analytic_TR_Eloss_Model.hh"
#include "cdi_CPEloss/config.h" // Defines GSL_FOUND
#include "units/PhysicalConstantsSI.hh"
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

namespace {

using namespace std;
using rtt_cdi_cpeloss::Analytic_CP_Eloss;

// Deuterium:
int32_t constexpr deuterium_zaid{1002};
double constexpr deuterium_mass{3.34358e-24};
// Tritium:
int32_t constexpr tritium_zaid{1003};
double constexpr tritium_mass{5.00736e-24};
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

  auto const model_in(std::make_shared<rtt_cdi_cpeloss::Analytic_KP_Alpha_Eloss_Model>(
      deuterium_zaid, deuterium_mass, alpha_zaid, alpha_mass));

  Analytic_CP_Eloss const eloss_mod(model_in, rtt_cdi::CPModelAngleCutoff::NONE);

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
  FAIL_IF_NOT(deuterium_zaid == eloss_mod.targ_zaid());
  FAIL_IF_NOT(alpha_zaid == eloss_mod.proj_zaid());

  // Check that accessor returns the correct model angle cutoff
  FAIL_IF_NOT(eloss_mod.getModelAngleCutoff() == rtt_cdi::CPModelAngleCutoff::NONE);

  // Get eloss values for some sample data:
  {
    double const T = 1.0;                       // keV
    double const ndens = 10.0 / deuterium_mass; // a / cc
    double const vel0 = 1.0;                    // cm / shk

    double const eloss_coeff = eloss_mod.getEloss(T, ndens, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 212.287, 1.0e-3));
  }

  // Point near the maximum on rho in [0, 20] g/cc, T in [0, 100] keV:
  {
    double const T = 15.0;                      // keV
    double const ndens = 20.0 / deuterium_mass; // a / cc
    double const vel0 = 1.0;                    // cm / shk

    double const eloss_coeff = eloss_mod.getEloss(T, ndens, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 78.5845, 1.0e-4));
  }

  // Point where the range fit goes negative -- eloss should return 1.0e25:
  // (equivalent to the particle range limiting to a very small number)
  {
    double const T = 90.0;                      // keV
    double const ndens = 17.5 / deuterium_mass; // a / cc
    double const vel0 = 1.0;                    // cm / shk

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
  // Temperature in keV
  double const temp_electron = 1.0;
  double const temp_deuteron = 1.0;
  double const temp_triton = 1.0;

  // Plasma Density g/cm^3
  double const dens_electron = 1.0;
  double const dens_deuteron = 0.5;
  double const dens_triton = 0.5;

  // Number densities
  double const molar_mass = (0.5 * (2.013553214 + 3.01550071913));
  double const number_dens_electron = dens_electron / molar_mass * rtt_units::AVOGADRO;
  double const number_dens_deuteron = dens_deuteron / molar_mass * rtt_units::AVOGADRO;
  double const number_dens_triton = dens_triton / molar_mass * rtt_units::AVOGADRO;

  // Build cpeloss model
  auto const elec_model(std::make_shared<rtt_cdi_cpeloss::Analytic_Spitzer_Eloss_Model>(
      electron_zaid, electron_mass, alpha_zaid, alpha_mass));
  auto const deut_model(std::make_shared<rtt_cdi_cpeloss::Analytic_Spitzer_Eloss_Model>(
      deuterium_zaid, deuterium_mass, alpha_zaid, alpha_mass));
  auto const trit_model(std::make_shared<rtt_cdi_cpeloss::Analytic_Spitzer_Eloss_Model>(
      tritium_zaid, tritium_mass, alpha_zaid, alpha_mass));

  // Build CPEloss object using model:
  rtt_cdi_cpeloss::Analytic_CP_Eloss const eloss_elec(elec_model,
                                                      rtt_cdi::CPModelAngleCutoff::NONE);
  rtt_cdi_cpeloss::Analytic_CP_Eloss const eloss_deut(deut_model,
                                                      rtt_cdi::CPModelAngleCutoff::NONE);
  rtt_cdi_cpeloss::Analytic_CP_Eloss const eloss_trit(trit_model,
                                                      rtt_cdi::CPModelAngleCutoff::NONE);

  // Check incident particle info
  FAIL_IF_NOT(eloss_elec.proj_zaid() == alpha_zaid);
  FAIL_IF_NOT(eloss_deut.proj_zaid() == alpha_zaid);
  FAIL_IF_NOT(eloss_trit.proj_zaid() == alpha_zaid);
  FAIL_IF_NOT(eloss_elec.targ_zaid() == electron_zaid);
  FAIL_IF_NOT(eloss_deut.targ_zaid() == deuterium_zaid);
  FAIL_IF_NOT(eloss_trit.targ_zaid() == tritium_zaid);

  FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_elec.proj_mass(), alpha_mass, 1.0e-8));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_deut.proj_mass(), alpha_mass, 1.0e-8));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_trit.proj_mass(), alpha_mass, 1.0e-8));

  FAIL_IF_NOT(rtt_cdi::get_z_number(eloss_elec.get_Analytic_Model()->proj_zaid()) == 2);

  // Get dedt for 3.5 MeV alpha-particle
  {
    double const incident_energy = 3.5 * rtt_units::electronChargeSI * 1.0e+6; // J
    double const incident_speed =
        std::sqrt(2.0 * incident_energy / (alpha_mass * 1.0e-3)) * 1.0e-8 * 1.0e+2; // cm/shk

    double dedt{0.0};
    dedt += eloss_elec.getEloss(temp_electron, number_dens_electron, incident_speed);
    dedt += eloss_deut.getEloss(temp_deuteron, number_dens_deuteron, incident_speed);
    dedt += eloss_trit.getEloss(temp_triton, number_dens_triton, incident_speed);

    const double ref_dedt{3.05650072e+05 * incident_speed};
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(dedt, ref_dedt, 1.0e-8));
  }

  // Get dedt for 1.0 MeV alpha-particle
  {
    double const incident_energy = 1.0 * rtt_units::electronChargeSI * 1.0e+6; // J
    double const incident_speed =
        std::sqrt(2.0 * incident_energy / (alpha_mass * 1.0e-3)) * 1.0e-8 * 1.0e+2; // cm/shk

    double dedt{0.0};
    dedt += eloss_elec.getEloss(temp_electron, number_dens_electron, incident_speed);
    dedt += eloss_deut.getEloss(temp_deuteron, number_dens_deuteron, incident_speed);
    dedt += eloss_trit.getEloss(temp_triton, number_dens_triton, incident_speed);

    const double ref_dedt{1.92752975e+05 * incident_speed};
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(dedt, ref_dedt, 1.0e-8));
  }

  // Get dedt for E = 4*kT
  {
    double const incident_energy = 4.0e-3 * rtt_units::electronChargeSI * 1.0e+6; // J
    double const incident_speed =
        std::sqrt(2.0 * incident_energy / (alpha_mass * 1.0e-3)) * 1.0e-8 * 1.0e+2; // cm/shk

    double dedt{0.0};
    dedt += eloss_elec.getEloss(temp_electron, number_dens_electron, incident_speed);
    dedt += eloss_deut.getEloss(temp_deuteron, number_dens_deuteron, incident_speed);
    dedt += eloss_trit.getEloss(temp_triton, number_dens_triton, incident_speed);

    const double ref_dedt{2.98891620e+05 * incident_speed};
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(dedt, ref_dedt, 1.0e-8));
  }

  if (ut.numFails == 0)
    PASSMSG("Spitzer CPEloss test passes.");
  else
    FAILMSG("Spitzer CPEloss test fails.");

  return;
}

//------------------------------------------------------------------------------------------------//
void TR_test(rtt_dsxx::UnitTest &ut) {

  auto const ion_model_in(std::make_shared<rtt_cdi_cpeloss::Analytic_TR_Eloss_Model>(
      deuterium_zaid, deuterium_mass, alpha_zaid, alpha_mass));
  auto const electron_model_in(std::make_shared<rtt_cdi_cpeloss::Analytic_TR_Eloss_Model>(
      electron_zaid, electron_mass, alpha_zaid, alpha_mass));

  Analytic_CP_Eloss const eloss_ion_mod(ion_model_in, rtt_cdi::CPModelAngleCutoff::NONE);
  Analytic_CP_Eloss const eloss_ele_mod(electron_model_in, rtt_cdi::CPModelAngleCutoff::NONE);

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

  FAIL_IF_NOT(ion_model_in->targ_zaid() == eloss_ion_mod.targ_zaid());
  FAIL_IF_NOT(electron_model_in->targ_zaid() == eloss_ele_mod.targ_zaid());
  FAIL_IF_NOT(ion_model_in->proj_zaid() == eloss_ion_mod.proj_zaid());

  // Check that accessor returns the correct model angle cutoff
  FAIL_IF_NOT(eloss_ion_mod.getModelAngleCutoff() == rtt_cdi::CPModelAngleCutoff::NONE);

  // Get eloss values for some sample ion data:
  {
    double const T = 1.0;                       // keV
    double const ndens = 10.0 / deuterium_mass; // a / cc
    double const vel0 = 1.0;                    // cm / shk

    double const eloss_coeff = eloss_ion_mod.getEloss(T, ndens, vel0);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 9.927570982586801983e+05, 1.0e-8));
  }

  // Get eloss values for some sample electron data:
  {
    double const T = 1.0;                       // keV
    double const ndens = 1.e-2 / electron_mass; // el / cc
    double const vel0 = 1.0;                    // cm / shk

    double const eloss_coeff = eloss_ele_mod.getEloss(T, ndens, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 4.993807286690149340e+05, 1.0e-8));
  }

  // Check for large value when 2kT > E
  {
    double const T = 1.e2;                      // keV
    double const ndens = 10.0 / deuterium_mass; // g / cc
    double const vel0 = 1.0;                    // cm / shk

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
void LiPetrasso_test(rtt_dsxx::UnitTest &ut) {
  // Particle masses in [g]
  double const mass_electron = 9.1093837139e-28;
  double const mass_deuteron = 3.34449489528201e-24;
  double const mass_triton = 5.007356050158e-24;
  double const mass_alpha = 6.64647845996649e-24;

  // Temperature in keV
  double const temp_electron = 1.0;
  double const temp_deuteron = 1.0;
  double const temp_triton = 1.0;

  // Plasma Density g/cm^3
  double const dens_electron = 1.0;
  double const dens_deuteron = 0.5;
  double const dens_triton = 0.5;

  // Number densities
  double const number_dens_electron =
      dens_electron / (0.5 * (2.013553214 + 3.01550071913)) * rtt_units::AVOGADRO;
  double const number_dens_deuteron =
      dens_deuteron / (0.5 * (2.013553214 + 3.01550071913)) * rtt_units::AVOGADRO;
  double const number_dens_triton =
      dens_triton / (0.5 * (2.013553214 + 3.01550071913)) * rtt_units::AVOGADRO;

  // Build cpeloss model
  auto const elec_model(std::make_shared<rtt_cdi_cpeloss::Analytic_LiPet_Eloss_Model>(
      -1, mass_electron, 2004, mass_alpha));
  auto const deut_model(std::make_shared<rtt_cdi_cpeloss::Analytic_LiPet_Eloss_Model>(
      1002, mass_deuteron, 2004, mass_alpha));
  auto const trit_model(std::make_shared<rtt_cdi_cpeloss::Analytic_LiPet_Eloss_Model>(
      1003, mass_triton, 2004, mass_alpha));

  // Build CPEloss object using model:
  rtt_cdi_cpeloss::Analytic_CP_Eloss const eloss_elec(elec_model,
                                                      rtt_cdi::CPModelAngleCutoff::NONE);
  rtt_cdi_cpeloss::Analytic_CP_Eloss const eloss_deut(deut_model,
                                                      rtt_cdi::CPModelAngleCutoff::NONE);
  rtt_cdi_cpeloss::Analytic_CP_Eloss const eloss_trit(trit_model,
                                                      rtt_cdi::CPModelAngleCutoff::NONE);

  // Build background
  const std::map<int, double> zaid_num_dens = {
      {-1, number_dens_electron}, {1002, number_dens_deuteron}, {1003, number_dens_triton}};
  std::vector<double> const temperatures{temp_electron, temp_deuteron, temp_triton};
  const double debye_length(rtt_cdi::calculate_debye_length(zaid_num_dens, temperatures));

  // Check incident particle info
  FAIL_IF_NOT(eloss_elec.proj_zaid() == 2004);
  FAIL_IF_NOT(eloss_deut.proj_zaid() == 2004);
  FAIL_IF_NOT(eloss_trit.proj_zaid() == 2004);
  FAIL_IF_NOT(eloss_elec.targ_zaid() == -1);
  FAIL_IF_NOT(eloss_deut.targ_zaid() == 1002);
  FAIL_IF_NOT(eloss_trit.targ_zaid() == 1003);

  FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_elec.proj_mass(), mass_alpha, 1.0e-8));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_deut.proj_mass(), mass_alpha, 1.0e-8));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_trit.proj_mass(), mass_alpha, 1.0e-8));

  FAIL_IF_NOT(rtt_cdi::get_z_number(eloss_elec.get_Analytic_Model()->proj_zaid()) == 2);

  // Get dedt for 3.5 MeV alpha-particle
  {
    double const incident_energy = 3.5 * rtt_units::electronChargeSI * 1.0e+6; // J
    double const incident_speed =
        std::sqrt(2.0 * incident_energy / (mass_alpha * 1.0e-3)) * 1.0e-8 * 1.0e+2; // cm/shk

    double dedt{0.0};
    dedt += eloss_elec.getEloss(temp_electron, number_dens_electron, incident_speed, debye_length);
    dedt += eloss_deut.getEloss(temp_deuteron, number_dens_deuteron, incident_speed, debye_length);
    dedt += eloss_trit.getEloss(temp_triton, number_dens_triton, incident_speed, debye_length);

    const double ref_dedt{2.56071981e+05 * incident_speed};
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(dedt, ref_dedt, 1.0e-8));
  }

  // Get dedt for 1.0 MeV alpha-particle
  {
    double const incident_energy = 1.0 * rtt_units::electronChargeSI * 1.0e+6; // J
    double const incident_speed =
        std::sqrt(2.0 * incident_energy / (mass_alpha * 1.0e-3)) * 1.0e-8 * 1.0e+2; // cm/shk

    double dedt{0.0};
    dedt += eloss_elec.getEloss(temp_electron, number_dens_electron, incident_speed, debye_length);
    dedt += eloss_deut.getEloss(temp_deuteron, number_dens_deuteron, incident_speed, debye_length);
    dedt += eloss_trit.getEloss(temp_triton, number_dens_triton, incident_speed, debye_length);

    const double ref_dedt{1.61944659e+05 * incident_speed};
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(dedt, ref_dedt, 1.0e-8));
  }

  // Get dedt for E = 2kT
  {
    double const incident_energy = 2.0e-3 * rtt_units::electronChargeSI * 1.0e+6; // J
    double const incident_speed =
        std::sqrt(2.0 * incident_energy / (mass_alpha * 1.0e-3)) * 1.0e-8 * 1.0e+2; // cm/shk

    double dedt{0.0};
    dedt += eloss_elec.getEloss(temp_electron, number_dens_electron, incident_speed, debye_length);
    dedt += eloss_deut.getEloss(temp_deuteron, number_dens_deuteron, incident_speed, debye_length);
    dedt += eloss_trit.getEloss(temp_triton, number_dens_triton, incident_speed, debye_length);

    const double ref_dedt{3.05252070e+05 * incident_speed};
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(dedt, ref_dedt, 1.0e-8));
  }

#ifdef GSL_FOUND // Modified dedt requires calls to GSL
  std::cout << "[Running modified Li-Petrasso tests.]\n";
  // Get modified dedt for 3.5 MeV alpha-particle
  {
    double const incident_energy = 3.5 * rtt_units::electronChargeSI * 1.0e+6;
    double const incident_speed =
        std::sqrt(2.0 * incident_energy / (mass_alpha * 1.0e-3)) * 1.0e-8 * 1.0e+2; // cm/shk

    double dedt{0.0};
    dedt +=
        elec_model->getElossMod(temp_electron, number_dens_electron, incident_speed, debye_length);
    dedt +=
        deut_model->getElossMod(temp_deuteron, number_dens_deuteron, incident_speed, debye_length);
    dedt += trit_model->getElossMod(temp_triton, number_dens_triton, incident_speed, debye_length);

    const double ref_dedt{3.39685179e+05 * incident_speed};
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(dedt, ref_dedt, 1.0e-8));
  }

  // Get modified dedt for 1.0 MeV alpha-particle
  {
    double const incident_energy = 1.0 * rtt_units::electronChargeSI * 1.0e+6;
    double const incident_speed =
        std::sqrt(2.0 * incident_energy / (mass_alpha * 1.0e-3)) * 1.0e-8 * 1.0e+2; // cm/shk

    double dedt{0.0};
    dedt +=
        elec_model->getElossMod(temp_electron, number_dens_electron, incident_speed, debye_length);
    dedt +=
        deut_model->getElossMod(temp_deuteron, number_dens_deuteron, incident_speed, debye_length);
    dedt += trit_model->getElossMod(temp_triton, number_dens_triton, incident_speed, debye_length);

    const double ref_dedt{2.46587420e+05 * incident_speed};
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(dedt, ref_dedt, 1.0e-8));
  }

  // Get modified dedt for E = 2kT
  {
    double const incident_energy = 2.0e-3 * rtt_units::electronChargeSI * 1.0e+6;
    double const incident_speed =
        std::sqrt(2.0 * incident_energy / (mass_alpha * 1.0e-3)) * 1.0e-8 * 1.0e+2; // cm/shk

    double dedt{0.0};
    dedt +=
        elec_model->getElossMod(temp_electron, number_dens_electron, incident_speed, debye_length);
    dedt +=
        deut_model->getElossMod(temp_deuteron, number_dens_deuteron, incident_speed, debye_length);
    dedt += trit_model->getElossMod(temp_triton, number_dens_triton, incident_speed, debye_length);

    const double ref_dedt{3.20107656e+05 * incident_speed};
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(dedt, ref_dedt, 1.0e-8));
  }
#else
  std::cout << "[Skipping modified Li-Petrasso tests - GSL not found.]\n";
#endif

  if (ut.numFails == 0)
    PASSMSG("LiPetrasso CPEloss test passes.");
  else
    FAILMSG("LiPetrasso CPEloss test fails.");

  return;
}

//------------------------------------------------------------------------------------------------//
void Constant_test(rtt_dsxx::UnitTest &ut) {

  auto const zero_model_in(std::make_shared<rtt_cdi_cpeloss::Analytic_Constant_Eloss_Model>(
      deuterium_zaid, deuterium_mass, alpha_zaid, alpha_mass, 0.));
  auto const one_model_in(std::make_shared<rtt_cdi_cpeloss::Analytic_Constant_Eloss_Model>(
      deuterium_zaid, deuterium_mass, alpha_zaid, alpha_mass, 1.));

  Analytic_CP_Eloss const eloss_zero_mod(zero_model_in, rtt_cdi::CPModelAngleCutoff::NONE);
  Analytic_CP_Eloss const eloss_one_mod(one_model_in, rtt_cdi::CPModelAngleCutoff::NONE);

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

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    KP_alpha_test(ut);
    Spitzer_test(ut);
    TR_test(ut);
    LiPetrasso_test(ut);
    Constant_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstAnalytic_CP_Eloss.cc
//------------------------------------------------------------------------------------------------//
