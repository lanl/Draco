//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/Analytic_LiPet_Eloss_Model.cc
 * \author Edward Norris and Kyle S. Beling
 * \date   May 21, 2025
 * \brief  Analytic_LiPet_Eloss_Model class definition.
 * \note   Copyright (C) 2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi_CPEloss/Analytic_LiPet_Eloss_Model.hh"
#include "dsxx/Assert.hh"
#include "cdi/CPCommon.hh"
#include "cdi_CPEloss/config.h" // Defines GSL_FOUND
#include "units/MathConstants.hh"
#include "units/PhysicalConstantsSI.hh"

#include <cmath>
#ifdef GSL_FOUND
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#endif

namespace rtt_cdi_cpeloss {
//------------------------------------------------------------------------------------------------//
/*! \brief Calculate the eloss in units of [keV/shk] for a specific target
 * 
 * This wraps calculate_eloss_modified() in a way that looks similar to the Analytic_CP_Eloss
 * getEloss() call. It's currently implemented for comparison purposes.
 *
 * \param ndensity number density of target in [cm^-3]
 * \param temperature temperature of target particle in [keV]
 * \param v0 incident speed in [cm/shk]
 * \param debye_length material Debye length in [cm]
 * \return stopping power in [keV/shk]
 */
double Analytic_LiPet_Eloss_Model::getElossMod(const double temperature, const double ndensity,
                                               const double v0, const double debye_length) const {
  Require(v0 >= 0.0);

  // define return eloss
  double const eloss = calculate_eloss_modified(temperature, ndensity, v0, debye_length);

  Ensure(eloss >= 0.0);
  return eloss;
}

//------------------------------------------------------------------------------------------------//
/*! \brief Calculate the eloss in units of [keV/shk] for a specific target
 * 
 * This stopping power model is taken from the original 1993 Li-Petrasso model plus the 2015 erratum. 
 * SI units are used to internally maintain consistency and comparability with Edward's
 * implementation in the CPUtilities library.
 *
 * \param temperature_in temperature of target particle in [keV]
 * \param number_density_in number density of target in [cm^-3]
 * \param incident_speed_in incident speed in [cm/shk]
 * \param debye_length_in material Debye length in [cm]
 * \return stopping power in [keV/shk]
 */
double Analytic_LiPet_Eloss_Model::calculate_eloss(const double temperature_in,
                                                   const double number_density_in,
                                                   const double incident_speed_in,
                                                   const double debye_length_in) const {
  // Convert number density, temperature, debye
  double const number_density = number_density_in * 1.0e+6;             // [cm^-3 -> m^-3]
  double const temperature = temperature_in * rtt_units::EV2K * 1.0e+3; // [keV -> K]
  double const debye_length = debye_length_in * 1.0e-2;                 // [cm -> m]

  // Get z number:
  const double incident_z{static_cast<double>(rtt_cdi::get_z_number(zaidp))};

  // Convert speed
  double const incident_speed = incident_speed_in * 1.0e+8 * 1.0e-2; // [cm/shk -> m/s]

  // Convert masses (g-> kg):
  const double target_mass{mt * 1.0e-3};
  const double incident_mass{mp * 1.0e-3};

  // Compute Coulomb log
  double const mr = incident_mass * target_mass / (incident_mass + target_mass);
  double const speed2 = std::pow(incident_speed, 2);
  double const u_rel =
      sqrt(speed2) + sqrt(2.0 * rtt_units::boltzmannSI * temperature / target_mass);
  double const t12 = std::pow(rtt_units::h_barSI / (2.0 * mr * u_rel), 2);
  double const p_perp = std::pow(
      incident_z * rtt_units::electronChargeSI * rtt_units::electronChargeSI /
          (4.0 * rtt_units::PI * rtt_units::permittivityOfFreeSpaceSI * mr * u_rel * u_rel),
      2);
  double const p_min = std::sqrt(p_perp + t12);
  double const lambda_D = debye_length / p_min;
  double const coulomb_log = 0.5 * std::log(1.0 + lambda_D * lambda_D);

  double const omega2 = number_density * rtt_units::electronChargeSI * rtt_units::electronChargeSI /
                        (target_mass * rtt_units::permittivityOfFreeSpaceSI);

  double const x = speed2 * target_mass / (2.0 * rtt_units::boltzmannSI * temperature);
  double const xcol = 2.0 * x; // 2015 correction

  double const maxwell_integral =
      std::erf(std::sqrt(x)) - 2.0 / std::sqrt(rtt_units::PI) * std::exp(-x) * std::sqrt(x);
  double const maxwell_func = 2.0 / std::sqrt(rtt_units::PI) * std::exp(-x) * std::sqrt(x);

  double const corrector = 1.0 / coulomb_log * (maxwell_integral + maxwell_func);
  double const G = maxwell_integral - target_mass / incident_mass * (maxwell_func - corrector);

  double theta = 0.0;
  if (xcol > 1.0) {
    theta = 1.0;
  }

  double dedx = incident_z * incident_z * rtt_units::electronChargeSI *
                rtt_units::electronChargeSI /
                (4.0 * rtt_units::PI * rtt_units::permittivityOfFreeSpaceSI) * omega2 / speed2 *
                (G * coulomb_log + theta * log(1.123 * sqrt(xcol)));

  // Convert units J/m -> keV/cm
  dedx /= (rtt_units::electronChargeSI * 1.0e+3 * 1.0e+2);
  // On return, convert from dedx [keV/cm] -> dedt [keV/shk]
  return dedx * incident_speed_in;
}

//------------------------------------------------------------------------------------------------//
/*! \brief Calculate the modified Li-Petrasso (Zylstra-Li-Petrasso) stopping power eloss in units 
 *         of [keV/shk]
 * 
 * This stopping power model is taken from Zylstra et al. 2019, "Modified parameterization of the 
 * Li-Petrasso charged particle stopping power theory".
 *
 * \param temperature_in temperature of target particle in [keV]
 * \param number_density_in number density of target in [cm^-3]
 * \param incident_speed_in incident speed in [cm/shk]
 * \param debye_length_in material Debye length in [cm]
 * \return stopping power in [keV/shk]
 */
#ifdef GSL_FOUND
double Analytic_LiPet_Eloss_Model::calculate_eloss_modified(const double temperature_in,
                                                            const double number_density_in,
                                                            const double incident_speed_in,
                                                            const double debye_length_in) const {
  // Convert number density, temperature, debye
  double const number_density = number_density_in * 1.0e+6;             // [cm^-3 -> m^-3]
  double const temperature = temperature_in * rtt_units::EV2K * 1.0e+3; // [keV -> K]
  double const debye_length = debye_length_in * 1.0e-2;                 // [cm -> m]

  // Get z number:
  const double incident_z{static_cast<double>(rtt_cdi::get_z_number(zaidp))};

  // Convert speed
  double const incident_speed = incident_speed_in * 1.0e+8 * 1.0e-2; // [cm/shk -> m/s]

  // Get masses (convert from g -> kg):
  const double target_mass{mt * 1.0e-3};
  const double incident_mass{mp * 1.0e-3};

  // Modified temperature for quantum degeneracy
  double const theta =
      temperature * (std::pow(8.0 * rtt_units::PI / (3.0 * number_density), 2.0 / 3.0) * 2.0 *
                     rtt_units::electronMassSI * rtt_units::boltzmannSI /
                     (rtt_units::planckSI * rtt_units::planckSI));
  double const param =
      rtt_units::boltzmannSI * temperature *
      (-1.5 * std::log(theta) + std::log(4.0 / (3.0 * std::sqrt(rtt_units::PI))) +
       (0.25054 * std::pow(theta, -1.858) +
        0.027 * std::pow(theta, -1.858 / 2.0) / (1.0 + 0.25054 * std::pow(theta, -0.868))));
  double const T_eff =
      temperature * gsl_sf_fermi_dirac_3half(param) / gsl_sf_fermi_dirac_half(param);

  // Compute Coulomb log
  double const mr = incident_mass * target_mass / (incident_mass + target_mass);
  double const speed2 = std::pow(incident_speed, 2);
  double const u_rel = sqrt(speed2) + sqrt(2.0 * rtt_units::boltzmannSI * T_eff / target_mass);
  double const t12 = std::pow(rtt_units::h_barSI / (2.0 * mr * u_rel), 2);
  double const p_perp = std::pow(
      incident_z * rtt_units::electronChargeSI * rtt_units::electronChargeSI /
          (4.0 * rtt_units::PI * rtt_units::permittivityOfFreeSpaceSI * mr * u_rel * u_rel),
      2);
  double const p_min = std::sqrt(p_perp + t12);
  double const lambda_D = debye_length / p_min;
  double const coulomb_log = 0.5 * std::log(1.0 + lambda_D * lambda_D);

  double const omega2 = number_density * rtt_units::electronChargeSI * rtt_units::electronChargeSI /
                        (target_mass * rtt_units::permittivityOfFreeSpaceSI);

  double const x = speed2 * target_mass / (2.0 * rtt_units::boltzmannSI * T_eff);

  double const maxwell_integral =
      std::erf(std::sqrt(x)) - 2.0 / std::sqrt(rtt_units::PI) * std::exp(-x) * std::sqrt(x);
  double const maxwell_func = 2.0 / std::sqrt(rtt_units::PI) * std::exp(-x) * std::sqrt(x);

  double const corrector = 1.0 / coulomb_log * (maxwell_integral + maxwell_func);
  double const G = maxwell_integral - target_mass / incident_mass * (maxwell_func - corrector);

  double const factor = std::sqrt(omega2 / speed2) * debye_length;
  double const parameterization = factor * gsl_sf_bessel_K0(factor) * gsl_sf_bessel_K1(factor);

  double dedx = incident_z * incident_z * rtt_units::electronChargeSI *
                rtt_units::electronChargeSI /
                (4.0 * rtt_units::PI * rtt_units::permittivityOfFreeSpaceSI) * omega2 / speed2 *
                (G * coulomb_log + parameterization);

  // Convert units J/m -> keV/cm
  dedx /= (rtt_units::electronChargeSI * 1.0e+3 * 1.0e+2);

  // On return, convert from dedx [keV/cm] -> dedt [keV/shk]
  return dedx * incident_speed_in;
}
#else
double Analytic_LiPet_Eloss_Model::calculate_eloss_modified(
    [[maybe_unused]] const double number_density_in, [[maybe_unused]] const double temperature_in,
    [[maybe_unused]] const double incident_speed_in,
    [[maybe_unused]] const double debye_length_in) const {
  Insist(false, "Analytic_LiPet_Eloss_Model::calculate_eloss_modified requires GSL!");
}
#endif

} // namespace rtt_cdi_cpeloss

//------------------------------------------------------------------------------------------------//
// End cdi_CPEloss/Analytic_LiPet_Eloss_Model.cc
//------------------------------------------------------------------------------------------------//
