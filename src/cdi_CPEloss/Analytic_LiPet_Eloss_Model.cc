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
#include "units/MathConstants.hh"
#include "units/PhysicalConstantsSI.hh"

#include <cmath>
#include <cstddef>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <vector>

namespace rtt_cdi_cpeloss {

//------------------------------------------------------------------------------------------------//
/*! \brief Calculate the total eloss in units of [MeV/cm] for some background composition
 *
 * This stopping power model is taken from the original 1993 Li-Petrasso model plus the 2015 erratum. 
 * SI units are used to maintain consistency and comparability with Edward's implementation in the CPUtilities library.
 *
 * \param incident_speed incident speed of the ion in [cm/shk]
 * \param zaid_list list of ion ZAID's
 * \param number_densities list of ion number densities in [cm^-3]
 * \param temperatures list of ion temperatures in [keV]
 * \param masses list of masses in [kg]
 * \return stopping power in [MeV/cm]
 */
double Analytic_LiPet_Eloss_Model::getEloss(const double incident_speed,
                                            const std::vector<int> &zaid_list,
                                            const std::vector<double> &number_densities,
                                            const std::vector<double> &temperatures,
                                            const std::vector<double> &masses) const {

  // Check that the correct number of entries are provided
  Require(zaid_list.size() == number_densities.size());
  Require(zaid_list.size() == temperatures.size());
  Require(zaid_list.size() == masses.size());

  // Convert units
  double const incident_speed_mps = incident_speed * 1.0e+8 * 1.0e-2; // [cm/shk -> m/s]

  // Get Debye length in [m]
  double const debye_length = calculate_debye_length(zaid_list, number_densities, temperatures);

  // Compute total stopping power
  double dedx = 0.0;
  for (size_t i = 0; i < zaid_list.size(); ++i) {
    double const number_density = number_densities[i] * 1.0e+6;            // [cm^-3 -> m^-3]
    double const temperature = temperatures[i] * rtt_units::EV2K * 1.0e+3; // [keV -> K]
    double const target_mass = masses[i];

    dedx +=
        calculate_eloss(incident_speed_mps, number_density, temperature, target_mass, debye_length);
  }

  // Convert units J/m -> MeV/cm
  dedx /= (rtt_units::electronChargeSI * 1.0e+6 * 1.0e+2);
  return dedx;
}

//------------------------------------------------------------------------------------------------//
/*! \brief Calculate the modified total eloss in units of [MeV/cm] for some background composition
 *
 * This stopping power model is taken from Zylstra et al. 2019, "Modified parameterization of the 
 * Li-Petrasso charged particle stopping power theory".
 *
 * \param incident_speed incident speed of the ion in [cm/shk]
 * \param zaid_list list of ion ZAID's
 * \param number_densities list of ion number densities in [cm^-3]
 * \param temperatures list of ion temperatures in [keV]
 * \param masses list of masses in [kg]
 * \return stopping power in [MeV/cm]
 */
double Analytic_LiPet_Eloss_Model::getElossMod(const double incident_speed,
                                               const std::vector<int> &zaid_list,
                                               const std::vector<double> &number_densities,
                                               const std::vector<double> &temperatures,
                                               const std::vector<double> &masses) const {

  // Check that the correct number of entries are provided
  Require(zaid_list.size() == number_densities.size());
  Require(zaid_list.size() == temperatures.size());
  Require(zaid_list.size() == masses.size());

  // Convert units
  double const incident_speed_mps = incident_speed * 1.0e+8 * 1.0e-2; // [cm/shk -> m/s]

  // Get Debye length in [m]
  double const debye_length = calculate_debye_length(zaid_list, number_densities, temperatures);

  // Compute total stopping power
  double dedx = 0.0;
  for (size_t i = 0; i < zaid_list.size(); ++i) {
    double const number_density = number_densities[i] * 1.0e+6;            // [cm^-3 -> m^-3]
    double const temperature = temperatures[i] * rtt_units::EV2K * 1.0e+3; // [keV -> K]
    double const target_mass = masses[i];

    dedx += calculate_eloss_modified(incident_speed_mps, number_density, temperature, target_mass,
                                     debye_length);
  }

  // Convert units J/m -> MeV/cm
  dedx /= (rtt_units::electronChargeSI * 1.0e+6 * 1.0e+2);
  return dedx;
}

//------------------------------------------------------------------------------------------------//
/*! \brief Calculate the Debye length in [m]
 *
 * lambda_D = sum( n * z^2 / kT )
 *
 * \param zaid_list list of ion ZAID's
 * \param number_densities list of ion number densities in [1/cm^3]
 * \param temperatures list of ion temperatures in [keV]
 * \return Debye length in [m]
 */
double
Analytic_LiPet_Eloss_Model::calculate_debye_length(const std::vector<int> &zaid_list,
                                                   const std::vector<double> &number_densities,
                                                   const std::vector<double> &temperatures) const {

  double debye_length = 0.0;
  for (size_t i = 0; i < zaid_list.size(); ++i) {
    int z_number = 0;
    if (zaid_list[i] == -1) {
      z_number = 1;
    } else {
      z_number = (zaid_list[i] - zaid_list[i] % 1000) / 1000;
    }

    // Density [cm^-3 -> m^-3]
    double const density = number_densities[i] * 1.0e+6;

    // Temperature [keV -> Kelvin]
    double const temperature = temperatures[i] * rtt_units::EV2K * 1.0e+3;

    debye_length += density * rtt_units::electronChargeSI * rtt_units::electronChargeSI * z_number *
                    z_number / temperature;
  }
  debye_length =
      std::sqrt(rtt_units::permittivityOfFreeSpaceSI * rtt_units::boltzmannSI / debye_length);

  return debye_length;
}

//------------------------------------------------------------------------------------------------//
/*! \brief Calculate the eloss in units of [J/m] for a specific target
 * 
 * This stopping power model is taken from the original 1993 Li-Petrasso model plus the 2015 erratum. 
 * SI units are used to maintain consistancy and comparability with Edward's implimentation in the CPUtilities library.
 *
 * \param incident_speed incident energy in [m/s]
 * \param number_density number density of target in [m^-3]
 * \param temperature temperature of target particle in [Kelvin]
 * \param target_mass mass of target particle in [kg]
 * \param debye_length material Debye length in [m]
 * \return stopping power in [J/m]
 */
double Analytic_LiPet_Eloss_Model::calculate_eloss(const double incident_speed,
                                                   const double number_density,
                                                   const double temperature,
                                                   const double target_mass,
                                                   const double debye_length) const {

  // Compute Coulomb log
  double const mr = incident_mass_ * target_mass / (incident_mass_ + target_mass);
  double const speed2 = std::pow(incident_speed, 2);
  double const u_rel =
      sqrt(speed2) + sqrt(2.0 * rtt_units::boltzmannSI * temperature / target_mass);
  double const t12 = std::pow(rtt_units::h_barSI / (2.0 * mr * u_rel), 2);
  double const p_perp = std::pow(
      incident_z_ * rtt_units::electronChargeSI * rtt_units::electronChargeSI /
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
  double const G = maxwell_integral - target_mass / incident_mass_ * (maxwell_func - corrector);

  double theta = 0.0;
  if (xcol > 1.0) {
    theta = 1.0;
  }

  double const dedx = incident_z_ * incident_z_ * rtt_units::electronChargeSI *
                      rtt_units::electronChargeSI /
                      (4.0 * rtt_units::PI * rtt_units::permittivityOfFreeSpaceSI) * omega2 /
                      speed2 * (G * coulomb_log + theta * log(1.123 * sqrt(xcol)));

  return dedx;
}

//------------------------------------------------------------------------------------------------//
/*! \brief Calculate the modified Li-Petrasso (Zylstra-Li-Petrasso) stopping power eloss in units of [J/m]
 * 
 * This stopping power model is taken from Zylstra et al. 2019, "Modified parameterization of the 
 * Li-Petrasso charged particle stopping power theory".
 *
 * \param incident_speed incident energy in [m/s]
 * \param number_density number density of target in [m^-3]
 * \param temperature temperature of target particle in [Kelvin]
 * \param target_mass mass of target particle in [kg]
 * \param debye_length material Debye length in [m]
 * \return stopping power in [J/m]
 */
double Analytic_LiPet_Eloss_Model::calculate_eloss_modified(const double incident_speed,
                                                            const double number_density,
                                                            const double temperature,
                                                            const double target_mass,
                                                            const double debye_length) const {

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
  double const mr = incident_mass_ * target_mass / (incident_mass_ + target_mass);
  double const speed2 = std::pow(incident_speed, 2);
  double const u_rel = sqrt(speed2) + sqrt(2.0 * rtt_units::boltzmannSI * T_eff / target_mass);
  double const t12 = std::pow(rtt_units::h_barSI / (2.0 * mr * u_rel), 2);
  double const p_perp = std::pow(
      incident_z_ * rtt_units::electronChargeSI * rtt_units::electronChargeSI /
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
  double const G = maxwell_integral - target_mass / incident_mass_ * (maxwell_func - corrector);

  double const factor = std::sqrt(omega2 / speed2) * debye_length;
  double const parameterization = factor * gsl_sf_bessel_K0(factor) * gsl_sf_bessel_K1(factor);

  double const dedx = incident_z_ * incident_z_ * rtt_units::electronChargeSI *
                      rtt_units::electronChargeSI /
                      (4.0 * rtt_units::PI * rtt_units::permittivityOfFreeSpaceSI) * omega2 /
                      speed2 * (G * coulomb_log + parameterization);

  return dedx;
}

} // namespace rtt_cdi_cpeloss

//------------------------------------------------------------------------------------------------//
// End cdi_CPEloss/Analytic_LiPet_Eloss_Model.cc
//------------------------------------------------------------------------------------------------//
