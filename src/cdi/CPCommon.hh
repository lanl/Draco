//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/CPCommon.hh
 * \author Kendra Long
 * \date   Mon Jan 19 13:41:01 2001
 * \brief  Datatypes needed in CPEloss classes
 * \note   Copyright (C) 2019-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_CPCommon_hh
#define rtt_cdi_CPCommon_hh

#include "dsxx/config.h" // defines [[nodiscard]]
#include "units/MathConstants.hh"
#include "units/PhysicalConstantsSI.hh"
#include <cstdint>
#include <map>
#include <vector>

namespace rtt_cdi {

//------------------------------------------------------------------------------------------------//
/*! \brief Retrieve the charge for a given ZAID
 *
 * \param zaid ion ZAID
 * \return absolute charge (aka "z number")
 */
static int get_z_number(int zaid) { return (zaid == -1) ? 1 : (zaid) / 1000; }
//------------------------------------------------------------------------------------------------//
/*! \brief Calculate the Debye length in [cm]
 *
 * lambda_D = sum( n * z^2 / kT )
 *
 * \param zaid_num_dens map of ion ZAIDs / number densities [1/cm^3]
 * \param temperatures list of ion temperatures in [keV]
 * \return Debye length in [cm]
 */
[[maybe_unused]] static double calculate_debye_length(const std::map<int, double> &zaid_num_dens,
                                                      const std::vector<double> &temperatures) {
  Check(zaid_num_dens.size() == temperatures.size());
  double debye_length = 0.0;
  size_t idx{0};
  for (auto &zd : zaid_num_dens) {
    // Skip isos with zero number density.
    if (zd.second > 0.0) {
      int const z_number{get_z_number(zd.first)};

      // Density [cm^-3 -> m^-3]
      double const density = zd.second * 1.0e+6;

      // Temperature [keV -> Kelvin]
      double const temperature = temperatures[idx] * rtt_units::EV2K * 1.0e+3;

      debye_length += density * rtt_units::electronChargeSI * rtt_units::electronChargeSI *
                      z_number * z_number / temperature;
    }
    idx++;
  }
  // Includes m -> cm conversion
  debye_length =
      std::sqrt(rtt_units::permittivityOfFreeSpaceSI * rtt_units::boltzmannSI / debye_length) *
      1.0e+2;

  return debye_length;
}
//------------------------------------------------------------------------------------------------//
/*! \brief Calculate the Debye length in [cm] for a single ion temperature
 *         This version wraps the full species-specific-temperature call.
 *
 * lambda_D = sum( n * z^2 / kT )
 *
 * \param zaid_num_dens map of ion ZAIDs / number densities [1/cm^3]
 * \param Ti ion temperature in [keV]
 * \param Te electron temperature in [keV]
 * \return Debye length in [cm]
 */
[[maybe_unused]] static double calculate_debye_length(const std::map<int, double> &zaid_num_dens,
                                                      const double Ti, const double Te) {
  // Set up temperature array:
  std::vector<double> temperatures(zaid_num_dens.size(), 0.0);
  size_t idx{0};
  for (auto zd : zaid_num_dens) {
    temperatures[idx] = (zd.first == -1) ? Te : Ti;
    idx++;
  }

  return (calculate_debye_length(zaid_num_dens, temperatures));
}

namespace constants {
//! Number of models contained in rtt_cdi::CPModel.
unsigned int const num_CPModels(3);

//! "Effectively infinite" stopping power (keV sh^-1)
constexpr double max_eloss = 1.e15;
} // namespace constants

enum class CPModelAngleCutoff {
  NONE = 0,  /*!< no angle cutoff */
  TNBURN = 1 /*!< cutoff angle formula from "TN Burn Project: Minimum Impact parameter Scattering
                  Angle Cutoff/Nuclear Interactions */
};

enum class CPModelType {
  UNASSIGNED_ETYPE = 0, /*!< unassigned type; used as a placeholder before deciding type */
  ANALYTIC_ETYPE = 1,   /*!< an Analytic eloss model */
  TABULAR_ETYPE = 2     /*!< a Tabular eloss model */
};

} // namespace rtt_cdi

#endif

//------------------------------------------------------------------------------------------------//
// End cdi/CPCommon.hh
//------------------------------------------------------------------------------------------------//
