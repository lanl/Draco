//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/Analytic_KP_Alpha_Eloss_Model.cc
 * \brief  Analytic_KP_Alpha_Eloss_Model implementation file.
 * \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi_CPEloss/Analytic_KP_Alpha_Eloss_Model.hh"
#include "dsxx/Assert.hh"
#include <cmath> //log

namespace rtt_cdi_cpeloss {

/*! \brief Calculate the eloss in units of shk^-1; T given in keV, ndens in a/cc, v0 in cm/shk
 *
 * The constants in this formula come directly from the fit in  Eq. (2) of
 *
 * * Kirkpatrick, R. C. and Wheeler, J. A. (1981). "The Physics of DT Ignition In Small Fusion
 *   Targets." Nuclear Fusion, 21(3):389–401.
 *
 * These constants are ONLY valid for alpha energy loss in DT gas.
 *
 * \param T material temperature in keV
 * \param ndens number density in atoms/cm^3
 * \param v incident particle speed in cm/shk
 * \param debye_length plasma Debye length (cm, unused)
 * \return eloss stopping power in keV shk^-1
 */
double
Analytic_KP_Alpha_Eloss_Model::calculate_eloss(const double T, const double ndens, const double v,
                                               [[maybe_unused]] const double debye_length) const {
  Require(T >= 0.0);
  Require(ndens >= 0.0);
  Require(v >= 0.0);

  // g/cc = at/cc * g/at
  double const rho{ndens * mt};

  double const range = 0.03 * T * (1.0 - 0.24 * std::log(1.0 + T)) *
                       (1.0 + 0.37 * std::log((1.0 + rho) / (1.0 + 0.01 * T * T)));

  // Exponent of energy deposition term, sans minus sign and delta_t:
  double const eloss = rho * v / range;

  // This analytic model can return negative elosses in some rho-T regimes; catch these and simply
  // return a large, positive number instead (equivalent to a very small, positive value of "range"
  // in the above formula)
  return (eloss >= 0.0) ? eloss : 1.0e25;
}

} // namespace rtt_cdi_cpeloss

//------------------------------------------------------------------------------------------------//
// End cdi_CPEloss/Analytic_KP_Alpha_Eloss_Model.cc
//------------------------------------------------------------------------------------------------//
