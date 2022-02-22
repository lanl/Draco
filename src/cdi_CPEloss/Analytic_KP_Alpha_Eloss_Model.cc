//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/Analytic_KP_Alpha_Eloss_Model.cc
 * \author Kendra P. Long, Ben R. Ryan
 * \date   Feb 21 2020
 * \brief  Analytic_KP_Alpha_Eloss_Model implementation file.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Analytic_KP_Alpha_Eloss_Model.hh"
#include "ds++/Assert.hh"

namespace rtt_cdi_cpeloss {

/*! \brief Calculate the eloss in units of shk^-1; T given in keV, rho in g/cc,
 *         v0 in cm/shk
 *
 * The constants in this formula come directly from the fit in  Eq. (2) of
 * Kirkpatrick, R. C. and Wheeler, J. A. (1981).
 * ``The Physics of DT Ignition In Small Fusion Targets.''
 * Nuclear Fusion, 21(3):389–401.
 *
 * These constants are ONLY valid for alpha energy loss in DT gas.
 *
 * \param T material temperature in keV
 * \param rho material density in g/cm^3
 * \param v incident particle speed in cm/shk
 * \return eloss stopping power in keV shk^-1
 *
 */
double Analytic_KP_Alpha_Eloss_Model::calculate_eloss(const double T, const double rho,
                                                      const double v) const {
  Require(T >= 0.0);
  Require(rho >= 0.0);
  Require(v >= 0.0);

  double range = 0.03 * T * (1.0 - 0.24 * std::log(1.0 + T)) *
                 (1.0 + 0.37 * std::log((1.0 + rho) / (1.0 + 0.01 * T * T)));

  // Exponent of energy deposition term, sans minus sign and delta_t:
  double eloss = rho * v / range;

  // This analytic model can return negative elosses in some rho-T regimes;
  // catch these and simply return a large, positive number instead
  // (equivlant to a very small, positive value of "range" in the above formula)
  return (eloss >= 0.0) ? eloss : 1.0e25;
}

} // namespace rtt_cdi_cpeloss

//------------------------------------------------------------------------------------------------//
// End cdi_CPEloss/Analytic_KP_Alpha_Eloss_Model.cc
//------------------------------------------------------------------------------------------------//
