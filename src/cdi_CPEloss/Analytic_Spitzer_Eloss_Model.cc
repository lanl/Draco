//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/Analytic_Spitzer_Eloss_Model.cc
 * \author Ben R. Ryan
 * \date   Feb 21 2020
 * \brief  Analytic_Spitzer_Eloss_Model implementation file.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Analytic_Spitzer_Eloss_Model.hh"
#include "ds++/Assert.hh"

namespace rtt_cdi_cpeloss {

/*! \brief Calculate the eloss in units of keV shk^-1; T given in keV,
 *         rho in g/cc, v0 in cm/shk
 *
 * CP energy loss class using standard free particle approach to stopping power
 * with Spitzer prescription for Coulomb logarithm.
 *
 * \param T material temperature in keV
 * \param rho material density in g/cm^3
 * \param v incident particle speed in cm/shk
 * \return eloss stopping power in keV shk^-1
 */
double Analytic_Spitzer_Eloss_Model::calculate_eloss(const double T, const double rho,
                                                     const double v) const {
  Require(T >= 0.0);
  Require(rho >= 0.0);
  Require(v >= 0.0);

  const double keV = 1e3 * pc.eV();
  const double eloss_unit_fac = 1.e-8 / keV;

  const double nt = rho / mt;
  const double Tt = T * keV / pc.k();
  const double vp = v * 1.e8;
  const double vt = sqrt(2. * pc.k() * Tt / mt);
  const double Ep = 1. / 2. * mp * vp * vp;

  const double prefac = 4. * pc.pi() * nt * qtabs * qtabs * qpabs * qpabs / (mt * vp * vp);

  // Spitzer approach: Argument of the Coulomb logarithm is the ratio of Debye
  // length to distance of closest approach
  const double debye_length = sqrt(pc.k() * Tt / (4. * pc.pi() * nt * qtabs * qtabs));
  const double vrel = sqrt(vp * vp + vt * vt);
  const double mreduced = (mt * mp) / (mt + mp);
  const double bmin = qtabs * qpabs / (mreduced * vrel * vrel);
  const double coulomb_log = log(std::max<double>(exp(1.), debye_length / bmin));

  double eloss = prefac * coulomb_log;

  if (Ep < 2. * pc.k() * Tt) {
    return rtt_cdi::constants::max_eloss;
  }

  return vp * eloss * eloss_unit_fac;
}

} // namespace rtt_cdi_cpeloss

//------------------------------------------------------------------------------------------------//
// End cdi_CPEloss/Analytic_Spitzer_Eloss_Model.cc
//------------------------------------------------------------------------------------------------//
