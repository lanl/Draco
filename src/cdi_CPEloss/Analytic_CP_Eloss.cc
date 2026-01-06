//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/Analytic_CP_Eloss.cc
 * \author Kendra P. Long
 * \date   Fri Aug  2 14:28:08 2019
 * \brief  Analytic_CP_Eloss member definitions.
 * \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi_CPEloss/Analytic_CP_Eloss.hh"
#include "dsxx/Assert.hh"
#include "cdi/CPCommon.hh"
#include "cdi/CPEloss.hh"
#include <utility> // move

namespace rtt_cdi_cpeloss {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for an analytic gray eloss model.
 *
 * This constructor builds an eloss model defined by the rtt_cdi_analytic::Analytic_Eloss_Model
 * derived class argument.
 *
 * The reaction type for this instance of the class is determined by the rtt_cdi::Reaction argument.
 *
 * \param[in] model_in shared_ptr to a derived rtt_cdi_analytic::Analytic_Eloss_Model object
 * \param[in] model_angle_cutoff_in rtt_cdi::CPModelAngleCutoff the angle separating the stopping
 *                 power approximation from analog scattering
 */
Analytic_CP_Eloss::Analytic_CP_Eloss(SP_Analytic_Model model_in,
                                     rtt_cdi::CPModelAngleCutoff model_angle_cutoff_in)
    : rtt_cdi::CPEloss(rtt_cdi::CPModelType::ANALYTIC_ETYPE, model_angle_cutoff_in),
      analytic_model(std::move(model_in)) {
  Ensure(analytic_model);
}

//------------------------------------------------------------------------------------------------//
// ELOSS INTERFACE FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return a scalar eloss given a scalar temperature, density, and particle speed.
 *
 * Given a scalar temperature/density/speed, return an eloss for the reaction type specified by the
 * constructor.  The analytic eloss model is specified in the constructor (Analytic_CP_Eloss()).
 *
 * \param[in] temperature material temperature(s) in keV
 * \param[in] ndensity target number density/densities in atoms/cm^3
 * \param[in] v0 incident particle speed in cm/shk
 * \param[in] debye_length plasma Debye length in cm
 * \return eloss (time coefficient) in keV/shk^-1
 */
double Analytic_CP_Eloss::getEloss(const double temperature, const double ndensity, const double v0,
                                   const double debye_length) const {
  Require(v0 >= 0.0);

  // define return eloss
  double const eloss = analytic_model->calculate_eloss(temperature, ndensity, v0, debye_length);

  Ensure(eloss >= 0.0);
  return eloss;
}

} // namespace rtt_cdi_cpeloss

//------------------------------------------------------------------------------------------------//
// End cdi_CPEloss/Analytic_CP_Eloss.cc
//------------------------------------------------------------------------------------------------//
