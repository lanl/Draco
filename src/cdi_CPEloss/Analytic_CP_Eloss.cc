//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/Analytic_CP_Eloss.cc
 * \author Kendra P. Long
 * \date   Fri Aug  2 14:28:08 2019
 * \brief  Analytic_CP_Eloss member definitions.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Analytic_CP_Eloss.hh"
#include "ds++/Assert.hh"

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
 * \param[in] target_in int32_t target particle
 * \param[in] projectile_in int32_t particle being transported
 * \param[in] model_angle_cutoff_in rtt_cdi::CPModelAngleCutoff the angle separating the stopping
 *                 power approximation from analog scattering
 */
Analytic_CP_Eloss::Analytic_CP_Eloss(SP_Analytic_Model model_in,
                                     rtt_cdi::CParticle target_in, // NOLINT
                                     rtt_cdi::CParticle projectile_in,
                                     rtt_cdi::CPModelAngleCutoff model_angle_cutoff_in)
    : rtt_cdi::CPEloss(target_in, projectile_in, rtt_cdi::CPModelType::ANALYTIC_ETYPE,
                       model_angle_cutoff_in),
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
 * \param[in] temperature material temperature in keV
 * \param[in] density material density in g/cm^3
 * \param[in] v0 incident particle speed in cm/shk
 * \return eloss (time coefficient) in shk^-1
 */
double Analytic_CP_Eloss::getEloss(double temperature, double density, double v0) const {
  Require(temperature >= 0.0);
  Require(density >= 0.0);
  Require(v0 >= 0.0);

  // define return eloss
  double eloss = analytic_model->calculate_eloss(temperature, density, v0);

  Ensure(eloss >= 0.0);
  return eloss;
}

} // namespace rtt_cdi_cpeloss

//------------------------------------------------------------------------------------------------//
// End cdi_CPEloss/Analytic_CP_Eloss.cc
//------------------------------------------------------------------------------------------------//
