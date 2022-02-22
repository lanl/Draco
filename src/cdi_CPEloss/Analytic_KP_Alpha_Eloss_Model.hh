//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/Analytic_KP_Alpha_Eloss_Model.hh
 * \author Kendra P. Long, Ben R. Ryan
 * \date   Feb 21 2020
 * \brief  Analytic_KP_Alpha_Eloss_Model class definition.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_CPEloss_Analytic_KP_Alpha_Eloss_Model_hh
#define rtt_cdi_CPEloss_Analytic_KP_Alpha_Eloss_Model_hh

#include "Analytic_Eloss_Model.hh"
#include "ds++/Assert.hh"

namespace rtt_cdi_cpeloss {
//================================================================================================//
/*!
 * \class Analytic_KP_Alpha_Eloss_Model
 * \brief Analytic_KP_Alpha_Eloss_Model derived class.
 *
 * CP energy loss class using analytic Kirkpatrick model for alpha particles in
 * DT. Energy loss rates based on the range fit calculated in:
 *
 * Kirkpatrick, R. C. and Wheeler, J. A. (1981).
 * ``The Physics of DT Ignition In Small Fusion Targets.''
 * Nuclear Fusion, 21(3):389–401.
 *
 * Equation (2) gives the range formula. We convert this to an energy loss rate
 * per unit time for ease of use in transport.
 */
//================================================================================================//

class Analytic_KP_Alpha_Eloss_Model : public Analytic_Eloss_Model {
public:
  //! Constructor
  Analytic_KP_Alpha_Eloss_Model(const rtt_cdi::CParticle &target,
                                const rtt_cdi::CParticle &projectile)
      : Analytic_Eloss_Model(target, projectile) {
    Require(projectile.get_zaid() == 2004);
  }

  //! Calculate the eloss rate in units of shk^-1;
  //! T given in keV, rho in g/cc, v0 in cm/shk
  double calculate_eloss(const double T, const double rho, const double v0) const override;
};

} // namespace rtt_cdi_cpeloss

#endif // rtt_cdi_CPEloss_Analytic_KP_Alpha_Eloss_Model_hh

//------------------------------------------------------------------------------------------------//
// End cdi_CPEloss/Analytic_KP_Alpha_Eloss_Model.hh
//------------------------------------------------------------------------------------------------//
