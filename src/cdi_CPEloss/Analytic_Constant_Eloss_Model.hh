//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/Analytic_Constant_Eloss_Model.hh
 * \author Ben R. Ryan
 * \date   Feb 21 2020
 * \brief  Analytic_Constant_Eloss_Model class definition.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_CPEloss_Analytic_Constant_Eloss_Model_hh
#define rtt_cdi_CPEloss_Analytic_Constant_Eloss_Model_hh

#include "Analytic_Eloss_Model.hh"
#include <algorithm>

namespace rtt_cdi_cpeloss {
//================================================================================================//
/*!
 * \class Analytic_Constant_Eloss_Model
 * \brief Analytic_Constant_Eloss_Model derived class.
 *
 * CP energy loss class using constant value provided to constructor.
 */
//================================================================================================//

class Analytic_Constant_Eloss_Model : public Analytic_Eloss_Model {
public:
  //! Constructor
  Analytic_Constant_Eloss_Model(const rtt_cdi::CParticle &target,
                                const rtt_cdi::CParticle &projectile, const double eloss_value_in)
      : Analytic_Eloss_Model(target, projectile), eloss_value(eloss_value_in) {}

  //! Calculate the eloss rate in units of shk^-1;
  //! T given in keV, rho in g/cc, v0 in cm/shk
  double calculate_eloss(const double /*T*/, const double /*rho*/,
                         const double /*v0*/) const override {
    return eloss_value;
  }

private:
  double eloss_value;
};

} // namespace rtt_cdi_cpeloss

#endif // rtt_cdi_CPEloss_Analytic_Constant_Eloss_Model_hh

//------------------------------------------------------------------------------------------------//
// End cdi_CPEloss/Analytic_Constant_Eloss_Model.hh
//------------------------------------------------------------------------------------------------//
