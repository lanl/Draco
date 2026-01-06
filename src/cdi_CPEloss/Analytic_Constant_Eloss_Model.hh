//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/Analytic_Constant_Eloss_Model.hh
 * \author Ben R. Ryan
 * \date   Feb 21 2020
 * \brief  Analytic_Constant_Eloss_Model class definition.
 * \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved. */
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
  Analytic_Constant_Eloss_Model(const int zaidt_, const double mt_, const int zaidp_,
                                const double mp_, const double eloss_value_in) noexcept
      : Analytic_Eloss_Model(zaidt_, mt_, zaidp_, mp_), eloss_value(eloss_value_in) {}

  //! Calculate the eloss rate in units of shk^-1; T given in keV, ndens in a/cc, v0 in cm/shk,
  //! debye in cm
  [[nodiscard]] double
  calculate_eloss([[maybe_unused]] const double T, [[maybe_unused]] const double ndens,
                  [[maybe_unused]] const double v0,
                  [[maybe_unused]] const double debye_length = 0.0) const override {
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
