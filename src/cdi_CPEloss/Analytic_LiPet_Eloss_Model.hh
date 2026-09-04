//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/Analytic_LiPet_Eloss_Model.hh
 * \author Edward Norris and Kyle S. Beling
 * \brief  Analytic_LiPet_Eloss_Model class definition.
 * \note   Copyright (C) 2025-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_CPEloss_Analytic_LiPet_Eloss_Model_hh
#define rtt_cdi_CPEloss_Analytic_LiPet_Eloss_Model_hh

#include "cdi_CPEloss/Analytic_Eloss_Model.hh"
#include <cmath>
#include <vector>

namespace rtt_cdi_cpeloss {

class Analytic_LiPet_Eloss_Model : public Analytic_Eloss_Model {

public:
  //! Constructor
  Analytic_LiPet_Eloss_Model(const int tzaid, const double tmass, const int pzaid,
                             const double pmass, const bool edward_rel_vel = false) noexcept
      : Analytic_Eloss_Model(tzaid, tmass, pzaid, pmass), edward_relative_velocity(edward_rel_vel) {
  }

  //! This doesn't overload an Analytic_Eloss_Model function - it wraps an alternate version
  //! of Li-Petrasso for comparison.
  [[nodiscard]] double getElossMod(const double temperature_in, const double number_density_in,
                                   const double incident_speed_in,
                                   const double debye_length_in) const;

private:
  //! Toggle to use Edward's relative velocity for standard LiPet dedt.
  bool edward_relative_velocity{false};

  //! Calculate the stopping power for some target
  [[nodiscard]] double calculate_eloss_modified(const double temperature_in,
                                                const double number_density_in,
                                                const double incident_speed_in,
                                                const double debye_length_in) const;

  //! Calculate the modified Li-Petrasso stopping power
  [[nodiscard]] double calculate_eloss( // NOLINT(misc-override-with-different-visibility)
      const double temperature_in, const double number_density_in, const double incident_speed_in,
      const double debye_length_in) const override;
};

} // namespace rtt_cdi_cpeloss

#endif // rtt_cdi_CPEloss_Analytic_LiPet_Eloss_Model_hh

//------------------------------------------------------------------------------------------------//
// End cdi_CPEloss/Analytic_LiPet_Eloss_Model.hh
//------------------------------------------------------------------------------------------------//
