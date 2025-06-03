//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/Analytic_LiPet_Eloss_Model.hh
 * \author Edward Norris and Kyle S. Beling
 * \date   May 21, 2025
 * \brief  Analytic_LiPet_Eloss_Model class definition.
 * \note   Copyright (C) 2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_CPEloss_Analytic_LiPet_Eloss_Model_hh
#define rtt_cdi_CPEloss_Analytic_LiPet_Eloss_Model_hh

#include <cmath>
#include <vector>

namespace rtt_cdi_cpeloss {

class Analytic_LiPet_Eloss_Model {
  int incident_zaid_ = 0;
  double incident_mass_ = 0.0; // kg
  int incident_z_ = 0;

public:
  //! Constructor
  Analytic_LiPet_Eloss_Model(const int zaid, const double mass)
      : incident_zaid_(zaid), incident_mass_(mass) {
    if (zaid == -1) {
      incident_z_ = 1;
    } else {
      incident_z_ = (zaid - zaid % 1000) / 1000;
    }
  }

  //! ZAID of incident particle
  [[nodiscard]] int incident_zaid() const { return this->incident_zaid_; }

  //! Mass of incident particle
  [[nodiscard]] double incident_mass() const { return this->incident_mass_; }

  //! z-number of incident particle
  [[nodiscard]] int z_number() const { return this->incident_z_; }

  //! Calculate the total eloss rate for some background composition
  [[nodiscard]] double getEloss(const double incident_energy, const std::vector<int> &zaid_list,
                                const std::vector<double> &number_densities,
                                const std::vector<double> &temperatures,
                                const std::vector<double> &masses) const;

  //! Calculate the modified parameterization total eloss rate for some background composition
  [[nodiscard]] double getElossMod(const double incident_energy, const std::vector<int> &zaid_list,
                                   const std::vector<double> &number_densities,
                                   const std::vector<double> &temperatures,
                                   const std::vector<double> &masses) const;

private:
  //! Calculate the total debye length for some background composition
  [[nodiscard]] double calculate_debye_length(const std::vector<int> &zaid_list,
                                              const std::vector<double> &number_densities,
                                              const std::vector<double> &temperatures) const;

  //! Calculate the stopping power for some target
  [[nodiscard]] double calculate_eloss(const double incident_energy, const double number_density,
                                       const double temperature, const double target_mass,
                                       const double debye_length) const;

  //! Calculate the modified Li-Petrasso stopping power
  [[nodiscard]] double calculate_eloss_modified(const double incident_energy,
                                                const double number_density,
                                                const double temperature, const double target_mass,
                                                const double debye_length) const;
};

} // namespace rtt_cdi_cpeloss

#endif // rtt_cdi_CPEloss_Analytic_LiPet_Eloss_Model_hh

//------------------------------------------------------------------------------------------------//
// End cdi_CPEloss/Analytic_Spitzer_Eloss_Model.hh
//------------------------------------------------------------------------------------------------//
