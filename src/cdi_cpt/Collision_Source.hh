//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_cpt/Collision_Source.hh
 * \author Ben R. Ryan
 * \date   2024 Oct 9
 * \brief  Collision_Source class definition.
 * \note   Copyright (C) 2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef cdi_cpt_Collision_Source_hh
#define cdi_cpt_Collision_Source_hh

#include "dsxx/Assert.hh"
#include "dsxx/DracoMath.hh"
#include "shacl/physics.hpp"
#include "cdi/CPCommon.hh"
#include "units/PhysicalConstexprs.hh"
#include <cmath>
#include <fstream>
#include <string>

namespace rtt_cdi_cpt {

class CollisionSource {
public:
  // Constructor
  CollisionSource(const rtt_cdi::CParticle target_in, const rtt_cdi::CParticle projectile_in,
                  const std::string &ace_library,
                  const std::vector<std::string> &ndatk_path_user = {},
                  const std::string &exsdir_path = "exsdir");
  // Constructor with MPI communicator
  CollisionSource(const rtt_cdi::CParticle target_in, const rtt_cdi::CParticle projectile_in,
                  const std::string &ace_library, MPI_Comm comm_in,
                  const std::vector<std::string> &ndatk_path_user = {},
                  const std::string &exsdir_path = "exsdir");

  // Get cross section in cm^2 for this reaction at specified temperature
  [[nodiscard]] double get_cross_section(const double incident_energy,
                                         const double temperature) const;

  // Get cross section in cm^2 for this reaction at specified temperature with interpolation cache
  [[nodiscard]] double
  get_cross_section(const double incident_energy, const double temperature,
                    const shacl::physics::InterpolationCache &interpCache) const;

  // Sample energy of target after collision
  template <typename RNG>
  [[nodiscard]] double sample_target_energy_delta(const double incident_energy,
                                                  const double temperature, const RNG rng) const {
    const int temp_idx = get_data_idx(temperature);
    const auto &proj_targ = projTargetDatas[temp_idx];
    const auto &elasticReaction = proj_targ.reaction(avalanche::ReactionID(avalanche::MT::Elastic));

    const auto energyAndCosine = shacl::physics::samplePrimaryReactionAngleEnergy(
        elasticReaction, incident_energy, awr, rng);

    // Leave it to consumer to correctly interpret negative target energy here
    return incident_energy - energyAndCosine.energy;
  }

  // Get target particle
  [[nodiscard]] rtt_cdi::CParticle get_target() const { return target; }

  // Get projectile particle
  [[nodiscard]] rtt_cdi::CParticle get_projectile() const { return projectile; }

  // Get atomic weight ratio of target
  [[nodiscard]] double get_awr() const { return awr; }

  // Get vector of temperatures of underlying data tables in keV
  [[nodiscard]] std::vector<double> get_temperatures() const { return temperatures; }

  // Get the vector of incident energies at a given temperature
  [[nodiscard]] const auto &get_incident_energies(const double temperature) const {
    return projTargetDatas[get_data_idx(temperature)].incidentEnergies();
  }

private:
  // Temperatures of data tables in keV
  std::vector<double> temperatures;

  // Projectile-target data for each data table
  std::vector<avalanche::ProjectileTarget> projTargetDatas;

  // Atomic weight ratio of target
  double awr{0.0};

  // Target particle
  const rtt_cdi::CParticle target;

  // Projectile particle
  const rtt_cdi::CParticle projectile;

  // Get index into data array for given temperature
  [[nodiscard]] int get_data_idx(const double temperature) const {
    Check(temperature >= 0.);
    size_t idx = 0;
    for (size_t n = 1; n < temperatures.size(); n++) {
      if (temperatures[n] < temperature) {
        idx++;
      } else {
        break;
      }
    }
    Check(idx < projTargetDatas.size());
    Check(idx < INT32_MAX); // is the cast safe?
    return static_cast<int>(idx);
  }
};

} // namespace rtt_cdi_cpt

#endif // cdi_cpt_Collision_Source

//------------------------------------------------------------------------------------------------//
// End cdi_cpt_Collision_Source.hh
//------------------------------------------------------------------------------------------------//
