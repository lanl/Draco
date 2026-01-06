//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_cpt/Collision_Source.hh
 * \author Ben R. Ryan
 * \date   2024 Oct 9
 * \brief  Collision_Source class definition.
 * \note   Copyright (C) 2025-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef cdi_cpt_Collision_Source_hh
#define cdi_cpt_Collision_Source_hh

// Defines SHACL_PHYSICS_FOUND:
#include "cdi_cpt/config.h" // NOLINT(misc-include-cleaner)
#include "dsxx/Assert.hh"
#include "dsxx/DracoMath.hh"
#ifdef SHACL_PHYSICS_FOUND
#include "shacl/physics.hpp"
#endif
#include "c4/c4_mpi.h"
#include "cdi/CPCommon.hh"
#include "units/PhysicalConstexprs.hh"
#include <cmath>
#include <fstream>
#include <string>

namespace rtt_cdi_cpt {

class CollisionSource {
public:
  // Constructor
  CollisionSource(const int32_t target_zaid_in, const int32_t projectile_zaid_in,
                  const std::string &ace_library,
                  const std::vector<std::string> &ndatk_path_user = {},
                  const std::string &exsdir_path = "exsdir");

// Need ifdef for scalar builds.
#ifdef C4_MPI
  // Constructor with MPI communicator
  CollisionSource(const int32_t target_zaid_in, const int32_t projectile_zaid_in,
                  const std::string &ace_library, MPI_Comm comm_in,
                  const std::vector<std::string> &ndatk_path_user = {},
                  const std::string &exsdir_path = "exsdir");
#endif

  // Get cross section in cm^2 for this reaction at specified temperature
  [[nodiscard]] double get_cross_section(const double incident_energy,
                                         const double temperature) const;

#ifdef SHACL_PHYSICS_FOUND
  // Get cross section in cm^2 for this reaction at specified temperature with interpolation cache
  [[nodiscard]] double
  get_cross_section(const double incident_energy, const double temperature,
                    const shacl::physics::InterpolationCache &interpCache) const;

  // Sample energy of target after collision
  template <typename RNG>
  [[nodiscard]] double sample_target_energy_delta(const double incident_energy,
                                                  const double temperature, RNG &rng) const {
    const int temp_idx = get_data_idx(temperature);
    const auto &proj_targ = projTargetDatas[temp_idx];
    const auto &elasticReaction = proj_targ.reaction(avalanche::ReactionID(avalanche::MT::Elastic));

    // This lambda wraps our ran() function in an operator-like call that shacl-physics expects.
    auto ran_operator = [&rng]() { return rng.ran(); };

    const auto energyAndCosine = shacl::physics::samplePrimaryReactionAngleEnergy(
        elasticReaction, incident_energy, awr, ran_operator);

    // Leave it to consumer to correctly interpret negative target energy here
    return incident_energy - energyAndCosine.energy;
  }

  // Get the vector of incident energies at a given temperature
  [[nodiscard]] const auto &get_incident_energies(const double temperature) const {
    return projTargetDatas[get_data_idx(temperature)].incidentEnergies();
  }
#else // !(SHACL_PHYSICS_FOUND)
  // Provide stubbed interface functions that throw if called without shacl::physics support.
  template <typename RNG>
  [[nodiscard]] double sample_target_energy_delta(const double /*incident_energy*/,
                                                  const double /*temperature*/,
                                                  RNG & /*rng*/) const {
    Insist(false,
           "CollisionSource::sample_target_energy_delta() only defined when SHACL_PHYSICS_FOUND!");
    return 0.0;
  }
#endif

  // Get target particle ZAID
  [[nodiscard]] int32_t get_target() const { return target_zaid; }

  // Get projectile particle ZAID
  [[nodiscard]] int32_t get_projectile() const { return projectile_zaid; }

  // Get atomic weight ratio of target
  [[nodiscard]] double get_awr() const { return awr; }

  // Get vector of temperatures of underlying data tables in keV
  [[nodiscard]] std::vector<double> get_temperatures() const { return temperatures; }

private:
  // Temperatures of data tables in keV
  std::vector<double> temperatures;

#ifdef SHACL_PHYSICS_FOUND
  // Projectile-target data for each data table
  std::vector<avalanche::ProjectileTarget> projTargetDatas;
#endif

  // Atomic weight ratio of target
  double awr{0.0};

  // Target particle
  const int32_t target_zaid{0};

  // Projectile particle
  const int32_t projectile_zaid{0};

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
#ifdef SHACL_PHYSICS_FOUND
    Check(idx < projTargetDatas.size());
#endif
    Check(idx < INT32_MAX); // is the cast safe?
    return static_cast<int>(idx);
  }
};

} // namespace rtt_cdi_cpt

#endif // cdi_cpt_Collision_Source

//------------------------------------------------------------------------------------------------//
// End cdi_cpt_Collision_Source.hh
//------------------------------------------------------------------------------------------------//
