//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_cpt/Collision_Source.cc
 * \author Ben R. Ryan
 * \date   2024 Oct 9
 * \brief  Collision_Source member definitions.
 * \note   Copyright (C) 2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi_cpt/Collision_Source.hh"
#include "dsxx/Assert.hh"
#include "dsxx/DracoStrings.hh"
#include "dsxx/Interpolate.hh"
#include "cdi/CPCommon.hh"
#include <cstdint> // uint32_t
#include <utility>
#include <vector>

namespace rtt_cdi_cpt {

avalanche::ParticleID GetParticleID(const rtt_cdi::CParticle &particle) {
  if (particle.get_zaid() == -1) {
    return avalanche::ParticleID::electron();
  } else {
    return avalanche::ParticleID(particle.get_zaid() * 100);
  }
}

CollisionSource::CollisionSource(const rtt_cdi::CParticle target_in,
                                 const rtt_cdi::CParticle projectile_in,
                                 const std::string ace_library,
                                 const std::vector<std::string> ndatk_path_user,
                                 const std::string exsdir_path)
    : target(target_in), projectile(projectile_in) {

  auto atarget = GetParticleID(target_in);
  // TODO(BRR) verify that projectile_in is correct for this data table?
  // auto aprojectile = GetParticleID(projectile_in);

  // If no custom path provided, use ndatk default
  std::vector<std::string> ndatk_path =
      ndatk_path_user.empty() ? ndatk::Exsdir::default_path() : ndatk_path_user;

  avalanche::AceLibrary ace_lib(exsdir_path, ace_library, ndatk_path);
  temperatures = ace_lib.temperatures(std::to_string(target_in.get_zaid()));
  Insist(!temperatures.empty(), "No temperature points found for this target!");

  // Create a vector of elastic reactions at each temperature
  for (size_t n = 0; n < temperatures.size(); n++) {
    // Use each temperature value of datafiles to extract correct table
    auto ace_path = ace_lib.path(std::to_string(target_in.get_zaid()), temperatures[n]);

    // Create projectile target data, using MPI if available to reduce file reads
    avalanche::ace::ProjectileTargetFactory proj_targ_fac;
    proj_targ_fac.secondaryParticleTypes({atarget});
#ifdef SHACL_STORM_MPI
    auto proj_targ = proj_targ_fac.create(ace_path, shacl::storm::mpi::Communicator{});
#else
    auto proj_targ = proj_targ_fac.create(ace_path);
#endif

    // Store atomic weight ratio (assume same value at all temperatures)
    awr = proj_targ.AWR();

    // Store projectile target data in vector
    projTargetDatas.push_back(std::move(proj_targ));

    // Convert this temperature to keV
    temperatures[n] *= 1000.;
  }
}

double CollisionSource::get_cross_section(const double incident_energy,
                                          const double temperature) const {
  const int temp_idx = get_data_idx(temperature);
  const auto &proj_targ = projTargetDatas[temp_idx];

  // Make temporary interpolation cache
  auto interpCache = shacl::physics::makeInterpolationCache(
      proj_targ.incidentEnergies(), incident_energy / 1000., avalanche::InterpolationType::LogLog);

  return get_cross_section(incident_energy, temperature, interpCache);
}

double
CollisionSource::get_cross_section(const double /*incident_energy*/, const double temperature,
                                   const shacl::physics::InterpolationCache &interpCache) const {
  const int temp_idx = get_data_idx(temperature);
  const auto &proj_targ = projTargetDatas[temp_idx];
  const auto &elasticReaction = proj_targ.reaction(avalanche::ReactionID(avalanche::MT::Elastic));

  // Get cross section in barns
  const double cross_section =
      shacl::physics::computeCrossSection(elasticReaction.crossSection(), interpCache);

  // Return cross section in cm^2
  return cross_section * 1.e-24;
}

} // namespace rtt_cdi_cpt

//------------------------------------------------------------------------------------------------//
// End cdi_cpt/Collision_Source.cc
//------------------------------------------------------------------------------------------------//
