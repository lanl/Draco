//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cpt_data/CPT_Collision.cc
 * \author Kyle S. Beling
 * \brief  Collision member definitions.
 * \note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cpt_data/CPT_Collision.hh"
#include "cpt_data/config.h" // Definition of JSON_FOUND // NOLINT(misc-include-cleaner)
#include "dsxx/Assert.hh"
#include <string>
#include <vector>

// Only implemented if JSON is found
#ifdef JSON_FOUND
#include "cpt_data/CPT_MultigroupCollapse.hh"
#include <limits>
#include <numeric>
#include <utility>
#endif

namespace rtt_cpt_data {

// Only implemented if JSON is found
#ifdef JSON_FOUND

//================================================================================================//
// CONSTRUCTOR
//================================================================================================//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for JSON reader.
 *
 * \param[in] collisionName string of the collision name
 * \param[in] collisionDict nlohmann::json dictionary containing the multigroup data for the
 *                          collision
 * \param[in] mtNumberIn MT number of the collision
 * \param[in] isHeavyCollisionIn true if the incident particle is heavier than the target
 * \param[in] coarseBounds user defined bounds that the data will be collapsed too
 */
CPT_Collision::CPT_Collision(std::string collisionName,
                             const nlohmann::json &collisionDict, // NOLINT(misc-include-cleaner)
                             const int mtNumberIn, const bool isHeavyCollisionIn,
                             const std::vector<double> &coarseBounds)
    : CPT_MultigroupCollapse{coarseBounds, collisionDict["Energy Grid"].get<std::vector<double>>()},
      collisionName_{std::move(collisionName)} {
  load_collision(collisionDict, mtNumberIn, isHeavyCollisionIn);
}
//------------------------------------------------------------------------------------------------//

//================================================================================================//
// MEMBER FUNCTIONS
//================================================================================================//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Load the collision data into memory
 *
 * \param[in] collisionIn nlohmann::json dictionary containing the multigroup data for the collision
 * \param[in] mtNumberIn MT number of the collision
 * \param[in] isHeavyCollision true if the incident particle is heavier than the target
 */
void CPT_Collision::load_collision(const nlohmann::json &collisionIn, const int mtNumberIn,
                                   const bool isHeavyCollisionIn) {
  // Build MT number string
  std::string const mtNumber = "MT = " + std::to_string(mtNumberIn);

  // Determine if collision is light or heavy
  std::string heavyCollision;
  if (isHeavyCollisionIn) {
    heavyCollision = "Heavy";
  } else {
    heavyCollision = "Light";
  }

  // Check if collision exists in data
  Insist(collisionIn.contains(mtNumber), "Key: \"MT number\" not found in collision file.");
  Require(collisionIn[mtNumber].contains(heavyCollision));

  // Get weights
  Insist(collisionIn.contains("Weights"), "Key: \"Weights\" not found in collision file.");
  auto const &weights = collisionIn["Weights"];

  // Get specific collision data
  nlohmann::json collision = collisionIn[mtNumber][heavyCollision];

  // Check collision data
  check_collision_data(collision);

  // Scattered and recoil particle ZAID numbers
  this->scatteredParticle_ = collision["Scattered Particle"];
  this->recoilingParticle_ = collision["Recoil Particle"];

  // Append list of products
  this->productList_.emplace_back(collision["Scattered Particle"]);
  this->productList_.emplace_back(collision["Recoil Particle"]);

  // Total collision cross section
  this->totalCrossSectionVector_ =
      this->collapse_group_average_values(collision["Sigma Total"], weights);

  // Transfer Matrix (isotropic scattering only)
  auto transferMatrix2d = this->tensor_to_matrix(collision["Transfer Matrix"], 0);
  this->scatteringMatrix_ = this->collapse_transfer_matrices(transferMatrix2d, weights);

  // Recoil Matrix (isotropic scattering only)
  auto recoilMatrix2d = this->tensor_to_matrix(collision["Recoil Matrix"], 0);
  this->recoilMatrix_ = this->collapse_transfer_matrices(recoilMatrix2d, weights);
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Check collision data for all necessary Keys
 */
void CPT_Collision::check_collision_data(const nlohmann::json &collision) {
  // Check if scattered and recoil particle exists
  Insist(collision.contains("Scattered Particle"),
         "Key: \"Scattered Particle\" missing from collision data.");
  Insist(collision.contains("Recoil Particle"),
         "Key: \"Recoil Particle\" missing from collision data.");

  // Check if necessary data exists
  Insist(collision.contains("Sigma Total"), "Key: \"Sigma Total\" missing from collision data.");
  Insist(collision.contains("Transfer Matrix"),
         "Key: \"Transfer Matrix\" missing from collision data.");
  Insist(collision.contains("Recoil Matrix"),
         "Key: \"Recoil Matrix\" missing from collision data.");
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return total cross section for an incident energy
 *
 * \param[in] incidentEnergy of charged particle [MeV]
 *
 * \return Total cross section [barns]
 */
double CPT_Collision::total_cross_section(const double incidentEnergy) const {
  return this->totalCrossSectionVector_[this->get_index(incidentEnergy)];
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return normalized probability distribution function for a product and incident energy
 *
 * \param[in] productZAID ZAID of the outgoing particle
 * \param[in] incidentEnergy of charged particle [MeV]
 *
 * \return Normalized PDF of collision product energy
 */
std::vector<double> CPT_Collision::pdf(const int productZAID, const double incidentEnergy) const {
  const int index = this->get_index(incidentEnergy);
  std::vector<double> distribution(this->numberEnergyGroups_, 0.0);

  if (productZAID == this->scattered_particle()) {
    for (int k = 0; k < this->numberEnergyGroups_; ++k) {
      distribution[k] = this->scatteringMatrix_[k][index];
    }
  } else if (productZAID == this->recoiling_particle()) {
    for (int k = 0; k < this->numberEnergyGroups_; ++k) {
      distribution[k] = this->recoilMatrix_[k][index];
    }
  } else {
    Insist(0, "productZAID not found in collision");
  }

  const double sum = std::accumulate(distribution.begin(), distribution.end(), 0.0);
  if (sum > std::numeric_limits<double>::epsilon()) {
    for (double &value : distribution)
      value /= sum;
  }
  return distribution;
}
//------------------------------------------------------------------------------------------------//

// Stubbed out member functions
#else

//! Returns collision name
[[nodiscard]] std::string CPT_Collision::name() const {
  Insist(0, "name() only available when the nlohmann::json library is found");
}

//! Returns scattered particle ZAID
[[nodiscard]] int CPT_Collision::scattered_particle() const {
  Insist(0, "scattered_particle() only available when the nlohmann::json library is found");
}

//! Returns recoil particle ZAID
[[nodiscard]] int CPT_Collision::recoiling_particle() const {
  Insist(0, "recoiling_particle() only available when the nlohmann::json library is found");
}

//! Returns ZAID list of products
[[nodiscard]] std::vector<int> CPT_Collision::product_list() const {
  Insist(0, "product_list() only available when the nlohmann::json library is found");
}

//! Energy-dependent total cross section
[[nodiscard]] double CPT_Collision::total_cross_section(double /*incidentEnergy*/) const {
  Insist(0, "total_cross_section() only available when the nlohmann::json library is found");
}

//! Returns PDF for a specific product and incident energy
[[nodiscard]] std::vector<double> CPT_Collision::pdf(const int /*productZAID*/,
                                                     const double /*incidentEnergy*/) const {
  Insist(0, "pdf() only available when the nlohmann::json library is found");
}

#endif

} // namespace rtt_cpt_data

//------------------------------------------------------------------------------------------------//
// End cpt_data/CPT_Collision.cc
//------------------------------------------------------------------------------------------------//
