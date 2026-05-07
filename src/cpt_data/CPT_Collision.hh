//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cpt_data/CPT_Collision.hh
 * \author Kyle S. Beling
 * \brief  Collision class definition.
 * \note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cpt_data_Collision
#define rtt_cpt_data_Collision

#include "cpt_data/CPT_MultigroupCollapse.hh"
#include "cpt_data/config.h" // Definition of JSON_FOUND

#include "dsxx/Assert.hh"

namespace rtt_cpt_data {
//================================================================================================//
/*!
 * \class CPT_Collision
 * \brief Collision class.
 *
 * The cross section data is stored in this class in the totalCrossSectionVector_, scatteringMatrix_,
 * and recoilMatrix_. Energies are in units of MeV and temperatures are in units of keV. Velocity-
 * averaged cross sections are in units of cm^3 sh^-1. Probability density functions sum to unity. 
 * Note: this data is only multigroup, the "equivalent" continuous energy and angle is availible 
 * through ACE; however, at the time of writing this class significant formatting and particle 
 * kinematic errors exist in the ACE data.
 */
//================================================================================================//

// Only implemented if JSON is found
#ifdef JSON_FOUND

class CPT_Collision : public CPT_MultigroupCollapse {

public:
  //! Constructor
  CPT_Collision(const std::string &collisionName, const nlohmann::json &collisionDict,
                const int mtNumberIn, const bool isHeavyCollisionIn,
                const std::vector<double> &coarseBounds);

  //! Disable default constructor
  CPT_Collision() = delete;

  //! Returns collision name
  [[nodiscard]] std::string name() const { return collisionName_; }

  //! Returns scattered particle ZAID
  [[nodiscard]] int scattered_particle() const { return this->scatteredParticle_; }

  //! Returns recoil particle ZAID
  [[nodiscard]] int recoiling_particle() const { return this->recoilingParticle_; }

  //! Returns ZAID list of products
  [[nodiscard]] std::vector<int> product_list() const { return this->productList_; }

  //! Energy-dependent total cross section
  [[nodiscard]] double total_cross_section(double incidentEnergy) const;

  //! Returns PDF for a specific product and incident energy
  [[nodiscard]] std::vector<double> pdf(const int productZAID, const double incidentEnergy) const;

private:
  //! Collision name
  std::string collisionName_{};

  //! Scattered particle
  int scatteredParticle_{};

  //! Recoiling particle
  int recoilingParticle_{};

  //! List of produced particles ZAIDs
  std::vector<int> productList_{};

  //! Total cross section
  std::vector<double> totalCrossSectionVector_{};

  //! Scattering matrix (only isotropic scattering)
  std::vector<std::vector<double>> scatteringMatrix_{};

  //! Recoil matrix (only isotropic scattering)
  std::vector<std::vector<double>> recoilMatrix_{};

  //! Load in collision data
  void load_collision(const nlohmann::json &collisionIn, const int mtNumberIn,
                      const bool isHeavyCollisionIn);

  //! Check collision data
  void check_collision_data(const nlohmann::json &collision);
};

#else

class CPT_Collision : public CPT_MultigroupCollapse {

public:
  //! Returns collision name
  [[nodiscard]] std::string name() const;

  //! Returns scattered particle ZAID
  [[nodiscard]] int scattered_particle() const;

  //! Returns recoil particle ZAID
  [[nodiscard]] int recoiling_particle() const;

  //! Returns ZAID list of products
  [[nodiscard]] std::vector<int> product_list() const;

  //! Energy-dependent total cross section
  [[nodiscard]] double total_cross_section(double /*incidentEnergy*/) const;

  //! Returns PDF for a specific product and incident energy
  [[nodiscard]] std::vector<double> pdf(const int /*productZAID*/,
                                        const double /*incidentEnergy*/) const;
};

#endif

} // namespace rtt_cpt_data

#endif // rtt_cpt_data_Collision

//------------------------------------------------------------------------------------------------//
// End cpt_data/CPT_Collision.hh
//------------------------------------------------------------------------------------------------//
