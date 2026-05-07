//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cpt_data/CPT_ElasticScattering.hh
 * \author Kyle S. Beling
 * \brief  CPT_ElasticScattering class definition.
 * \note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cpt_data_elastic_scattering
#define rtt_cpt_data_elastic_scattering

#include "cpt_data/CPT_Base.hh"
#include "cpt_data/CPT_Collision.hh"

namespace rtt_cpt_data {
//================================================================================================//
/*!
 * \class CPT_ElasticScattering
 * \brief CPT_ElasticScattering derived class.
 *
 * Wraps all elastic scattering collisions that occur in the specified reaction network into one 
 * vector object. A user specifies the collision of interrest and the class will return the
 * appropriate Collision object.
 */
//================================================================================================//

class CPT_ElasticScattering : public CPT_Base {

public:
  //! Constructors
  CPT_ElasticScattering(const std::string &dictionaryPathIn, const std::string &datasetIn,
                        std::string reactionNetworkNameIn, std::vector<double> groupBoundsIn);

  CPT_ElasticScattering(const std::string &dictionaryPathIn, std::string reactionNetworkNameIn,
                        std::vector<double> groupBoundsIn);

  CPT_ElasticScattering(const std::string &dictionaryPathIn, const std::string &datasetIn,
                        const std::vector<std::string> &listOfCollisions,
                        std::vector<double> groupBoundsIn);

  //! Default destructor
  ~CPT_ElasticScattering() = default;

  //! Disable default constructor
  CPT_ElasticScattering() = delete;

  //! Disable copy constructor
  CPT_ElasticScattering(const CPT_ElasticScattering &rhs) = delete;

  //! Disable move constructor
  CPT_ElasticScattering(CPT_ElasticScattering &&rhs) noexcept = delete;

  //! Disable assignment
  CPT_ElasticScattering &operator=(const CPT_ElasticScattering &rhs) = delete;

  //! Disable move assignment
  CPT_ElasticScattering &operator=(CPT_ElasticScattering &&rhs) noexcept = delete;

  [[nodiscard]] std::vector<std::string> list_of_collisions() const {
    return this->listOfCollisions_;
  }

  [[nodiscard]] std::string reaction_network_name() const { return this->reactionNetworkName_; }

private:
  //! Reaction network name
  std::string reactionNetworkName_{};

  //! Coarse energy groups
  std::vector<double> coarseGroupBounds_{};

  //! List of reactions
  std::vector<std::string> listOfCollisions_{};

  //! Reaction data classes
  std::map<std::string, rtt_cpt_data::CPT_Collision> collisionClasses_{};

  // Only implemented if JSON is found
#ifdef JSON_FOUND
  void load_json();
  void load_json(const std::vector<std::string> &listOfCollisions);
#endif

public:
  [[nodiscard]] rtt_cpt_data::CPT_Collision get_collision(const std::string &collisionName) const {
    return this->collisionClasses_.at(collisionName);
  }
};

} // namespace rtt_cpt_data

#endif // rtt_cpt_data_reaction_in_flight

//------------------------------------------------------------------------------------------------//
// End cpt_data/CPT_ElasticScattering.hh
//------------------------------------------------------------------------------------------------//
