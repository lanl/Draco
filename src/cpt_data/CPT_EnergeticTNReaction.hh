//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cpt_data/CPT_EnergeticTNReaction.hh
 * \author Kyle S. Beling
 * \brief  CPT_EnergeticTNReaction class definition.
 * \note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cpt_data_reaction_in_flight
#define rtt_cpt_data_reaction_in_flight

#include "cpt_data/CPT_Base.hh"
#include "cpt_data/CPT_Collision.hh"

namespace rtt_cpt_data {
//================================================================================================//
/*!
 * \class CPT_EnergeticTNReaction
 * \brief CPT_EnergeticTNReaction derived class.
 *
 * Wraps all reaction in flight collisions that occur in the specified reaction network into one 
 * vector object. A user specifies the reaction of interrest and the class will return the
 * appropriate Collision object.
 */
//================================================================================================//

class CPT_EnergeticTNReaction : public CPT_Base {

public:
  //! Constructors
  CPT_EnergeticTNReaction(const std::string &dictionaryPathIn, const std::string &datasetIn,
                          std::string reactionNetworkNameIn, std::vector<double> groupBoundsIn);

  CPT_EnergeticTNReaction(const std::string &dictionaryPathIn, std::string reactionNetworkNameIn,
                          std::vector<double> groupBoundsIn);

  CPT_EnergeticTNReaction(const std::string &dictionaryPathIn, const std::string &datasetIn,
                          const std::vector<std::string> &reactionList,
                          std::vector<double> groupBoundsIn);

  //! Default destructor
  ~CPT_EnergeticTNReaction() = default;

  //! Disable default constructor
  CPT_EnergeticTNReaction() = delete;

  //! Disable copy constructor
  CPT_EnergeticTNReaction(const CPT_EnergeticTNReaction &rhs) = delete;

  //! Disable move constructor
  CPT_EnergeticTNReaction(CPT_EnergeticTNReaction &&rhs) noexcept = delete;

  //! Disable assignment
  CPT_EnergeticTNReaction &operator=(const CPT_EnergeticTNReaction &rhs) = delete;

  //! Disable move assignment
  CPT_EnergeticTNReaction &operator=(CPT_EnergeticTNReaction &&rhs) noexcept = delete;

  //! Returns reaction network name
  [[nodiscard]] std::string reaction_network_name() const { return this->reactionNetworkName_; }

  //! Returns list of RIF reactions
  [[nodiscard]] std::vector<std::string> list_of_reactions() const {
    return this->listOfReactions_;
  }

private:
  //! Reaction network name
  std::string reactionNetworkName_{};

  //! Coarse energy groups
  std::vector<double> coarseGroupBounds_{};

  //! List of reactions
  std::vector<std::string> listOfReactions_{};

  //! Reaction data classes
  std::map<std::string, rtt_cpt_data::CPT_Collision> reactionClasses_{};

  // Only implemented if JSON is found
#ifdef JSON_FOUND
  void load_json();
  void load_json(const std::vector<std::string> &reactionList);
#endif

public:
  [[nodiscard]] rtt_cpt_data::CPT_Collision get_reaction(const std::string &reactionString) const {
    return this->reactionClasses_.at(reactionString);
  }
};

} // namespace rtt_cpt_data

#endif // rtt_cpt_data_reaction_in_flight

//------------------------------------------------------------------------------------------------//
// End cpt_data/CPT_EnergeticTNReaction.hh
//------------------------------------------------------------------------------------------------//
