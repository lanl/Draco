//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cpt_data/CPT_EnergeticTNReaction.cc
 * \author Kyle S. Beling
 * \date   Nov 19 2025
 * \brief  CPT_EnergeticTNReaction member definitions.
 * \note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cpt_data/CPT_EnergeticTNReaction.hh"
#include "cpt_data/CPT_Base.hh"

#include "dsxx/Assert.hh"
#include <string>
#include <utility>
#include <vector>

// cSpell:ignore NOLINTNEXTLINE, ccolor

namespace rtt_cpt_data {

//================================================================================================//
// CONSTRUCTORS
//================================================================================================//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Normal constructor for RectionInflight class.
 *
 * \param[in] dictionaryPathIn path to dictionary file
 * \param[in] datasetIn name of requested dataset
 * \param[in] reactionNetworkNameIn name of requested reaction network
 * \param[in] groupBoundsIn energy boundaries of multigroup bins (MeV)
 */
CPT_EnergeticTNReaction::CPT_EnergeticTNReaction(const std::string &dictionaryPathIn,
                                                 const std::string &datasetIn,
                                                 std::string reactionNetworkNameIn,
                                                 std::vector<double> groupBoundsIn)
    : CPT_Base(dictionaryPathIn, datasetIn, "reaction_in_flight"),
      reactionNetworkName_(std::move(reactionNetworkNameIn)),
      coarseGroupBounds_(std::move(groupBoundsIn)) {
#ifdef JSON_FOUND
  load_json();
#endif
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for RectionInflight class with no dataset specified.
 *
 * \param[in] dictionaryPathIn path to dictionary file
 * \param[in] reactionNetworkNameIn name of requested reaction network
 * \param[in] groupBoundsIn energy boundaries of multigroup bins (MeV)
 */
CPT_EnergeticTNReaction::CPT_EnergeticTNReaction(const std::string &dictionaryPathIn,
                                                 std::string reactionNetworkNameIn,
                                                 std::vector<double> groupBoundsIn)
    : CPT_Base(dictionaryPathIn, "reaction_in_flight"),
      reactionNetworkName_(std::move(reactionNetworkNameIn)),
      coarseGroupBounds_(std::move(groupBoundsIn)) {
#ifdef JSON_FOUND
  load_json();
#endif
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for CPT_EnergeticTNReaction class with a list of RIF reactions specified.
 *
 * \param[in] dictionaryPathIn path to dictionary file
 * \param[in] datasetIn name of requested dataset
 * \param[in] reactionList name of requested reaction network
 * \param[in] groupBoundsIn energy boundaries of multigroup bins (MeV)
 */
CPT_EnergeticTNReaction::CPT_EnergeticTNReaction(const std::string &dictionaryPathIn,
                                                 const std::string &datasetIn,
                                                 const std::vector<std::string> &reactionList,
                                                 std::vector<double> groupBoundsIn)
    : CPT_Base(dictionaryPathIn, datasetIn, "reaction_in_flight"), reactionNetworkName_("Custom"),
      coarseGroupBounds_(std::move(groupBoundsIn)), listOfReactions_(reactionList) {

  // Ensure reactionList is not empty
  Require(reactionList.size() > 0);

#ifdef JSON_FOUND
  load_json(reactionList);
#endif
}
//------------------------------------------------------------------------------------------------//

//================================================================================================//
// MEMBER FUNCTIONS
//================================================================================================//

// Protect actual JSON calls with JSON_FOUND macro:
#ifdef JSON_FOUND

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Load and collapse JSON TN reaction data into memory for each reaction in the specified 
 *        reaction network.
 */
void CPT_EnergeticTNReaction::load_json() {
  // Open and read reaction network info into variable: reactionNetworkInfo
  std::ifstream reactionNetworkFile(this->reactionNetworkFileName_);
  nlohmann::json reactionNetworkInfo;
  reactionNetworkFile >> reactionNetworkInfo;
  reactionNetworkFile.close();

  // Open and read TN reaction info into variable: reactionNetworkInfo
  std::ifstream reactionsFile(this->tnReactionFileName_);
  nlohmann::json tnReactionInfo;
  reactionsFile >> tnReactionInfo;
  reactionsFile.close();

  // Get TN reactions of interest
  auto tnReactions = reactionNetworkInfo[this->reactionNetworkName_];

  // Loop through reactions in the reaction network
  for (const auto &reactionName : tnReactions) {
    // Add reaction to list_of_reactions
    this->listOfReactions_.emplace_back(static_cast<std::string>(reactionName));

    // Check if reaction is in tnReactionInfo dictionary
    Require(tnReactionInfo.contains(reactionName));

    // Check if the necessary keys exist in the tnReactionInfo dictinary exist
    Require(tnReactionInfo[reactionName].contains("File"));
    Require(tnReactionInfo[reactionName].contains("MT"));
    Require(tnReactionInfo[reactionName].contains("Heavy"));

    // Open reaction data file
    std::string const fileName = tnReactionInfo[reactionName]["File"];
    int const reactionNumber = tnReactionInfo[reactionName]["MT"];
    bool const isHeavyReaction = tnReactionInfo[reactionName]["Heavy"];

    std::string const file_path =
        this->dictionaryPath_ + "/" +
        static_cast<std::string>(this->dictionary_["Data Sets"][this->dataset_]["Location"]) + "/" +
        fileName;

    // Ensure file exists
    Require(rtt_dsxx::fileExists(file_path));

    // Open and read individual reaction file
    std::ifstream reactionFile(file_path);
    nlohmann::json reactionDictionary;
    reactionFile >> reactionDictionary;
    reactionFile.close();

    // Build reaction data
    rtt_cpt_data::CPT_Collision rifReaction(static_cast<std::string>(reactionName),
                                            reactionDictionary, reactionNumber, isHeavyReaction,
                                            this->coarseGroupBounds_);

    // Add reaction to map
    reactionClasses_.insert(
        std::make_pair(static_cast<std::string>(reactionName), std::move(rifReaction)));
  }
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Load and collapse JSON TN reaction data into memory for each reaction in the user 
 *        provided variable reactionList
 *
 * \param[in] reactionList contains a list of user specified TN reactions to parsed
 */
void CPT_EnergeticTNReaction::load_json(const std::vector<std::string> &reactionList) {
  // Open and read TN reaction info into variable: reactionNetworkInfo
  std::ifstream reactionsFile(this->tnReactionFileName_);
  nlohmann::json tnReactionInfo;
  reactionsFile >> tnReactionInfo;
  reactionsFile.close();

  // Loop through reactions in the reaction network
  for (const auto &reactionName : reactionList) {
    // Open reaction data file
    std::string const fileName = tnReactionInfo[reactionName]["File"];
    int const reactionNumber = tnReactionInfo[reactionName]["MT"];
    bool const isHeavyReaction = tnReactionInfo[reactionName]["Heavy"];

    // Check if reaction is in tnReactionInfo dictionary
    Require(tnReactionInfo.contains(reactionName));

    // Check if the necessary keys exist in the tnReactionInfo dictinary exist
    Require(tnReactionInfo[reactionName].contains("File"));
    Require(tnReactionInfo[reactionName].contains("MT"));
    Require(tnReactionInfo[reactionName].contains("Heavy"));

    std::string const file_path =
        this->dictionaryPath_ + "/" +
        static_cast<std::string>(this->dictionary_["Data Sets"][this->dataset_]["Location"]) + "/" +
        fileName;

    // Ensure file exists
    Require(rtt_dsxx::fileExists(file_path));

    // Open and read individual reaction file
    std::ifstream reactionFile(file_path);
    nlohmann::json reactionDictionary;
    reactionFile >> reactionDictionary;
    reactionFile.close();

    // Build reaction data
    rtt_cpt_data::CPT_Collision rifReaction(static_cast<std::string>(reactionName),
                                            reactionDictionary, reactionNumber, isHeavyReaction,
                                            this->coarseGroupBounds_);

    // Add reaction to map
    reactionClasses_.insert(
        std::make_pair(static_cast<std::string>(reactionName), std::move(rifReaction)));
  }
}
//------------------------------------------------------------------------------------------------//

#else

#endif

} // namespace rtt_cpt_data

//------------------------------------------------------------------------------------------------//
// End cpt_data/CPT_EnergeticTNReaction.cc
//------------------------------------------------------------------------------------------------//
