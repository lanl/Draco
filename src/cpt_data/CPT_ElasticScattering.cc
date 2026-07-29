//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cpt_data/CPT_ElasticScattering.cc
 * \author Kyle S. Beling
 * \brief  CPT_ElasticScattering member definitions.
 * \note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cpt_data/CPT_ElasticScattering.hh"
#include "cpt_data/CPT_Base.hh"
#include "dsxx/Assert.hh"
#include <string>
#include <utility>
#include <vector>

#ifdef JSON_FOUND
#include "cpt_data/CPT_Collision.hh"
#include <filesystem>
#include <fstream>
#endif

namespace rtt_cpt_data {

//================================================================================================//
// CONSTRUCTORS
//================================================================================================//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for JSON reader.
 *
 * \param[in] dictionaryPathIn path to dictionary file
 * \param[in] datasetIn name of requested temperature dataset
 * \param[in] reactionNetworkNameIn name of requested reaction network
 * \param[in] groupBoundsIn energy boundaries of multigroup bins (MeV)
 */
CPT_ElasticScattering::CPT_ElasticScattering(const std::string &dictionaryPathIn,
                                             const std::string &datasetIn,
                                             std::string reactionNetworkNameIn,
                                             std::vector<double> groupBoundsIn)
    : CPT_Base(dictionaryPathIn, datasetIn, "elastic_scattering"),
      reactionNetworkName_(std::move(reactionNetworkNameIn)),
      coarseGroupBounds_(std::move(groupBoundsIn)) {
#ifdef JSON_FOUND
  load_json();
#endif
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for JSON reader.
 *
 * \param[in] dictionaryPathIn path to dictionary file
 * \param[in] reactionNetworkNameIn name of requested reaction network
 * \param[in] groupBoundsIn energy boundaries of multigroup bins (MeV)
 */
CPT_ElasticScattering::CPT_ElasticScattering(const std::string &dictionaryPathIn,
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
 * \brief Constructor for JSON reader.
 *
 * \param[in] dictionaryPathIn path to dictionary file
 * \param[in] datasetIn name of requested temperature dataset
 * \param[in] listOfCollisions list of collisions to be parsed
 * \param[in] groupBoundsIn energy boundaries of multigroup bins (MeV)
 */
CPT_ElasticScattering::CPT_ElasticScattering(const std::string &dictionaryPathIn,
                                             const std::string &datasetIn,
                                             const std::vector<std::string> &listOfCollisions,
                                             std::vector<double> groupBoundsIn)
    : CPT_Base(dictionaryPathIn, datasetIn, "elastic_scattering"), reactionNetworkName_("Custom"),
      coarseGroupBounds_(std::move(groupBoundsIn)), listOfCollisions_(listOfCollisions) {
  // Ensure listOfCollisions is not empty
  Require(listOfCollisions.size() > 0);

#ifdef JSON_FOUND
  load_json(listOfCollisions);
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
 * \brief Load and collapse JSON elastic scattering data into memory for each reaction in the  
 *        specified reaction network.
 */
void CPT_ElasticScattering::load_json() {
  // Open and read reaction network info into variable: reactionNetworkInfo
  std::ifstream reactionNetworkFile(this->reactionNetworkFileName_);
  nlohmann::json reactionNetworkInfo; // NOLINT(misc-include-cleaner)
  reactionNetworkFile >> reactionNetworkInfo;
  reactionNetworkFile.close();

  // Open and read Elastic scattering info into variable: elasticScatteringInfo
  std::ifstream elasticScatteringFile(this->elasticScatteringFileName_);
  nlohmann::json elasticScatteringInfo;
  elasticScatteringFile >> elasticScatteringInfo;
  elasticScatteringFile.close();

  // Get TN reactions of interest
  auto tnReactions = reactionNetworkInfo[this->reactionNetworkName_];

  // Loop through reactions in the reaction network
  for (const auto &reactionName : tnReactions) {
    // Add collision to list_of_reactions
    std::string const collisionName = static_cast<std::string>(reactionName).substr(0, 3);
    this->listOfCollisions_.emplace_back(collisionName);

    // Check if collision is in elasticScatteringInfo dictionary
    Require(elasticScatteringInfo.contains(collisionName));

    // Check if the necessary keys exist in the elasticScatteringInfo dictinary exist
    Require(elasticScatteringInfo[collisionName].contains("File"));
    Require(elasticScatteringInfo[collisionName].contains("MT"));
    Require(elasticScatteringInfo[collisionName].contains("Heavy"));

    // Open collision data file
    std::string const fileName = elasticScatteringInfo[collisionName]["File"];
    int const collisionNumber = elasticScatteringInfo[collisionName]["MT"];
    bool const isHeavyCollision = elasticScatteringInfo[collisionName]["Heavy"];

    std::string const file_path =
        this->dictionaryPath_ + "/" +
        static_cast<std::string>(this->dictionary_["Data Sets"][this->dataset_]["Location"]) + "/" +
        fileName;

    // Ensure file exists
    Require(std::filesystem::exists(file_path));

    // Open and read individual collision file
    std::ifstream collisionFile(file_path);
    nlohmann::json collisionDictionary;
    collisionFile >> collisionDictionary;
    collisionFile.close();

    // Build reaction data
    rtt_cpt_data::CPT_Collision elasticCollision(collisionName, collisionDictionary,
                                                 collisionNumber, isHeavyCollision,
                                                 this->coarseGroupBounds_);

    // Add reaction to map
    collisionClasses_.insert(std::make_pair(collisionName, std::move(elasticCollision)));
  }
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Load and collapse JSON elastic scattering data into memory for each reaction in the 
 *        user provided variable listOfCollisions
 *
 * \param[in] listOfCollisions contains a list of user specified elastic scattering collisions
 *            to parsed
 */
void CPT_ElasticScattering::load_json(const std::vector<std::string> &listOfCollisions) {
  // Open and read Elastic scattering info into variable: elasticScatteringInfo
  std::ifstream elasticScatteringFile(this->elasticScatteringFileName_);
  nlohmann::json elasticScatteringInfo;
  elasticScatteringFile >> elasticScatteringInfo;
  elasticScatteringFile.close();

  // Loop through reactions in the reaction network
  for (const auto &reactionName : listOfCollisions) {
    // Add collision to list_of_reactions
    std::string const collisionName = static_cast<std::string>(reactionName).substr(0, 3);

    // Check if collision is in elasticScatteringInfo dictionary
    Require(elasticScatteringInfo.contains(collisionName));

    // Check if the necessary keys exist in the elasticScatteringInfo dictionary exist
    Require(elasticScatteringInfo[collisionName].contains("File"));
    Require(elasticScatteringInfo[collisionName].contains("MT"));
    Require(elasticScatteringInfo[collisionName].contains("Heavy"));

    // Open collision data file
    std::string const fileName = elasticScatteringInfo[collisionName]["File"];
    int const collisionNumber = elasticScatteringInfo[collisionName]["MT"];
    bool const isHeavyCollision = elasticScatteringInfo[collisionName]["Heavy"];

    std::string const file_path =
        this->dictionaryPath_ + "/" +
        static_cast<std::string>(this->dictionary_["Data Sets"][this->dataset_]["Location"]) + "/" +
        fileName;
    Require(std::filesystem::exists(file_path));

    // Open and read individual collision file
    std::ifstream collisionFile(file_path);
    nlohmann::json collisionDictionary;
    collisionFile >> collisionDictionary;
    collisionFile.close();

    // Build reaction data
    rtt_cpt_data::CPT_Collision elasticCollision(collisionName, collisionDictionary,
                                                 collisionNumber, isHeavyCollision,
                                                 this->coarseGroupBounds_);

    // Add reaction to map
    collisionClasses_.insert(std::make_pair(collisionName, std::move(elasticCollision)));
  }
}
//------------------------------------------------------------------------------------------------//

#else

#endif

} // namespace rtt_cpt_data

//------------------------------------------------------------------------------------------------//
// End cpt_data/CPT_ElasticScattering.cc
//------------------------------------------------------------------------------------------------//
