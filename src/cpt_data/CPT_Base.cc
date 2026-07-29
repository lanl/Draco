//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cpt_data/CPT_Base.cc
 * \author Kyle S. Beling
 * \brief  CPT_Base member definitions.
 * \note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cpt_data/CPT_Base.hh"
#include "dsxx/Assert.hh"
#include <string>
#include <utility>

#ifdef JSON_FOUND
#include <filesystem>
#include <fstream>
#endif

namespace rtt_cpt_data {

//================================================================================================//
// CONSTRUCTORS
//================================================================================================//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for generic JSON reader
 *
 * \param[in] dictionaryPathIn path to dictionary file
 * \param[in] datasetIn name of requested temperature dataset
 * \param[in] libraryIn name of requested library (provided by inherited class)
 */
CPT_Base::CPT_Base(std::string dictionaryPathIn, std::string datasetIn, std::string libraryIn)
    : dictionaryPath_(std::move(dictionaryPathIn)), dataset_(std::move(datasetIn)),
      library_(std::move(libraryIn)) {

#ifdef JSON_FOUND
  Require(std::filesystem::exists(dictionaryPath_));

  Require(dictionaryPath_.length() > 0);
  Require(dataset_.length() > 0);
  Require(library_.length() > 0);

  // Open and read dictionary
  std::string const dictionary_file_name = this->dictionaryPath_ + "/dictionary.json";
  std::ifstream dictionary_file(dictionary_file_name);
  dictionary_file >> dictionary_;
  dictionary_file.close();

  // Check dictionary data
  check_dictionary_data();

#else
  Insist(0, "CPT Data library is only available when JSON is found.");
#endif
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for generic JSON reader with no specified dataset
 *
 * \param[in] dictionaryPathIn path to dictionary file
 * \param[in] libraryIn name of requested library (provided by inherited class)
 */
CPT_Base::CPT_Base(std::string dictionaryPathIn, std::string libraryIn)
    : dictionaryPath_(std::move(dictionaryPathIn)), dataset_("Default"),
      library_(std::move(libraryIn)) {
#ifdef JSON_FOUND
  Require(std::filesystem::exists(dictionaryPath_));

  Require(dictionaryPath_.length() > 0);
  Require(dataset_.length() > 0);
  Require(library_.length() > 0);

  // Open and read dictionary
  std::string const dictionary_file_name = this->dictionaryPath_ + "/dictionary.json";
  std::ifstream dictionary_file(dictionary_file_name);
  dictionary_file >> dictionary_;
  dictionary_file.close();

  // Check dictionary data
  check_dictionary_data();

#else
  Insist(0, "CPT Data library is only available when JSON is found.");
#endif
}
//------------------------------------------------------------------------------------------------//

//================================================================================================//
// MEMBER FUNCTIONS
//================================================================================================//

// Only implemented if JSON is found
#ifdef JSON_FOUND

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Check the dictionary.json file for the necessary keys and that the required sub-dictionary 
 *        files exist.
 */
void CPT_Base::check_dictionary_data() {
  // Check if dictionary_ contains a valid format
  Insist(this->dictionary_.contains("Data Sets"), "Key: \"Data Set\" not found in dictionary.json");
  Insist(this->dictionary_.contains("Reaction Networks"),
         "Key: \"Reaction Networks\" not found in dictionary.json");
  Insist(this->dictionary_.contains("Elastic Scattering"),
         "Key: \"Elastic Scattering\" not found in dictionary.json");
  Insist(this->dictionary_.contains("TN Reactions"),
         "Key: \"TN Reactions\" not found in dictionary.json");

  // Check if dictionary_ contains the user specified dataset_
  Insist(this->dictionary_["Data Sets"].contains(this->dataset_),
         "Specified dataset not found in availible \"Data Sets\"");

  // Check reaction network file exists
  this->reactionNetworkFileName_ =
      this->dictionaryPath_ + "/" +
      static_cast<std::string>(this->dictionary_["Reaction Networks"]["File"]);
  Insist(std::filesystem::exists(this->reactionNetworkFileName_),
         "Reaction Networks file not found.");

  // Check TN reaction info file exists
  this->tnReactionFileName_ = this->dictionaryPath_ + "/" +
                              static_cast<std::string>(this->dictionary_["TN Reactions"]["File"]);
  Insist(std::filesystem::exists(this->tnReactionFileName_), "TN Reactions file not found.");

  // Check elastic scattering info file exists
  this->elasticScatteringFileName_ =
      this->dictionaryPath_ + "/" +
      static_cast<std::string>(this->dictionary_["Elastic Scattering"]["File"]);
  Insist(std::filesystem::exists(this->elasticScatteringFileName_),
         "Elastic Scattering file not found.");
}
//------------------------------------------------------------------------------------------------//

#else

#endif

} // namespace rtt_cpt_data

//------------------------------------------------------------------------------------------------//
// End cpt_data/CPT_Base.cc
//------------------------------------------------------------------------------------------------//
