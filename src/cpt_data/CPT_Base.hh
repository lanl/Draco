//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cpt_data/CPT_Base.hh
 * \author Kyle S. Beling
 * \date   Nov 19 2025
 * \brief  CPT_Base class definition.
 * \note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cpt_data_Base
#define rtt_cpt_data_Base

#include "cpt_data/config.h" // Definition of JSON_FOUND

#include "dsxx/Assert.hh"
#include "dsxx/path.hh"
#include <algorithm>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace rtt_cpt_data {
//================================================================================================//
/*!
 * \class CPT_Base
 * \brief Base class for wrapping JSON calls to light ion data.
 *
 * Base level class that opens and reads the dictionary.json file and tells the inheriting classes
 * where the data files are located.
 */
//================================================================================================//

class CPT_Base {

#ifdef JSON_FOUND
  //! Check dictionary.json object
  void check_dictionary_data();
#endif

protected:
  //! Path to dictionary file, which indexes a JSON dataset
  const std::string dictionaryPath_{};

  //! Name of dataset (Cutoff, BFP, etc.)
  const std::string dataset_{};

  //! Name of the library (RIF reactions, elastic scattering, icident neutron)
  const std::string library_{};

  //! Location and names of reaction, TN, and elastic scattering
  std::string reactionNetworkFileName_{}, tnReactionFileName_{}, elasticScatteringFileName_{};

  //! Constructors
  CPT_Base(std::string dictionaryPathIn, std::string datasetIn, std::string libraryIn);
  CPT_Base(std::string dictionaryPathIn, std::string libraryIn);

#ifdef JSON_FOUND
  nlohmann::json dictionary_;
#endif

public:
  //! Disable default constructor
  CPT_Base() = delete;

  //! Destructor
  ~CPT_Base() = default;

  //! Disable copy constructor
  CPT_Base(const CPT_Base &rhs) = delete;

  //! Disable move constructor
  CPT_Base(CPT_Base &&rhs) noexcept = delete;

  //! Disable assignment
  CPT_Base &operator=(const CPT_Base &rhs) = delete;

  //! Disable move assignment
  CPT_Base &operator=(CPT_Base &&rhs) noexcept = delete;

  //! Returns the location of the dictionary file
  [[nodiscard]] inline std::string get_dictionary_path() const & { return dictionaryPath_; }

  //! Returns the name of the dataset
  [[nodiscard]] inline std::string get_dataset() const & { return dataset_; }

  //! Returns the name of the library
  [[nodiscard]] inline std::string get_library() const & { return library_; }

#ifdef JSON_FOUND
  //! Returns the dictionary.json object
  [[nodiscard]] inline nlohmann::json get_dictionary() const & { return dictionary_; }
#endif
};

} // namespace rtt_cpt_data

#endif // cpt_data_Base

//------------------------------------------------------------------------------------------------//
// End cpt_data/CPT_Base.hh
//------------------------------------------------------------------------------------------------//
