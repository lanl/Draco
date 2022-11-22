//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/NDI_Base.cc
 * \author Ben R. Ryan
 * \date   2020 Feb 4
 * \brief  NDI_Base member definitions.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "NDI_Base.hh"
#include "ds++/DracoStrings.hh"
#include "ds++/Query_Env.hh"
#include "ds++/SystemCall.hh"
#include "ds++/dbc.hh"
#include <algorithm>

namespace rtt_cdi_ndi {

//------------------------------------------------------------------------------------------------//
// Auxillary functions
//------------------------------------------------------------------------------------------------//

#ifdef NDI_FOUND
/*!
 * \brief Warn if NDI library version doesn't match GENDIR version to 2 digits.
 *
 * \param[in] gendir path to the gendir file
 *
 * This function is a no-op if \c NDI_FOUND is false.
 */
void NDI_Base::warn_ndi_version_mismatch(std::string const &gendir) {

  std::string gendir_ver = rtt_dsxx::extract_version(gendir, 2);
  std::string ndi_ver = rtt_dsxx::extract_version(NDI_ROOT_DIR, 2);
  if (gendir_ver.size() == 0)
    gendir_ver = "version unknown";
  if (ndi_ver.size() == 0)
    ndi_ver = "version unknown";
  if (gendir_ver != ndi_ver) {
    using DT = Term::DracoTerminal;
    std::cout << "\n"
              << Term::ccolor(DT::error) << "WARNING: In the cdi_ndi/NDI_Base constructor, the NDI "
              << "library version (" << ndi_ver << ") is different than the NDI GENDIR version ("
              << gendir_ver << "). \n"
              << Term::ccolor(DT::reset) << std::endl;
  }
}
#else
/*!
 * \brief Warn if NDI library version doesn't match GENDIR version to 2 digits.  No-op when
 *        NDI_FOUND is false. */
void NDI_Base::warn_ndi_version_mismatch(std::string const & /*gendir*/) {}
#endif

//------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//------------------------------------------------------------------------------------------------//

#ifndef NDI_FOUND

//================================================================================================//
// Stubbed implementation when NDI is unavailable
//================================================================================================//

//! Constructor for generic NDI reader- throws when NDI not available
NDI_Base::NDI_Base(std::string dataset_in, std::string library_in)
    : dataset(std::move(dataset_in)), library(std::move(library_in)) {
  Insist(0, "NDI default gendir path only available when NDI is found.");
}

#else

//================================================================================================//
// Normal implementation
//================================================================================================//

/*!
 * \brief Constructor for generic NDI reader, to be inherited by readers for generic gendir file
 *        path and dataset.
 *
 * This base constructor only sets some data members based on constructor input.  For more details
 * on NDI, see https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html
 *
 * \param[in] dataset_in name of requested dataset (provided by inherited class)
 * \param[in] library_in name of requested NDI data library
 */
NDI_Base::NDI_Base(std::string dataset_in, std::string library_in)
    : gendir(rtt_dsxx::get_env_val<std::string>("NDI_GENDIR_PATH").second),
      dataset(std::move(dataset_in)), library(std::move(library_in)) {

  Require(rtt_dsxx::fileExists(gendir));

  Require(gendir.length() > 0);
  Require(dataset.length() > 0);
  Require(library.length() > 0);
  warn_ndi_version_mismatch(gendir);
}

#endif

//================================================================================================//
// Normal implementation
//================================================================================================//

/*!
 * \brief Constructor for generic NDI reader, to be inherited by readers for specific gendir file
 *        path and dataset.
 *
 * This base constructor only sets some data members based on constructor input.  For more details
 * on NDI, see https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html
 *
 * \param[in] gendir_in path to non-standard NDI gendir file
 * \param[in] dataset_in name of requested dataset (provided by inherited class)
 * \param[in] library_in name of requested NDI data library
 */
NDI_Base::NDI_Base(std::string gendir_in, std::string dataset_in, std::string library_in)
    : gendir(std::move(gendir_in)), dataset(std::move(dataset_in)), library(std::move(library_in)) {

#ifdef NDI_FOUND
  Require(rtt_dsxx::fileExists(gendir));

  Require(gendir.length() > 0);
  Require(dataset.length() > 0);
  Require(library.length() > 0);
  warn_ndi_version_mismatch(gendir);
#else
  Insist(0, "NDI default gendir path only available when NDI is found.");
#endif
}

} // namespace rtt_cdi_ndi

//------------------------------------------------------------------------------------------------//
// End cdi_ndi/NDI_Base.cc
//------------------------------------------------------------------------------------------------//
