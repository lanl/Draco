//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/NDI_Base.hh
 * \author Ben R. Ryan
 * \date   2019 Nov 4
 * \brief  NDI_Base class definition.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef cdi_ndi_NDI_Base_hh
#define cdi_ndi_NDI_Base_hh

#include "cdi_ndi/config.h" // Definition of NDI_FOUND

#include "ds++/Assert.hh"
#include "ds++/path.hh"
#include <algorithm>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace rtt_cdi_ndi {

//================================================================================================//
/*!
 * \class NDI_Base
 *
 * \brief Base class for wrapping NDI calls to NDI data.
 *
 * Reads data, constructs internal storage amenable to radiation calculations, and provides
 * accessors. Instantiated only through a data type-specific derived class. Energies and
 * temperatures are in units of keV. Velocity- averaged cross sections are in units of cm^3 sh^-1.
 * Probability density functions sum to unity. Unit conversions from NDI data are done when data is
 * initially read in. For more details on NDI, see
 * https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html Currently only multigroup data is
 * supported, continuous energy data is probably best added through a refactor.
 */
//================================================================================================//

class NDI_Base {

protected:
  //! Path to gendir file, which indexes an NDI dataset
  std::string gendir;

  //! Type of data to read (NDI supports multigroup_neutron, multigroup_photon, multigroup_multi,
  //! tn, tnreactions, and dosimetry_neutrons)
  const std::string dataset;

  //! Name of library in which to find reaction
  const std::string library;

  //! Constructors
  NDI_Base(std::string dataset_in, std::string library_in);
  NDI_Base(std::string gendir_in, std::string dataset_in, std::string library_in);

public:
  //! Disable default constructor
  NDI_Base() = delete;

  //! Destructor
  ~NDI_Base() = default;

  //! Disable copy constructor
  NDI_Base(const NDI_Base &rhs) = delete;

  //! Disable move constructor
  NDI_Base(NDI_Base &&rhs) noexcept = delete;

  //! Disable assignment
  NDI_Base &operator=(const NDI_Base &rhs) = delete;

  //! Disable move assignment
  NDI_Base &operator=(NDI_Base &&rhs) noexcept = delete;

  //! Get the name of the gendir file
  inline std::string get_gendir() const & { return gendir; }

  //! Get the dataset
  inline std::string get_dataset() const & { return dataset; }

  //! Get the library
  inline std::string get_library() const & { return library; }

  // >> Non-interacting helper functions.

  //! Helper function to format a warning message.
  static void warn_ndi_version_mismatch(std::string const &gendir);
};

} // namespace rtt_cdi_ndi

#endif // cdi_ndi_NDI_Base_hh

//------------------------------------------------------------------------------------------------//
// End cdi_ndi/NDI_Base.hh
//------------------------------------------------------------------------------------------------//
