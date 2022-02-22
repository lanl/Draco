//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/NDI_AtomicMass.hh
 * \author Ben R. Ryan
 * \date   2020 Mar 6
 * \brief  NDI_AtomicMass class definition.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef cdi_ndi_NDI_AtomicMass_hh
#define cdi_ndi_NDI_AtomicMass_hh

#include "cdi_ndi/config.h" // definition of NDI_FOUND
#include "ds++/Assert.hh"
#include "ds++/Query_Env.hh"
#include "ds++/path.hh"
#include "units/PhysicalConstexprs.hh"
#include <string>

namespace rtt_cdi_ndi {

//================================================================================================//
/*!
 * \class NDI_AtomicMass
 *
 * \brief Class for getting atomic mass weights by ZAID from NDI data using a method due to
 *        T. Saller. For more details on NDI, see
 *        https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html Currently only multigroup data is
 *        supported, continuous energy data is probably best added through a refactor.
 * \example cdi_ndi/test/tstNDI_AtomicMass.cc
 *
 * Upon contruction, warn if the NDI library version is different that the NDI gendir version.
 * Assume versions are compatible for differences in the patch version.
 */
//================================================================================================//
class NDI_AtomicMass {
public:
  //! Default constructor
  explicit NDI_AtomicMass(
      std::string gendir_path_in = rtt_dsxx::get_env_val<std::string>("NDI_GENDIR_PATH").second);

  //! Retrieve atomic mass weight for isotope with given ZAID
  double get_amw(const int zaid) const;

private:
  //! Path to gendir file
  std::string const gendir_path;

  //! Unit system
  rtt_units::PhysicalConstexprs<rtt_units::CGS> const pc;
};

} // namespace rtt_cdi_ndi

#endif // cdi_ndi_NDI_AtomicMass_hh

//------------------------------------------------------------------------------------------------//
// End cdi_ndi/NDI_AtomicMass.hh
//------------------------------------------------------------------------------------------------//
