//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/NDI_AtomicMass.cc
 * \author Ben R. Ryan
 * \date   2020 Mar 6
 * \brief  NDI_AtomicMass class declaration.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "NDI_AtomicMass.hh"
#include "NDI_Base.hh"
#include "ds++/DracoStrings.hh"
#include "ds++/SystemCall.hh"
#include <array>

namespace rtt_cdi_ndi {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for NDI atomic mass weight reader, using custom path to NDI gendir file.
 * \param[in] gendir_path_in path to gendir file
 *
 * Print a warning if the gendir version and the ndi library version are not compatible.
 */
NDI_AtomicMass::NDI_AtomicMass(std::string gendir_path_in)
    : gendir_path(std::move(gendir_path_in)), pc() {
  Insist(rtt_dsxx::fileExists(gendir_path),
         "Specified NDI library is not available. gendir_path = " + gendir_path);
  NDI_Base::warn_ndi_version_mismatch(gendir_path);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Get atomic mass weight of an isotope with given ZAID. Use method due to T. Saller that
 *        invokes multigroup_neutron dataset which includes atomic weights.
 *
 * \pre This function requires gendir_path to be valid.  If it isn't valid this NDI_atomicMass
 *      object will fail to contruct at run time.
 *
 * \param[in] zaid ZAID of isotope for which to return the atomic mass.
 * \return mass of isotope in grams
 */
double NDI_AtomicMass::get_amw(const int zaid) const {
  // The below NDI method doesn't work for electrons
  if (zaid == -1) {
    return pc.electronMass();
  }

#ifdef NDI_FOUND
  int gendir_handle = -1;
  int ndi_error = -9999;
  ndi_error = NDI2_open_gendir(&gendir_handle, gendir_path.c_str());
  Require(ndi_error == 0);
  Insist(gendir_handle != -1, "gendir_handle still has default value!");

  ndi_error = NDI2_set_option_gendir(gendir_handle, NDI_LIB_TYPE_DEFAULT, "multigroup_neutron");
  Require(ndi_error == 0);

  ndi_error = NDI2_set_option_gendir(gendir_handle, NDI_LIBRARY_DEFAULT, "mendf71x");
  Require(ndi_error == 0);

  std::string zaid_formatted = std::to_string(zaid) + ".";

  int size = NDI2_get_size_x(gendir_handle, NDI_AT_WGT, zaid_formatted.c_str(), &ndi_error);
  Require(ndi_error == 0);
  Insist(size == 1, "NDI returned more or fewer than one atomic weight?");

  std::array<double, 1> arr{0.0};
  ndi_error =
      NDI2_get_float64_vec_x(gendir_handle, NDI_AT_WGT, zaid_formatted.c_str(), arr.data(), size);
  Require(ndi_error == 0);

  ndi_error = NDI2_close_gendir(gendir_handle);
  Require(ndi_error == 0);

  return arr[0] * pc.amu();
#else
  // NDI gendir not available.
  Insist(rtt_dsxx::fileExists(gendir_path),
         "Specified NDI library is not available. gendir_path = " + gendir_path);
  return 0.0;
#endif
}

} // namespace rtt_cdi_ndi

//------------------------------------------------------------------------------------------------//
// End cdi_ndi/NDI_AtomicMass.cc
//------------------------------------------------------------------------------------------------//
