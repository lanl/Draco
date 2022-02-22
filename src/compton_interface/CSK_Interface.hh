//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   compton_interface/CSK_Interface.hh
 * \author Kendra Keady
 * \date   Mon Feb 27 2017
 * \brief  Header file for compton_interface CSK_generator interface
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_compton_interface_CSK_Interface_hh
#define rtt_compton_interface_CSK_Interface_hh

#include "compton_interface/config.h"

#ifdef COMPTON_FOUND

// C++ standard library dependencies
#include <iostream>
#include <memory>
#include <vector>

namespace rtt_compton_interface {

//================================================================================================//
/*!
 * \class CSK_Interface
 *
 * \brief Provides access to relativistic CSK_Interface scattering angle and
 *        multigroup frequency distributions from the CSK_generator project.
 *
 * This interface class allows the client to:
 * 1) access (interpolate) data from existing multigroup CSK_generator libraries
 * 2) build new multigroup libraries from existing CSK_generator pointwise libraries
 * 3) obtain auxiliary information for existing multigroup libraries (electron temperature bounds,
 *    frequency group structures, etc)
 *
 * This class is designed to be used with the CSK_generator library and headers.  If this is not
 * found at the CMake configure step, the lib_compton_interface portion of draco will not be built.
 *
 * \b User's \b environment
 *
 * CMake searches for the CSK_generator library/include headers during the configuration step. The
 * script that does this is located at:
 *
 * \c /draco/config/FindCOMPTON.cmake
 *
 * \example compton_interface/test/tCSK_Interface.cc
 *
 * This unit test demonstrates the two methods for constructing a CSK_Interface object, and
 * exercises all routines for interpolation and data access.
 */
//================================================================================================//
class CSK_Interface {

private:
  //! Shared pointer to an electron interpolation object:
  std::unique_ptr<etemp_interp> ei;

  //! Shared pointer to an llnl interpolation object:
  std::unique_ptr<llnl_interp> llnli;

public:
  //! Constructor for an existing multigroup library
  explicit CSK_Interface(const std::string &filehandle, const bool llnl_style = false);

  //! Constructor to build a multigroup library from an existing pointwise file
  CSK_Interface(const std::string &file, const std::vector<double> &group_bounds,
                const std::string &opac_type, const std::string &wt_func, const bool induced,
                const bool det_bal = false, const size_t nxi = 0);

  //! Dtor - defined in .cc file to prevent incomplete_type issues
  ~CSK_Interface();

  //! Interpolation of all csk opacity data to a certain electron temperature:
  std::vector<std::vector<std::vector<std::vector<double>>>>
  interpolate_csk(const double etemp, const bool limit_grps = true) const;

  //! Interpolation of all nu_ratio data to an electron temperature:
  std::vector<std::vector<double>> interpolate_nu_ratio(const double etemp,
                                                        const bool limit_grps = true) const;

  //! Versions for use with on-the-fly temperature interpolation:
  //! Interpolate (E)xpected (R)elative (E)nergy (C)hange:
  double interpolate_erec(const double Tm, const double freq) const;
  //! Interpolation total sigma_s for a given temperature and freq:
  double interpolate_sigc(const double Tm, const double freq) const;

  //! Versions for use with already-interped-in-etemp data:
  //! Interpolate (E)xpected (R)elative (E)nergy (C)hange:
  double interpolate_cell_erec(const int64_t cell, const double freq) const;
  //! Interpolate total sigma_s for a given cell index and freq:
  double interpolate_cell_sigc(const int64_t cell, const double freq) const;

  //! Interpolate all CSK data in temperature (stored in CSK for each cell):
  void interpolate_precycle(const std::vector<double> &Tms, const std::vector<double> &dens) const;

  //! Retrieve group structure for the given library (in kev):
  std::vector<double> get_group_bounds() const { return ei->get_group_bds_kev(); }

  //! Retrieve min electron temperature for the given library:
  double get_min_etemp() const { return ei->get_min_etemp(); }

  //! Retrieve max electron temperature for the given library:
  double get_max_etemp() const { return ei->get_max_etemp(); }

  //! Retrieve number of groups in the given multigroup structure:
  size_t get_num_groups() const { return ei->get_n_grps(); }

  //! Retrieve number of angular moments/evaluation points in the lib data:
  size_t get_num_xi() const { return ei->get_n_xi_pts(); }

  //! Retrieve electron temperature eval points (diagnostic use)
  std::vector<double> get_etemp_pts() const { return ei->get_etemp_pts(); }
};
} // namespace rtt_compton_interface

#endif // COMPTON_FOUND

#endif // rtt_compton_interface_CSK_Interface_hh

//------------------------------------------------------------------------------------------------//
// End compton_interface/CSK_Interface.hh
//------------------------------------------------------------------------------------------------//
