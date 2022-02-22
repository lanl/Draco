//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   compton_interface/CSK_Interface.cc
 * \author Kendra Keady
 * \date   Tues Feb 21 2017
 * \brief  Implementation file for compton CSK_generator interface
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

// headers provided in draco:
#include "compton_interface/CSK_Interface.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"

#ifdef COMPTON_FOUND

namespace rtt_compton_interface {

// ------------ //
// Constructors //
// ------------ //

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for an existing multigroup libfile.
 *
 * This calls CSK_generator methods to read the data file and store everything in a CSK_Interface
 * data object, a smart pointer to which is then passed to (and held by) the CSK_generator
 * etemp_interp class.
 *
 * \param[in] filehandle The name of the Compton multigroup file
 * \param[in] llnl_style Defaults to false. True indicates that data uses LLNL format.
 */
CSK_Interface::CSK_Interface(const std::string &filehandle, const bool llnl_style) {

  // Make a compton_interface file object to read the multigroup data
  compton_file Cfile(false);
  if (llnl_style) {
    // initialize the etemp/frequency interpolated with the library data:
    llnli = std::make_unique<llnl_interp>(Cfile.read_llnl_data(filehandle));
    // Make sure the SP exists...
    Ensure(llnli);
  } else {
    // initialize the electron temperature interpolator with the mg compton data
    ei = std::make_unique<etemp_interp>(Cfile.read_mg_csk_data(filehandle));
    // Make sure the SP exists...
    Ensure(ei);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for an existing pointwise file and a multigroup structure.
 *
 * In the .hh file, we default the number of angular evals (n_xi) to zero This causes the CSK
 * routines to use the full angular fidelity of the library if no n_xi argument is passed
 *
 * This calls CSK_generator methods to read the pointwise library and construct a multigroup
 * CSK_Interface data object, a smart pointer to which is then passed to (and held by) the
 * CSK_generator etemp_interp class.
 *
 * \param[in] filehandle The name of the pointwise lib to build MG data from
 * \param[in] grp_bds    A vector containing the multigroup bounds (in keV)
 * \param[in] opac_type The type of opacity to build. Valid options for CSK v0.3 are "jayenne" (for
 *                       IMC-style opacities) or "capsaicin" (for Sn-style opacities). Any other
 *                       string will cause CSK to throw an exception
 * \param[in] wt_func The frequency weighting function used to numerically integrate the
 *                       opacities. Valid options for CSK v0.2 are "flat", "wien" or "planck." Any
 *                       other string will cause CSK to throw an exception.
 * \param[in] induced    Bool to toggle consideration of induced effects off/on
 * \param[in] det_bal    Bool to toggle detailed balance enforcement off/on
 * \param[in] nxi        The number of angular points/Legendre moments desired
 */
CSK_Interface::CSK_Interface(const std::string &filehandle, const std::vector<double> &grp_bds,
                             const std::string &opac_type, const std::string &wt_func,
                             const bool induced, const bool det_bal, const size_t nxi) {

  // Check input validity
  Require(std::ifstream(filehandle).good());
  Require(grp_bds.size() > 0);

  if (rtt_c4::node() == 0) {
    std::cout << "*********************************************************\n"
              << "WARNING! Building a multigroup library from scratch might\n"
              << " take a LOOOOOOONG time! (Don't say I didn't warn you.)  \n"
              << "*********************************************************\n"
              << std::endl;
  }

  // do quick sanity check
  if (det_bal) {
    // if we're enforcing detailed balance, we need induced and wt_func to be set to
    // <0,"wien">||<1,"planck">; these are the only valid cases
    Insist(((!induced && wt_func == std::string("wien")) ||
            (induced && wt_func == std::string("planck"))),
           "CSK_Interface error: Detailed balance enforcement (det_bal = 1) \n"
           "only valid for induced=0 w/wien -OR- induced=1 w/planck!");
  }

  // make a group_data struct to pass to the lib builder:
  multigroup::Group_data grp_data = {multigroup::Library_type::EXISTING,
                                     multigroup::string_to_opac_type(opac_type),
                                     multigroup::string_to_wt_func(wt_func),
                                     induced,
                                     det_bal,
                                     filehandle,
                                     nxi,
                                     grp_bds};

  // Construct a multigroup library builder:
  multigroup_lib_builder MG_builder(grp_data, rtt_c4::node());

  // build the library:
  MG_builder.build_library();

  // initialize the electron temperature interpolator with the mg compton data
  ei = std::make_unique<etemp_interp>(MG_builder.package_data());

  // Make sure the SP exists...
  Ensure(ei);
}

// Default destructor.
CSK_Interface::~CSK_Interface() = default;

// ------------ //
//  Interfaces  //
// ------------ //

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate opacity data to a given SCALED electron temperature \f$ (T / m_e) \f$
 *
 * This method interpolates MG CSK_Interface opacity data to a given electron temperature. It
 * returns the interpolated values for ALL g, g', and angular points in the specified multigroup
 * structure.
 *
 * \param[in] etemp The SCALED electron temperature (temp / electron rest-mass)
 * \param[in] limit_grps When true, CSK attempts to reduce the energy domain considered.  This can
 *                  reduce required memory and compute resources (default = true).
 * \return   n_opac x n_grp x n_grp x n_xi interpolated opacity values
 */
std::vector<std::vector<std::vector<std::vector<double>>>>
CSK_Interface::interpolate_csk(const double etemp, const bool limit_grps) const {

  // Be sure the passed electron temperature is within the bounds of the lib!
  Require(etemp >= ei->get_min_etemp());
  Require(etemp <= ei->get_max_etemp());

  // call the appropriate routine in the electron interp object
  return ei->interpolate_csk(etemp, limit_grps);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate nu_ratio data to a given SCALED electron temperature (T / m_e)
 *
 * This method interpolates MG Compton nu_ratio data to a given electron temperature. It returns the
 * interpolated values for ALL g and g' points in the specified multigroup structure.
 *
 * \param[in] etemp The SCALED electron temperature (temp / electron rest-mass)
 * \param[in] limit_grps When true, CSK attempts to reduce the energy domain considered.  This can
 *                  reduce required memory and compute resources (default = true).
 * \return    n_grp x n_grp interpolated nu_ratio values
 */
std::vector<std::vector<double>> CSK_Interface::interpolate_nu_ratio(const double etemp,
                                                                     const bool limit_grps) const {

  // Be sure the passed electron temperature is within the bounds of the lib!
  Require(etemp >= ei->get_min_etemp());
  Require(etemp <= ei->get_max_etemp());

  // call the appropriate routine in the electron interp object
  return ei->interpolate_nu_ratio(etemp, limit_grps);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate EREC data to a given electron temperature / frequency
 *
 * This method uses data and routines in CSK to interpolate a value of the expected relative energy
 * change for some temperature / frequency
 *
 * \param[in] Tm   Electron temperature (temp / electron rest-mass)
 * \param[in] freq Incident frequency (keV)
 * \return    The interpolated relative energy change (Delta-E / E)
 */
double CSK_Interface::interpolate_erec(const double Tm, const double freq) const {
  // call the appropriate routine in the electron interp object
  // (unscaled -- it'll be scaled in the library
  return llnli->interpolate_erec(Tm, freq);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate Compton opacity data to a given electron temperature / frequency
 *
 * This method uses data and routines in CSK to interpolate a value of the Compton scattering
 * opacity for some temperature / frequency. The returned value will have units of cm^2/g, and must
 * be scaled by density for direct use in transport
 *
 * \param[in] Tm   Electron temperature (temp / electron rest-mass)
 * \param[in] freq Incident frequency (keV)
 * \return    The interpolated opacity (cm^2 / g)
 */
double CSK_Interface::interpolate_sigc(const double Tm, const double freq) const {
  // call the appropriate routine in the electron interp object
  // (unscaled -- it'll be scaled in the library
  return llnli->interpolate_sigc(Tm, freq);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate EREC data in a cell for a given frequency
 *
 * This method uses data and routines in CSK to interpolate a value of the expected relative energy
 * change for some cell index / frequency.  This call is only valid if the interpolate_precycle()
 * function has been called, which interpolates all opacity data to the cell temperatures
 * (otherwise, CSK will throw an error).
 *
 * \param[in] cell Cell index
 * \param[in] freq Incident frequency (keV)
 * \return    The interpolated relative energy change (Delta-E / E)
 */
double CSK_Interface::interpolate_cell_erec(const int64_t cell, const double freq) const {
  // call the appropriate routine in the electron interp object (unscaled -- it'll be scaled in the
  // library)
  Require(llnli->pre_interped());
  return llnli->interpolate_erec(cell, freq);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate Compton opacity data in a cell for a given frequency
 *
 * This method uses data and routines in CSK to interpolate a value of the Compton scattering
 * opacity for some cell index / frequency. The returned value will have units of cm^-1. This call
 * is only valid if the interpolate_precycle() function has been called, which interpolates all EREC
 * data to the cell temperatures (otherwise, CSK will throw an error).
 *
 * \param[in] cell Cell index
 * \param[in] freq Incident frequency (keV)
 * \return    The interpolated opacity (cm^-1)
 */
double CSK_Interface::interpolate_cell_sigc(const int64_t cell, const double freq) const {
  // call the appropriate routine in the electron interp object (unscaled -- it'll be scaled in the
  // library)
  Require(llnli->pre_interped());
  return llnli->interpolate_sigc(cell, freq);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate opacity and EREC data to cell temperatures
 *
 * This function passes the cell temperatures and densities to CSK before a transport cycle. The
 * opacity and EREC data is then "pre-interpolated" in electron temperature, so it can later be
 * referenced by cell index.
 *
 * \param[in] Tms  Cell electron temperature (keV)
 * \param[in] dens Cell densities (g/cc)
 */
void CSK_Interface::interpolate_precycle(const std::vector<double> &Tms,
                                         const std::vector<double> &dens) const {
  llnli->preinterp_in_temp(Tms, dens);
}
} // namespace rtt_compton_interface

#endif

//------------------------------------------------------------------------------------------------//
// End compton_interface/CSK_Interface.cc
//------------------------------------------------------------------------------------------------//
