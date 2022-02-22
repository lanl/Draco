//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   compton_tools/Compton_Native.hh
 * \author Andrew Till
 * \date   11 May 2020
 * \brief  Implementation file for native Compton binary-read and temperature interpolation
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_compton_tools_Compton_Native_hh
#define rtt_compton_tools_Compton_Native_hh

#include <cstdint>
#include <string>
#include <vector>
//#include "compton_tools/Sparse_Compton_Matrix.hh"

namespace rtt_compton_tools {

enum class Eval : size_t { in_lin = 0, out_nonlin = 1, nl_diff = 2 };

//================================================================================================//
/*!
 * \class Compton_Native
 *
 * \brief Provides access to relativistic Compton scattering angle and multigroup frequency
 *        distributions from CSK data files
 *
 * This interface class allows the client to:
 * - access (interpolate) data from existing multigroup CSK libraries
 * - obtain auxiliary information from existing multigroup libraries (electron temperature bounds,
 *   frequency group structures, etc)
 */

/*!
 * \example compton_tools/test/tCompton_Native.cc
 *
 * This unit test demonstrates the method for constructing a Compton object, and exercises all
 * routines for interpolation and data access.
*/
//================================================================================================//
class Compton_Native {

public:
  //----------------------------------------------------------------------------------------------//
  // SETUP FUNCTIONS
  //----------------------------------------------------------------------------------------------//

  // \brief Constructor
  explicit Compton_Native(const std::string &filename);

  //----------------------------------------------------------------------------------------------//
  // TEMPERATURE INTERPOLATION FUNCTIONS
  //----------------------------------------------------------------------------------------------//

  // Interpolate CSK data in temperature for the linear inscattering and return a flattened dense
  // matrix with ordering (slowest) [Legendre moment, group to, group from] (fastest). Must specify
  // the largest desired Legendre moment.
  void interp_dense_inscat(std::vector<double> &inscat, double Te_keV,
                           size_t num_moments_truncate) const;

  // Interpolate CSK data in temperature for the linear outscattering and sum over outgoing group to
  // return a vector with index [group from]. \c outscat must be pre-sized to the number of groups.
  void interp_linear_outscat(std::vector<double> &outscat, double Te_keV) const;

  // Adds the nonlinear difference to the \c outscat vector phi is a group-sized vector for the
  // scalar radiation intensity when the radiation is in equilibrium, sum_g phi_g = scale \c outscat
  // must be pre-sized and pre-filled with the \c outscat data by calling interp_linear_outscat
  void interp_nonlin_diff_and_add(std::vector<double> &outscat, double Te_keV,
                                  const std::vector<double> &phi, double scale) const;

  //----------------------------------------------------------------------------------------------//
  // UNIMPLEMENTED, POTENTIAL FUTURE FUNCTIONS
  //----------------------------------------------------------------------------------------------//

#if 0
  // Interpolate in temperature and multiply against a vector in-place
  // x := diag(leftscale) * csk[in_lin](Te_keV) * diag(rightscale) * x
  void interp_matvec(std::vector<double> &x,
                     const std::vector<double> &leftscale,
                     const std::vector<double> &rightscale, double Te_keV,
                     bool zeroth_moment_only = false) const;

  // Interpolate in temperature and transpose-multiply against a vector in-place
  // xT := xT * diag(leftscale) * csk[in_lin](Te_keV) * diag(rightscale)
  void interp_matvec_transpose(std::vector<double> &xT,
                               const std::vector<double> &leftscale,
                               const std::vector<double> &rightscale,
                               double Te_keV,
                               bool zeroth_moment_only = false) const;

  // Interpolate in temperature and return sparse in-scattering
  // inscat := diag(leftscale) * csk[in_lin](Te_keV) * diag(rightscale)
  void interp_sparse_inscat(Sparse_Compton_Matrix &inscat,
                            const std::vector<double> &leftscale,
                            const std::vector<double> &rightscale,
                            double Te_keV,
                            bool zeroth_moment_only = false) const;


  // Interpolate in temperature and return dense in-scattering
  // inscat := diag(leftscale) * csk[in_lin](Te_keV) * diag(rightscale)
  // (internally must transpose array)
  // ordering of inscat is 1D array (slow) [moment, group-to, group-from] (fast)
  void interp_dense_inscat(std::vector<double> &inscat,
                           const std::vector<double> &leftscale,
                           const std::vector<double> &rightscale, double Te_keV,
                           bool zeroth_moment_only = false) const;

  // Interpolate in temperature and return (linear) outscattering vector
  // outscat := onesT * diag(leftscale) * csk[out_lin](Te_keV) * diag(rightscale)
  // where onesT is a transposed array of ones; no induced effects
  void interp_linear_outscat(std::vector<double> &outscat,
                             const std::vector<double> &leftscale,
                             const std::vector<double> &rightscale,
                             double Te_keV) const;

  // Interpolate in temperature and return induced-scattering vector
  // nldiff := flux_scale * onesT * diag(leftscale) * csk[nl_diff](Te_keV) * diag(rightscale) * flux
  // flux should be of length number of groups
  // flux_scale should make flux sum to `a c T_r^4`
  void interp_nonlinear_diff(std::vector<double> &nldiff,
                             const std::vector<double> &leftscale,
                             const std::vector<double> &rightscale,
                             const std::vector<double> &flux, double flux_scale,
                             double Te_keV) const;
#endif

  //----------------------------------------------------------------------------------------------//
  // GETTER AND SIZE FUNCTIONS
  //----------------------------------------------------------------------------------------------//

  // Accessor functions
  size_t get_num_temperatures() const { return num_temperatures_; }
  size_t get_num_groups() const { return num_groups_; }
  size_t get_num_leg_moments() const { return num_leg_moments_; }
  size_t get_num_evals() const { return num_evals_; }
  size_t get_num_points() const { return num_points_; }
  size_t get_highest_leg_moment() const { return num_leg_moments_ - 1U; }
  const std::vector<double> &get_Ts() const { return Ts_; }
  const std::vector<double> &get_Egs() const { return Egs_; }

  //----------------------------------------------------------------------------------------------//
  // Size checks for valid state
  bool check_class_invariants() const {
    bool all_good = (num_temperatures_ > 0U) && (num_groups_ > 0U) && (num_leg_moments_ > 0U) &&
                    (num_evals_ >= 2U) && (num_evals_ <= 3U) &&
                    (num_points_ == num_evals_ + num_leg_moments_ - 1U) &&
                    (Ts_.size() == num_temperatures_) && (Egs_.size() == num_groups_ + 1U) &&
                    (first_groups_.size() == (num_temperatures_ * num_groups_)) &&
                    (indexes_.size() == (num_temperatures_ * num_groups_ + 1U)) &&
                    (data_.size() == derivs_.size()) &&
                    (data_.size() >= num_temperatures_ * num_groups_ * num_points_) &&
                    (data_.size() <= num_temperatures_ * num_groups_ * num_groups_ * num_points_) &&
                    true;
    return all_good;
  }

private:
  //----------------------------------------------------------------------------------------------//
  // PRIVATE DATA
  //----------------------------------------------------------------------------------------------//

  size_t num_temperatures_{0U}; //<! Number of temperature evaluations for CSK data
  size_t num_groups_{0U};       //<! Number of energy groups for CSK data
  size_t num_leg_moments_{0U};  //<! Number of Legendre moments for CSK data
  size_t num_evals_{0U};        //<! Number of "evaluations" (linear/nonlinear, in/out scattering)

  // a point is a (Legendre moment, evaluation) pair
  // first evaluation has all Legendre moments and all others have only the 0th moment
  size_t num_points_{0U};

  // Temperature grid for CSK data (keV)
  std::vector<double> Ts_{};

  // Energy grid (MG energy boundaries) for CSK data (keV)
  std::vector<double> Egs_{};

  // sparse data storage

  // first group-to with nonzero value
  // 1D array of [temperature, group-from]
  std::vector<size_t> first_groups_{};

  // cumulative sum of row offsets into data_ and derivs_
  // 1D array of [temperature, group-from]
  std::vector<size_t> indexes_{};

  // CSK data
  // 1D array of [eval, moment, temperature, group-from, group-to]
  std::vector<double> data_{};

  // temperature derivatives of CSK data
  // 1D array of [eval, moment, temperature, group-from, group-to]
  std::vector<double> derivs_{};

  //----------------------------------------------------------------------------------------------//
  // PRIVATE HELPER FUNCTIONS
  //----------------------------------------------------------------------------------------------//

  // broadcast data over MPI
  void broadcast_MPI(int errcode);

  // read the binary file
  int read_binary(const std::string &filename);
};

} // namespace rtt_compton_tools

#endif // rtt_compton_tools_Compton_Native_hh

//------------------------------------------------------------------------------------------------//
// End compton_tools/Compton_Native.hh
//------------------------------------------------------------------------------------------------//
