//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   compton_tools/Compton_Native.cc
 * \author Andrew Till
 * \date   11 May 2020
 * \brief  Implementation file for native Compton binary-read and temperature interpolation
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "compton_tools/Compton_Native.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include <array>
#include <fstream>

using UINT64 = uint64_t;
using FP = double;
using vec_d = std::vector<double>;

namespace rtt_compton_tools {

//------------------------------------------------------------------------------------------------//
// FREE-FLOATING HELPER FUNCTIONS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief find location in a sorted list and min/max value to be within list
 *
 * \param[in]     xs vector of monotonically increasing and unique values
 * \param[in,out] x  value whose location in xs needs to be found; modified so that xs[index] <= x
 *                   <= xs[index+1]
 * \return Index in [0, xs.size()-2] so that xs[index] and xs[index+1] are valid
 *
 * Use a binary search to find the location of x in xs
 */
size_t find_index(const vec_d &xs, double &x) {
  Require(xs.size() >= 2);

  // Do binary search for index
  auto loc = std::upper_bound(xs.begin(), xs.end(), x);
  // loc - begin guaranteed to be in [0,len]
  size_t index = (loc - xs.cbegin());

  // Move index and x to be interior
  const size_t len = xs.size();
  const double xmin = xs[0];
  const double xmax = xs[len - 1U];
  index = std::min(len - 1U, std::max(static_cast<size_t>(1), index)) - 1U;
  x = std::max(xmin, std::min(xmax, x));
  return index;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief compute Hermite polynomial for a given value and the left/right grid points
 *
 * \param[in] x value of independent variable
 * \param[in] xL grid point to the left of x (xL <= x)
 * \param[in] xR grid point to the right of x (x <= xR)
 * \return length-4 array of Hermite polynomials x for generic interpolation in x
 *
 * If vL/vR are function values at xL/xR and dL/dR are function derivatives at xL/xR, then
 * \code
 * H[0] * vL + H[1] * vR + H[2] * dL + H[3] * dR
 * \endcode
 * interpolates the value of the function at x, where this function computes and returns H
 */
template <typename T_FP> std::array<T_FP, 4> hermite(T_FP x, T_FP xL, T_FP xR) {
  Require(xL < xR);
  Require(xL <= x);
  Require(x <= xR);

  // Spacing of interval
  const T_FP dx = xR - xL;
  // Left/right basis functions
  const T_FP bL = (xR - x) / dx;
  const T_FP bR = (x - xL) / dx;

  // Hermite functions
  std::array<T_FP, 4> hermite_eval{bL * bL * (3 - 2 * bL), bR * bR * (3 - 2 * bR),
                                   -dx * bL * bL * (bL - 1), dx * bR * bR * (bR - 1)};
  return hermite_eval;
}

//------------------------------------------------------------------------------------------------//
// SETUP FUNCTIONS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for Compton_Native
 *
 * \param[in] filename Name of the binary CSK data file to be read
 *
 * Rank 0 reads the binary CSK data file, which fills in the class' data members.  The data members
 * are then broadcast to other MPI ranks to finish their construction.
 */
Compton_Native::Compton_Native(const std::string &filename) {
  Require(filename.length() > 0U);
  int rank = rtt_c4::node();
  constexpr int bcast_rank = 0;
  int errcode = 0;

  if (rank == bcast_rank) {
    errcode = read_binary(filename);
  }
  broadcast_MPI(errcode);

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
// BINARY READING AND BROADCASTING FUNCTIONS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Helper member function that broadcasts CSK data from rank 0 to all ranks and sets data
 *
 * \param[in] errcode If non-zero, everyone aborts.
 *
 * Uses rtt_c4's broadcast to send arrays and vectors to all ranks with the sparse CSK data,
 * temperature / energy grids, and sizes
 */
void Compton_Native::broadcast_MPI(int errcode) {
  rtt_c4::global_max(errcode);
  if (errcode == 1) {
    // Have all ranks throw a bad file exception
    std::ifstream f("");
    f.exceptions(f.failbit);
  }
  Insist(errcode == 0, "Non-zero errorcode. Exiting.");

  int rank = rtt_c4::node();
  constexpr size_t bcast_rank = 0;

  // Broadcast sizes
  size_t data_size = data_.size();
  std::array<size_t, 6> pack = {num_temperatures_, num_groups_, num_leg_moments_,
                                num_evals_,        num_points_, data_size};
  rtt_c4::broadcast(&pack[0], pack.size(), bcast_rank);
  size_t p = 0;
  num_temperatures_ = pack[p++];
  num_groups_ = pack[p++];
  num_leg_moments_ = pack[p++];
  num_evals_ = pack[p++];
  num_points_ = pack[p++];
  data_size = pack[p++];
  Insist(data_size < (1U << 31U), "Unrepresentable with int needed in the broadcast");

  // Broadcast grids
  if (rank != bcast_rank)
    Ts_.resize(num_temperatures_);
  rtt_c4::broadcast(Ts_.begin(), Ts_.end(), Ts_.begin());

  if (rank != bcast_rank)
    Egs_.resize(num_groups_ + 1U);
  rtt_c4::broadcast(Egs_.begin(), Egs_.end(), Egs_.begin());

  // Broadcast sparse data structures
  if (rank != bcast_rank)
    first_groups_.resize(num_temperatures_ * num_groups_);
  rtt_c4::broadcast(first_groups_.begin(), first_groups_.end(), first_groups_.begin());

  if (rank != bcast_rank)
    indexes_.resize(num_temperatures_ * num_groups_ + 1U);
  rtt_c4::broadcast(indexes_.begin(), indexes_.end(), indexes_.begin());

  // Broadcast data itself
  if (rank != bcast_rank)
    data_.resize(data_size);
  rtt_c4::broadcast(data_.begin(), data_.end(), data_.begin());

  if (rank != bcast_rank)
    derivs_.resize(data_size);
  rtt_c4::broadcast(derivs_.begin(), derivs_.end(), derivs_.begin());

  return;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Helper member function to read a binary CSK file and set class data
 *
 * \param[in] filename Path to CSK binary file
 * \return errcode Zero if read is successful, otherwise non-zero
 *
 * Reads a binary CSK file by interpreting the characters as 64-bit unsigned integers and doubles
 */
int Compton_Native::read_binary(const std::string &filename) {

  // Read
  auto fin = std::ifstream(filename, std::ios::in | std::ios::binary);

  // Check for valid file stream
  try {
    fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  } catch (std::fstream::failure & /*error*/) {
    std::cout << "Error: Compton data file '" << filename << "' not found!\n";
    return 1;
  }

  // Ensure valid type
  // (using vector of char to avoid using C-style strings)
  std::array<char, 6> expected = {' ', 'c', 's', 'k', ' ', '\0'};
  std::array<char, 6> actual{' ', ' ', ' ', ' ', ' ', ' '};
  for (char &c : actual)
    fin.read(&c, sizeof(char));
  if (!std::equal(expected.begin(), expected.end(), actual.begin())) {
    std::cerr << "Expecting binary file " << filename << " to start with '";
    for (char c : expected)
      std::cerr << c;
    std::cerr << "' but got '";
    for (char c : actual)
      std::cerr << c;
    std::cerr << "'" << std::endl;
    return 2;
  }

  UINT64 binary_ordering;
  UINT64 version_major;
  UINT64 version_minor;
  fin.read(reinterpret_cast<char *>(&version_major), sizeof(UINT64));
  fin.read(reinterpret_cast<char *>(&version_minor), sizeof(UINT64));
  fin.read(reinterpret_cast<char *>(&binary_ordering), sizeof(UINT64));
  if (version_major != 1U || binary_ordering > 1U) {
    std::cerr << "Expecting a CSK binary file (version 1) with ordering 0 or 1 but got "
              << version_major << " with ordering " << binary_ordering << std::endl;
    return 3;
  }

  constexpr size_t n = 7U;
  std::array<UINT64, n> szs{0};
  for (size_t i = 0; i < n; ++i)
    fin.read(reinterpret_cast<char *>(&szs[i]), sizeof(szs[i]));
  size_t j = 0;
  size_t tsz = szs[j++];
  size_t gsz = szs[j++];
  size_t lsz = szs[j++];
  size_t esz = szs[j++];
  size_t fgsz = szs[j++];
  size_t isz = szs[j++];
  size_t dsz = szs[j++];

  num_temperatures_ = tsz;
  num_groups_ = gsz;
  num_leg_moments_ = lsz;
  num_evals_ = esz;
  // point is (leg moment, eval) pair
  // first eval has all leg moments; others have only the 0th moment
  num_points_ = num_evals_ + num_leg_moments_ - 1U;
  size_t egsz = gsz + 1;

#if 0
    std::cout << "DBG num_temperatures_ " << num_temperatures_ << '\n';
    std::cout << "DBG num_groups_ " << num_groups_ << '\n';
    std::cout << "DBG num_leg_moments_ " << num_leg_moments_ << '\n';
    std::cout << "DBG num_evals_ " << num_evals_ << '\n';
    std::cout << "DBG num_points_ " << num_points_ << '\n';
    std::cout << "DBG len(first_groups_) " << fgsz << '\n';
    std::cout << "DBG len(indexes_) " << isz << '\n';
    std::cout << "DBG len(data_/derivs_) " << dsz << '\n';
#endif

  Ts_.resize(tsz);
  for (size_t i = 0; i < tsz; ++i) {
    // Convert from FP (type in the binary file) to double
    FP tmp;
    fin.read(reinterpret_cast<char *>(&tmp), sizeof(FP));
    Ts_[i] = tmp;
  }

  Egs_.resize(egsz);
  for (size_t i = 0; i < egsz; ++i) {
    FP tmp;
    fin.read(reinterpret_cast<char *>(&tmp), sizeof(FP));
    Egs_[i] = tmp;
  }

  first_groups_.resize(fgsz);
  for (size_t i = 0; i < fgsz; ++i) {
    // Convert from UINT (type in binary file) to size_t
    UINT64 tmp;
    fin.read(reinterpret_cast<char *>(&tmp), sizeof(UINT64));
    first_groups_[i] = tmp;
  }

  indexes_.resize(isz);
  for (size_t i = 0; i < isz; ++i) {
    UINT64 tmp;
    fin.read(reinterpret_cast<char *>(&tmp), sizeof(UINT64));
    indexes_[i] = tmp;
  }

  data_.resize(dsz);
  for (size_t i = 0; i < dsz; ++i) {
    FP tmp;
    fin.read(reinterpret_cast<char *>(&tmp), sizeof(FP));
    data_[i] = tmp;
  }

  derivs_.resize(dsz);
  for (size_t i = 0; i < dsz; ++i) {
    FP tmp;
    fin.read(reinterpret_cast<char *>(&tmp), sizeof(FP));
    derivs_[i] = tmp;
  }

  fin.close();

#if 0
    std::cout << '\n';
    std::cout << "-------------------------------------------------\n";

    std::cout << "\nDBG Ts";
    for (size_t i = 0; i < tsz; ++i)
      std::cout << ", " << Ts_[i] / 510.998;
    std::cout << '\n';

    std::cout << "\nDBG Egs";
    for (size_t i = 0; i < egsz; ++i)
      std::cout << ", " << Egs_[i] / 510.998;
    std::cout << '\n';

    std::cout << "\nDBG first_groups";
    for (size_t i = 0; i < fgsz; ++i)
      std::cout << ", " << first_groups_[i];
    std::cout << '\n';

    std::cout << "\nDBG indexes";
    for (size_t i = 0; i < isz; ++i)
      std::cout << ", " << indexes_[i];
    std::cout << '\n';

    std::cout << "\nDBG data";
    for (size_t p = 0; p < num_points_; ++p) {
      std::cout << "\nPoint " << p << "\n";
      size_t szp = dsz / num_points_;
      for (size_t ii = 0; ii < fgsz; ++ii) {
        size_t istrt = indexes_[ii] + p * szp;
        size_t iend = indexes_[ii + 1] + p * szp;
        std::cout << "  index " << ii;
        for (size_t i = istrt; i < iend; ++i) {
          std::cout << std::setprecision(12); // << std::scientific;
          std::cout << ", " << data_[i] / 0.075116337052433;
        }
        std::cout << '\n';
      }
    }
    std::cout << std::defaultfloat;

    std::cout << "-------------------------------------------------\n";
    std::cout << '\n';
#endif

  return 0;
}

//------------------------------------------------------------------------------------------------//
// TEMPERATURE INTERPOLATION FUNCTIONS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate csk data in temperature and return dense linear inscattering matrix
 *
 * \param[in,out] inscat The flattened (1D), dense inscattering matrix as a vector.  Order of inscat
 *                       is (slow) [moment, group-to, group-from] (fast) Does NOT need to be the
 *                       right size prior to calling
 * \param[in] Te_keV The electron temperature in keV at which the interpolation is desired
 * \param[in] num_moments_truncate The maximum number of Legendre moments to use. Function will use
 *                       the minimum of this variable and number of moments in the data
 */
void Compton_Native::interp_dense_inscat(vec_d &inscat, double Te_keV,
                                         size_t num_moments_truncate) const {
  // Finds index and nudges Teff such that
  // Ts_[index] <= Teff <= Ts_[index+1] and 0 <= index <= Ts_.size()-2;
  double Teff = Te_keV;
  size_t iT = rtt_compton_tools::find_index(Ts_, Teff);

  // Fill Hermite function
  std::array<double, 4> hermite = rtt_compton_tools::hermite<double>(Teff, Ts_[iT], Ts_[iT + 1U]);

  // Precompute some sparse indexes
  const size_t sz = indexes_[indexes_.size() - 1];
  const size_t end_leg = std::min(num_moments_truncate, num_leg_moments_);
  const size_t eval_offset = 0; // in_lin

  // Resize and fill with zeros
  inscat.resize(end_leg * num_groups_ * num_groups_);
  std::fill(inscat.begin(), inscat.end(), 0.0);

  // Apply Hermite function
  for (size_t k = 0; k < end_leg; ++k) {
    for (size_t gfrom = 0; gfrom < num_groups_; ++gfrom) {
      const size_t offset_jj = gfrom + num_groups_ * num_groups_ * k;
      // Get contributions from both Ts_[iT] and Ts_[iT+1]
      for (size_t n = 0; n < 2U; ++n) {
        const size_t i = gfrom + num_groups_ * (iT + n);
        const size_t first_gto = first_groups_[i];
        const size_t num_entries = indexes_[i + 1U] - indexes_[i];
        const size_t offset_ii = indexes_[i] + sz * k + eval_offset;
        for (size_t dg = 0; dg < num_entries; ++dg) {
          const size_t gto = dg + first_gto;
          const size_t ii = dg + offset_ii;
          const size_t jj = gto * num_groups_ + offset_jj;
          inscat[jj] += hermite[0U + n] * data_[ii] + hermite[2U + n] * derivs_[ii];
        }
      }
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate csk data in temperature and return linear outscattering vector
 *
 * \param[in,out] outscat The 1D linear outscattering array at the desired temperature. Has been
 *               summed over outgoing group so only index is [group-from]. Does NOT need to be the
 *               right size prior to calling.
 * \param[in] Te_keV The electron temperature in keV at which the interpolation is desired
 */
void Compton_Native::interp_linear_outscat(vec_d &outscat, double Te_keV) const {
  // Finds index and nudges Teff st Ts_[index] <= Teff <= Ts_[index+1]
  // and 0 <= index <= Ts_.size()-2;
  double Teff = Te_keV;
  size_t iT = rtt_compton_tools::find_index(Ts_, Teff);

  // Fill Hermite function
  std::array<double, 4> hermite = rtt_compton_tools::hermite<double>(Teff, Ts_[iT], Ts_[iT + 1U]);

  // Precompute some sparse indexes
  const size_t sz = indexes_[indexes_.size() - 1];
  const size_t eval_offset = sz * num_leg_moments_; // out_lin

  // Resize and fill with zeros
  outscat.resize(num_groups_);
  std::fill(outscat.begin(), outscat.end(), 0.0);

  // Apply Hermite function
  for (size_t gfrom = 0; gfrom < num_groups_; ++gfrom) {
    // Get contributions from both Ts_[iT] and Ts_[iT+1U]
    for (size_t n = 0; n < 2U; ++n) {
      const size_t i = gfrom + num_groups_ * (iT + n);
      const size_t num_entries = indexes_[i + 1U] - indexes_[i];
      const size_t offset = indexes_[i] + eval_offset;
      for (size_t dg = 0; dg < num_entries; ++dg) {
        const size_t ii = dg + offset;
        outscat[gfrom] += hermite[0U + n] * data_[ii] + hermite[2U + n] * derivs_[ii];
      }
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate csk data in temperature and add nonlinear difference to outscattering
 *
 * \param[in,out] outscat The 1D net outscattering array at the desired temperature
 *               Has been summed over outgoing group so only index is [group-from]
 *               MUST be the right size (# groups) and initialized with data prior to calling
 * \param[in] Te_keV The electron temperature in keV at which the interpolation is desired
 * \param[in] phi The multigroup radiation field of size number of groups
 * \param[in] scale The scale for phi. when the radiation is in equilibrium, sum_g phi_g = scale
 *
 * Adds difference (nonlinear outscattering minus nonlinear inscattering) to the outscattering
 * vector. The contribution is nonlinear because it depends on phi, the radiation field.
 * The use of scale allows phi to be passed in with arbitrary normalization (4pi, c, a, etc.).
 */
void Compton_Native::interp_nonlin_diff_and_add(vec_d &outscat, double Te_keV, const vec_d &phi,
                                                double scale) const {
  // Adds to existing outscat vector
  Require(outscat.size() == num_groups_);
  Require(phi.size() == num_groups_);

  // Finds index and nudges Teff st Ts_[index] <= Teff <= Ts_[index+1]
  // and 0 <= index <= Ts_.size()-2;
  double Teff = Te_keV;
  size_t iT = rtt_compton_tools::find_index(Ts_, Teff);

  // Fill Hermite function
  std::array<double, 4> hermite = rtt_compton_tools::hermite<double>(Teff, Ts_[iT], Ts_[iT + 1U]);

  // Precompute some sparse indexes
  const size_t sz = indexes_[indexes_.size() - 1];
  const size_t eval_offset = sz * (num_leg_moments_ + 1U); // nl_diff

  // Precompute constants
  const double invscale = scale > 0.0 ? 1.0 / scale : 0.0;

  // Apply Hermite function
  for (size_t gfrom = 0; gfrom < num_groups_; ++gfrom) {
    // Get contributions from both Ts_[iT] and Ts_[iT+1U]
    for (size_t n = 0; n < 2U; ++n) {
      const size_t i = gfrom + num_groups_ * (iT + n);
      const size_t first_gto = first_groups_[i];
      const size_t num_entries = indexes_[i + 1U] - indexes_[i];
      const size_t offset = indexes_[i] + eval_offset;
      for (size_t dg = 0; dg < num_entries; ++dg) {
        const size_t gto = dg + first_gto;
        const double mag = invscale * phi[gto];
        const size_t ii = dg + offset;
        const double val = hermite[0U + n] * data_[ii] + hermite[2U + n] * derivs_[ii];
        outscat[gfrom] += mag * val;
      }
    }
  }
}

//------------------------------------------------------------------------------------------------//
// UNIMPLEMENTED, POTENTIAL FUTURE FUNCTIONS
//------------------------------------------------------------------------------------------------//

#if 0
void Compton_Native::interp_matvec(vec_d &x, const vec_d &leftscale,
                             const vec_d &rightscale, double Te_keV,
                             bool zeroth_moment_only) const {
  // TODO: Redo interface? Not in-place??

  const vec_d &L = leftscale;
  const vec_d &R = rightscale;
  double const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  size_t const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;
  // TODO: implement
}

//------------------------------------------------------------------------------------------------//

void Compton_Native::interp_matvec_transpose(vec_d &xT, const vec_d &leftscale,
                                       const vec_d &rightscale, double Te_keV,
                                       bool zeroth_moment_only) const {
  // TODO: Redo interface? Not in-place??

  const vec_d &L = leftscale;
  const vec_d &R = rightscale;
  double const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  size_t const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;
  // TODO: implement
}

//------------------------------------------------------------------------------------------------//

void Compton_Native::interp_sparse_inscat(Sparse_Compton_Matrix &inscat,
                                    const vec_d &leftscale, const vec_d &rightscale,
                                    double Te_keV,
                                    bool zeroth_moment_only) const {
  // TODO: Allow offset to fill in data?

  const vec_d &L = leftscale;
  const vec_d &R = rightscale;
  double const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  size_t const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;
  // TODO: implement
}

//------------------------------------------------------------------------------------------------//

void Compton_Native::interp_dense_inscat(vec_d &inscat, const vec_d &leftscale,
                                   const vec_d &rightscale, double Te_keV,
                                   bool zeroth_moment_only) const {
  // Ordering of inscat is 1D array (slow) [moment, group-to, group-from] (fast)
  const vec_d &L = leftscale;
  const vec_d &R = rightscale;
  size_t const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;

  inscat.resize(num_groups_ * num_groups_ * end_leg);
  std::fill(inscat.begin(), inscat.end(), 0.0);

  // TODO: implement
}

//------------------------------------------------------------------------------------------------//

void Compton_Native::interp_linear_outscat(vec_d &outscat, const vec_d &leftscale,
                                     const vec_d &rightscale,
                                     double Te_keV) const {
  outscat.resize(num_groups_);
  std::fill(outscat.begin(), outscat.end(), 0.0);

  const vec_d &L = leftscale;
  const vec_d &R = rightscale;
  double const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  // TODO: implement
}

//------------------------------------------------------------------------------------------------//

void Compton_Native::interp_nonlinear_diff(vec_d &nldiff, const vec_d &leftscale,
                                     const vec_d &rightscale, const vec_d &flux,
                                     double flux_scale, double Te_keV) const {
  // Need "const" that changes based on DBC level?
  // Require(leftscale.len() == num_groups_);
  // Require(rightscale.len() == num_groups_);
  // Require(flux.len() == num_groups_);
  // Require(flux_scale > 0.0);
  // Require(Te_keV >= 0.0);

  nldiff.resize(num_groups_);
  std::fill(nldiff.begin(), nldiff.end(), 0.0);

  const vec_d &L = leftscale;
  const vec_d &R = rightscale;
  double const fscale = flux_scale;
  double const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  // TODO: implement

  // Ensure(nldiff.len() == num_groups_);
}
#endif
//------------------------------------------------------------------------------------------------//

} // namespace rtt_compton_tools

//------------------------------------------------------------------------------------------------//
// End compton_tools/Compton_Native.cc
//------------------------------------------------------------------------------------------------//
