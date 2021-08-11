//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/kde.hh
 * \author Mathew Cleveland
 * \brief  Define class kernel density estimator class
 * \note   Copyright (C) 2020-2021 Triad National Security, LLC., All rights reserved. 
 */
//------------------------------------------------------------------------------------------------//

#ifndef kde_kde_hh
#define kde_kde_hh

#include "quick_index.hh"
#include "c4/global.hh"
#include <array>
#include <vector>

namespace rtt_kde {

//================================================================================================//
/*!
 * \class kde
 * \brief kernel density estimator class for generated smoothed reconstructions of point wise PDF
 *        data
 * 
 * \tparam coord enumeration specifying the KDE coordinate system to use.
 *
 * Returns a KDE reconstruction of a multidimensional distribution
 */
//================================================================================================//
class kde {
public:
  //! Constructor
  kde(const std::array<bool, 6> reflect_boundary_ = {false, false, false, false, false, false})
      : reflect_boundary(reflect_boundary_) {}

  //! Reconstruct distribution
  std::vector<double> reconstruction(const std::vector<double> &distribution,
                                     const std::vector<std::array<double, 3>> &one_over_band_width,
                                     const quick_index &qindex,
                                     const double discontinuity_cutoff = 1.0) const;

  //! Reconstruct distribution in logarithmic space
  std::vector<double>
  log_reconstruction(const std::vector<double> &distribution,
                     const std::vector<std::array<double, 3>> &one_over_band_width,
                     const quick_index &qindex, const double discontinuity_cutoff = 1.0) const;

  //! Apply conservation to reconstructed distribution
  void apply_conservation(const std::vector<double> &original_distribution,
                          std::vector<double> &new_distribution, const bool domain_decompsed) const;
  // STATICS

  //! Epanechikov Kernel
  inline double epan_kernel(const double x) const;

  //! Transform the solution into log space
  inline double log_transform(const double value, const double bias) const;

  //! Move the solution back from log space
  inline double log_inv_transform(const double log_value, const double bias) const;

  //! Calculate radius
  inline double calc_radius(const std::array<double, 3> &center,
                            const std::array<double, 3> &location) const;

  //! Calculate Arch Length between two locations at a specific radius
  inline double calc_arch_length(const std::array<double, 3> &center, const double radius,
                                 const std::array<double, 3> &location_1,
                                 const std::array<double, 3> &location_2) const;
  //! Setup a spherical reconstruction with a
  void set_sphere_center(const std::array<double, 3> &sph_center, const double min_radius,
                         const double max_radius) {
    Insist(max_radius > min_radius, "Spherical KDE max radius must be larger then min radius");

    use_spherical_reconstruction = true;
    sphere_center = sph_center;
    sphere_min_radius = min_radius;
    sphere_max_radius = max_radius;
  }

protected:
  // IMPLEMENTATION

private:
  //! Private function to calculate kernel weight
  double calc_weight(const std::array<double, 3> &r0, const std::array<double, 3> &one_over_h0,
                     const std::array<double, 3> &r, const std::array<double, 3> &one_over_h,
                     const quick_index &qindex, const double &discontinuity_cutoff) const {
    return use_spherical_reconstruction
               ? calc_spherical_weight(r0, one_over_h0, r, one_over_h, qindex, discontinuity_cutoff)
               : calc_cartesian_weight(r0, one_over_h0, r, one_over_h, qindex,
                                       discontinuity_cutoff);
  }

  void calc_win_min_max(const quick_index &qindex, const std::array<double, 3> &position,
                        const std::array<double, 3> &one_over_bandwidth, std::array<double, 3> &min,
                        std::array<double, 3> &max) const;

  double calc_spherical_weight(const std::array<double, 3> &r0,
                               const std::array<double, 3> &one_over_h0,
                               const std::array<double, 3> &r,
                               const std::array<double, 3> &one_over_h, const quick_index &qindex,
                               const double &discontinuity_cutoff) const;

  double calc_cartesian_weight(const std::array<double, 3> &r0,
                               const std::array<double, 3> &one_over_h0,
                               const std::array<double, 3> &r,
                               const std::array<double, 3> &one_over_h, const quick_index &qindex,
                               const double &discontinuity_cutoff) const;

  // DATA
  //! reflecting boundary conditions [lower_x, upper_x, lower_y, upper_y, lower_z, upper_z]
  const std::array<bool, 6> reflect_boundary;
  //! Spherical Mesh Reconstruction Data
  std::array<double, 3> sphere_center{0.0, 0.0, 0.0};
  double sphere_min_radius{0.0};
  double sphere_max_radius{0.0};
  bool use_spherical_reconstruction{false};
};

} // end namespace rtt_kde

#include "kde.i.hh"

#endif // kde_kde_hh

//------------------------------------------------------------------------------------------------//
// end of kde/kde.hh
//------------------------------------------------------------------------------------------------//
