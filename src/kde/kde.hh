//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/kde.hh
 * \author Mathew Cleveland
 * \brief  Define class kernel density estimator class
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
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
  // NESTED CLASSES AND TYPEDEFS

  // CREATORS

  // ACCESSORS

  // SERVICES

  //! Reconstruct distribution
  std::vector<double> reconstruction(const std::vector<double> &distribution,
                                     const std::vector<std::array<double, 3>> &one_over_band_width,
                                     const quick_index &qindex,
                                     const double discontinuity_cutoff = 0.0) const;

  //! Reconstruct distribution in logarithmic space
  std::vector<double>
  log_reconstruction(const std::vector<double> &distribution,
                     const std::vector<std::array<double, 3>> &one_over_band_width,
                     const quick_index &qindex, const double discontinuity_cutoff = 0.0) const;

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

protected:
  // IMPLEMENTATION

private:
  // NESTED CLASSES AND TYPEDEFS

  // IMPLEMENTATION

  // DATA
};

} // end namespace rtt_kde

#include "kde.i.hh"

#endif // kde_kde_hh

//------------------------------------------------------------------------------------------------//
// end of kde/kde.hh
//------------------------------------------------------------------------------------------------//
