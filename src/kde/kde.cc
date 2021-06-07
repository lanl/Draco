//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/kde.cc
 * \author Mathew Cleveland
 * \date   November 10th 2020
 * \brief  Explicitly defined KDE functions for various dimensions and coordinate
 *         KDE or Kernel Density Estimators are unbiased statical based
 *         reconstruction.  They can significantly increase the convergence
 *         rate of statical distributions. The KDE performs a reconstruction by
 *         evaluating a mean over some discrete kernel shape. In this DRACO
 *         implementation the mean is evaluated based on the sample locations
 *         that are bound by the kernel shape.  A renormalization is used to
 *         ensure the proper mean is returned given there is no guarantee the
 *         full kernel (which integrates exactly to 1) will be integrated fully
 *         in space. This renormalization also avoids the need for boundary
 *         fix-ups which are typically used in KDE applications to account for
 *         the kernel extending beyond the bounds of the spatial domain. Other
 *         approaches that could be considered are quadrature based approaches
 *         that fully sample the Kernel space reducing the need for the
 *         normalization.
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "kde.hh"
#include <cmath>
#include <iostream>
#include <numeric>

namespace rtt_kde {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief KDE reconstruction 
 * 
 * \pre The local reconstruction data is passed into this function which
 * includes the original data distribution, its spatial position, and the
 * optimal bandwidth to be used at each point.
 *
 * \param[in] distribution original data to be reconstructed
 * \param[in] one_over_bandwidth inverse bandwidth size to be used at each data location
 * \param[in] qindex quick_index class to be used for data access.
 * \param[in] discontinuity_cutoff size of value discrepancies to exclude from the reconstruction
 * \return final local KDE function distribution reconstruction
 *
 * \post the local reconstruction of the original data is returned.
 */
std::vector<double>
kde::reconstruction(const std::vector<double> &distribution,
                    const std::vector<std::array<double, 3>> &one_over_bandwidth,
                    const quick_index &qindex, const double discontinuity_cutoff) const {
  size_t dim = qindex.dim;
  Require(dim < 3 && dim > 0);
  const size_t local_size = distribution.size();
  // be sure that the quick_index matches this data size
  Require(qindex.locations.size() == local_size);
  Require(one_over_bandwidth.size() == local_size);

  // used for the zero accumulation conservation
  double global_conservation = std::accumulate(distribution.begin(), distribution.end(), 0.0);
  std::vector<double> result(local_size, 0.0);
  std::vector<double> abs_result(local_size, 0.0);
  std::vector<double> normal(local_size, 0.0);
  if (qindex.domain_decomposed) {

    rtt_c4::global_sum(global_conservation);

    std::vector<double> ghost_distribution = qindex.collect_ghost_data(distribution);
    std::vector<std::array<double, 3>> ghost_one_over_bandwidth =
        qindex.collect_ghost_data(one_over_bandwidth);

    // now apply the kernel to the local ranks
    for (size_t i = 0; i < local_size; i++) {
      const std::array<double, 3> r0 = qindex.locations[i];
      const std::array<double, 3> one_over_h = one_over_bandwidth[i];
      std::array<double, 3> win_min{0.0, 0.0, 0.0};
      std::array<double, 3> win_max{0.0, 0.0, 0.0};
      for (size_t d = 0; d < dim; d++) {
        Check(one_over_h[d] > 0.0);
        win_min[d] = r0[d] - 1.0 / one_over_h[d];
        win_max[d] = r0[d] + 1.0 / one_over_h[d];
      }
      const std::vector<size_t> coarse_bins = qindex.window_coarse_index_list(win_min, win_max);
      // fetch local contribution
      for (auto cbItr = coarse_bins.begin(); cbItr < coarse_bins.end(); cbItr++) {
        // skip bins that aren't present in the map (for constness)
        auto mapItr = qindex.coarse_index_map.find(*cbItr);
        if (mapItr != qindex.coarse_index_map.end()) {
          // loop over local data
          for (auto iItr = mapItr->second.begin(); iItr != mapItr->second.end(); iItr++) {
            double weight = 1.0;
            for (size_t d = 0; d < dim; d++) {
              const double r = qindex.locations[*iItr][d];
              const double u = (r0[d] - r) * one_over_h[d];
              const double scale =
                  one_over_h[d] / one_over_bandwidth[*iItr][d] < discontinuity_cutoff ? 0.0 : 1.0;
              weight *= scale * epan_kernel(u) * one_over_h[d];
            }
            result[i] += distribution[*iItr] * weight;
            normal[i] += weight;
          }
        }
        auto gmapItr = qindex.local_ghost_index_map.find(*cbItr);
        if (gmapItr != qindex.local_ghost_index_map.end()) {
          // loop over ghost data
          for (auto gItr = gmapItr->second.begin(); gItr != gmapItr->second.end(); gItr++) {
            double weight = 1.0;
            for (size_t d = 0; d < dim; d++) {
              const double r = qindex.local_ghost_locations[*gItr][d];
              const double u = (r0[d] - r) * one_over_h[d];
              const double scale =
                  one_over_h[d] / ghost_one_over_bandwidth[*gItr][d] < discontinuity_cutoff ? 0.0
                                                                                            : 1.0;
              weight *= scale * epan_kernel(u) * one_over_h[d];
            }
            result[i] += ghost_distribution[*gItr] * weight;
            normal[i] += weight;
          }
        }
      }
    }
  } else { // local reconstruction only

    // now apply the kernel to the local ranks
    for (size_t i = 0; i < local_size; i++) {
      const std::array<double, 3> r0 = qindex.locations[i];
      const std::array<double, 3> one_over_h = one_over_bandwidth[i];
      std::array<double, 3> win_min{0.0, 0.0, 0.0};
      std::array<double, 3> win_max{0.0, 0.0, 0.0};
      for (size_t d = 0; d < dim; d++) {
        Check(one_over_h[d] > 0.0);
        win_min[d] = r0[d] - 1.0 / one_over_h[d];
        win_max[d] = r0[d] + 1.0 / one_over_h[d];
      }
      const std::vector<size_t> coarse_bins = qindex.window_coarse_index_list(win_min, win_max);
      for (auto cbItr = coarse_bins.begin(); cbItr < coarse_bins.end(); cbItr++) {
        // skip bins that aren't present in the map (can't use [] operator with constness)
        auto mapItr = qindex.coarse_index_map.find(*cbItr);
        if (mapItr != qindex.coarse_index_map.end()) {
          // loop over local data
          for (auto iItr = mapItr->second.begin(); iItr != mapItr->second.end(); iItr++) {
            double weight = 1.0;
            for (size_t d = 0; d < dim; d++) {
              const double r = qindex.locations[*iItr][d];
              const double u = (r0[d] - r) * one_over_h[d];
              const double scale =
                  one_over_h[d] / one_over_bandwidth[*iItr][d] < discontinuity_cutoff ? 0.0 : 1.0;
              weight *= scale * epan_kernel(u) * one_over_h[d];
            }
            result[i] += distribution[*iItr] * weight;
            normal[i] += weight;
          }
        }
      }
    }
  }

  // normalize the integrated weight contributions
  for (size_t i = 0; i < local_size; i++) {
    Check(normal[i] > 0.0);
    result[i] /= normal[i];
    if (!rtt_dsxx::soft_equiv(result[i], distribution[i], 1e-12))
      abs_result[i] = fabs(result[i]);
  }

  double reconstruction_conservation = std::accumulate(result.begin(), result.end(), 0.0);
  double abs_reconstruction_conservation =
      std::accumulate(abs_result.begin(), abs_result.end(), 0.0);

  if (qindex.domain_decomposed) {
    // accumulate global contribution
    rtt_c4::global_sum(reconstruction_conservation);
    rtt_c4::global_sum(abs_reconstruction_conservation);
  }

  if (abs_reconstruction_conservation > 0.0) {
    const double res = global_conservation - reconstruction_conservation;
    for (size_t i = 0; i < local_size; i++)
      result[i] += res * abs_result[i] / abs_reconstruction_conservation;
  }

  return result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief KDE reconstruction done in logarithmic data space
 * 
 * \pre The local reconstruction data is passed into this function which
 * includes the original data distribution, its spatial position, and the
 * optimal bandwidth to be used at each point. The original data distribution
 * is transformed into log space prior and post reconstruction. This is helpful
 * for strongly peaked data and should be exact for exponential distributions.
 *
 * \param[in] distribution original data to be reconstructed
 * \param[in] one_over_bandwidth inverse bandwidth size to be used at each data location
 * \param[in] qindex quick_index class to be used for data access.
 * \param[in] discontinuity_cutoff size of value discrepancies to exclude from the reconstruction
 * \return final local KDE function distribution reconstruction
 *
 * \post the local reconstruction of the original data is returned.
 */
std::vector<double>
kde::log_reconstruction(const std::vector<double> &distribution,
                        const std::vector<std::array<double, 3>> &one_over_bandwidth,
                        const quick_index &qindex, const double discontinuity_cutoff) const {
  size_t dim = qindex.dim;
  Require(dim < 3 && dim > 0);
  const size_t local_size = distribution.size();
  Require(qindex.locations.size() == local_size);
  Require(one_over_bandwidth.size() == local_size);

  // used for the zero accumulation conservation
  double global_conservation = std::accumulate(distribution.begin(), distribution.end(), 0.0);
  std::vector<double> result(local_size, 0.0);
  std::vector<double> abs_result(local_size, 0.0);
  std::vector<double> normal(local_size, 0.0);
  double min_value = *std::min_element(distribution.begin(), distribution.end());
  double log_bias = fabs(min_value) * (1.0 + 1e-12);
  if (qindex.domain_decomposed) {

    rtt_c4::global_sum(global_conservation);
    rtt_c4::global_min(min_value);

    std::vector<double> ghost_distribution = qindex.collect_ghost_data(distribution);
    std::vector<std::array<double, 3>> ghost_one_over_bandwidth =
        qindex.collect_ghost_data(one_over_bandwidth);

    log_bias = fabs(min_value) * (1.0 + 1e-12);
    log_bias = std::max(log_bias, 1e-12);
    // now apply the kernel to the local ranks
    for (size_t i = 0; i < local_size; i++) {
      const std::array<double, 3> r0 = qindex.locations[i];
      const std::array<double, 3> one_over_h = one_over_bandwidth[i];
      std::array<double, 3> win_min{0.0, 0.0, 0.0};
      std::array<double, 3> win_max{0.0, 0.0, 0.0};
      for (size_t d = 0; d < dim; d++) {
        Check(one_over_h[d] > 0.0);
        win_min[d] = r0[d] - 1.0 / one_over_h[d];
        win_max[d] = r0[d] + 1.0 / one_over_h[d];
      }
      const std::vector<size_t> coarse_bins = qindex.window_coarse_index_list(win_min, win_max);
      // fetch local contribution
      for (auto cbItr = coarse_bins.begin(); cbItr < coarse_bins.end(); cbItr++) {
        // skip bins that aren't present in the map (can't use [] operator with constness)
        auto mapItr = qindex.coarse_index_map.find(*cbItr);
        if (mapItr != qindex.coarse_index_map.end()) {
          // loop over local data
          for (auto iItr = mapItr->second.begin(); iItr != mapItr->second.end(); iItr++) {
            double weight = 1.0;
            for (size_t d = 0; d < dim; d++) {
              const double r = qindex.locations[*iItr][d];
              const double u = (r0[d] - r) * one_over_h[d];
              const double scale =
                  one_over_h[d] / one_over_bandwidth[*iItr][d] < discontinuity_cutoff ? 0.0 : 1.0;
              weight *= scale * epan_kernel(u) * one_over_h[d];
            }
            result[i] += log_transform(distribution[*iItr], log_bias) * weight;
            normal[i] += weight;
          }
        }
        auto gmapItr = qindex.local_ghost_index_map.find(*cbItr);
        if (gmapItr != qindex.local_ghost_index_map.end()) {
          // loop over ghost data
          for (auto gItr = gmapItr->second.begin(); gItr != gmapItr->second.end(); gItr++) {
            double weight = 1.0;
            for (size_t d = 0; d < dim; d++) {
              const double r = qindex.local_ghost_locations[*gItr][d];
              const double u = (r0[d] - r) * one_over_h[d];
              const double scale =
                  one_over_h[d] / ghost_one_over_bandwidth[*gItr][d] < discontinuity_cutoff ? 0.0
                                                                                            : 1.0;
              weight *= scale * epan_kernel(u) * one_over_h[d];
            }
            result[i] += log_transform(ghost_distribution[*gItr], log_bias) * weight;
            normal[i] += weight;
          }
        }
      }
    }
  } else { // local reconstruction only

    log_bias = std::max(log_bias, 1e-12);
    // now apply the kernel to the local ranks
    for (size_t i = 0; i < local_size; i++) {
      const std::array<double, 3> r0 = qindex.locations[i];
      const std::array<double, 3> one_over_h = one_over_bandwidth[i];
      std::array<double, 3> win_min{0.0, 0.0, 0.0};
      std::array<double, 3> win_max{0.0, 0.0, 0.0};
      for (size_t d = 0; d < dim; d++) {
        Check(one_over_h[d] > 0.0);
        win_min[d] = r0[d] - 1.0 / one_over_h[d];
        win_max[d] = r0[d] + 1.0 / one_over_h[d];
      }
      const std::vector<size_t> coarse_bins = qindex.window_coarse_index_list(win_min, win_max);
      // fetch local contribution
      for (auto cbItr = coarse_bins.begin(); cbItr < coarse_bins.end(); cbItr++) {
        // skip bins that aren't present in the map (can't use [] operator with constness)
        auto mapItr = qindex.coarse_index_map.find(*cbItr);
        if (mapItr != qindex.coarse_index_map.end()) {
          // loop over local data
          for (auto iItr = mapItr->second.begin(); iItr != mapItr->second.end(); iItr++) {
            double weight = 1.0;
            for (size_t d = 0; d < dim; d++) {
              const double r = qindex.locations[*iItr][d];
              const double u = (r0[d] - r) * one_over_h[d];
              const double scale =
                  one_over_h[d] / one_over_bandwidth[*iItr][d] < discontinuity_cutoff ? 0.0 : 1.0;
              weight *= scale * epan_kernel(u) * one_over_h[d];
            }
            result[i] += log_transform(distribution[*iItr], log_bias) * weight;
            normal[i] += weight;
          }
        }
      }
    }
  }

  // normalize the integrated weight contributions
  for (size_t i = 0; i < local_size; i++) {
    Check(normal[i] > 0.0);
    result[i] /= normal[i];
    result[i] = log_inv_transform(result[i], log_bias);
    // ZERO IS ZERO AND THE LOG TRANSFORM CAN MAKE THE ZEROS NOT MATCH... SO FIX IT LIKE THIS
    if (rtt_dsxx::soft_equiv(result[i], 0.0) && rtt_dsxx::soft_equiv(distribution[i], 0.0))
      result[i] = distribution[i];
    if (!rtt_dsxx::soft_equiv(result[i], distribution[i], 1e-12))
      abs_result[i] = fabs(result[i]);
  }

  double reconstruction_conservation = std::accumulate(result.begin(), result.end(), 0.0);
  double abs_reconstruction_conservation =
      std::accumulate(abs_result.begin(), abs_result.end(), 0.0);

  if (qindex.domain_decomposed) {
    // accumulate global contribution
    rtt_c4::global_sum(reconstruction_conservation);
    rtt_c4::global_sum(abs_reconstruction_conservation);
  }

  if (abs_reconstruction_conservation > 0.0) {
    const double res = global_conservation - reconstruction_conservation;
    for (size_t i = 0; i < local_size; i++)
      result[i] += res * abs_result[i] / abs_reconstruction_conservation;
  }

  return result;
}
} // end namespace rtt_kde

//------------------------------------------------------------------------------------------------//
// end of kde/kde.cc
//------------------------------------------------------------------------------------------------//
