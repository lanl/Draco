//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/kde.cc
 * \author Mathew Cleveland
 * \date   November 10th 2020
 * \brief  Explicitly defined KDE functions for various dimensions and coordinate KDE or Kernel
 *         Density Estimators are unbiased statical based reconstruction.  They can significantly
 *         increase the convergence rate of statical distributions. The KDE performs a
 *         reconstruction by evaluating a mean over some discrete kernel shape. In this DRACO
 *         implementation the mean is evaluated based on the sample locations that are bound by the
 *         kernel shape.  A renormalization is used to ensure the proper mean is returned given
 *         there is no guarantee the full kernel (which integrates exactly to 1) will be integrated
 *         fully in space. This renormalization also avoids the need for boundary fix-ups which are
 *         typically used in KDE applications to account for the kernel extending beyond the bounds
 *         of the spatial domain. Other approaches that could be considered are quadrature based
 *         approaches that fully sample the Kernel space reducing the need for the normalization.
 *
 * \note   Copyright (C) 2021-2021 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "kde.hh"
#include <cmath>
#include <iostream>
#include <numeric>

namespace rtt_kde {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Calculate Weight
 * 
 * \pre Calculate the effective weight in Cartesian and Spherical reconstructions from a given location to the current
 * kernel 
 *
 * \param[in] r0 current kernel center location
 * \param[in] one_over_h0 current kernel width
 * \param[in] r data location
 * \param[in] one_over_h kernel width at this data location
 * \param[in] qindex quick indexing class
 * \param[in] discontinuity_cutoff maximum size of value discrepancies to include in the
 * reconstruction
 *
 * \return weight contribution to the current kernel
 *
 * \post the local reconstruction of the original data is returned.
 */
double kde::calc_weight(const std::array<double, 3> &r0, const std::array<double, 3> &one_over_h0,
                        const std::array<double, 3> &r, const std::array<double, 3> &one_over_h,
                        const quick_index &qindex, const double &discontinuity_cutoff) const {
  Require(one_over_h0[0] > 0.0);
  Require(qindex.dim > 1 ? one_over_h0[1] > 0.0 : true);
  Require(qindex.dim > 2 ? one_over_h0[2] > 0.0 : true);
  Require(one_over_h[0] > 0.0);
  Require(qindex.dim > 1 ? one_over_h[1] > 0.0 : true);
  Require(qindex.dim > 2 ? one_over_h[2] > 0.0 : true);
  // do not allow spherical reflection of the radial direction
  Require(qindex.spherical ? !reflect_boundary[0] : true);
  Require(qindex.spherical ? !reflect_boundary[1] : true);
  double weight = 1.0;
  // set the radius at which all arch_lengths should be computed
  // this value is ignored by calc_orthogonal_distance in non-spherical geometry
  double arch_radius = r0[0];
  std::array<double, 3> distance = qindex.calc_orthogonal_distance(r0, r, arch_radius);
  std::array<double, 3> low_reflect_r0_distance =
      qindex.calc_orthogonal_distance(qindex.bounding_box_min, r0, arch_radius);
  std::array<double, 3> low_reflect_r_distance =
      qindex.calc_orthogonal_distance(qindex.bounding_box_min, r, arch_radius);
  std::array<double, 3> high_reflect_r0_distance =
      qindex.calc_orthogonal_distance(r0, qindex.bounding_box_max, arch_radius);
  std::array<double, 3> high_reflect_r_distance =
      qindex.calc_orthogonal_distance(r, qindex.bounding_box_max, arch_radius);
  for (size_t d = 0; d < qindex.dim; d++) {
    const double u = distance[d] * one_over_h0[d];
    const double scale =
        fabs(one_over_h0[d] - one_over_h[d]) / std::max(one_over_h0[d], one_over_h[d]) >
                discontinuity_cutoff
            ? 0.0
            : 1.0;
    // Apply Boundary Condition Weighting
    double bc_weight = 1.0;
    const bool low_reflect = reflect_boundary[d * 2];
    const bool high_reflect = reflect_boundary[d * 2 + 1];
    if (low_reflect) {
      const double low_u =
          (low_reflect_r0_distance[d] + low_reflect_r_distance[d]) * one_over_h0[d];
      bc_weight += epan_kernel(low_u);
    }
    if (high_reflect) {
      const double high_u =
          (high_reflect_r0_distance[d] + high_reflect_r_distance[d]) * one_over_h0[d];
      bc_weight += epan_kernel(high_u);
    }
    weight *= scale * bc_weight * epan_kernel(u) * one_over_h0[d];
  }
  Ensure(!(weight < 0.0));
  return weight;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief KDE reconstruction 
 * 
 * \pre The local reconstruction data is passed into this function which includes the original data
 * distribution, its spatial position, and the optimal bandwidth to be used at each point.
 *
 * \param[in] distribution original data to be reconstructed
 * \param[in] one_over_bandwidth inverse bandwidth size to be used at each data location
 * \param[in] qindex quick_index class to be used for data access.
 * \param[in] discontinuity_cutoff maximum size of value discrepancies to include in the
 * reconstruction
 * \return final local KDE function distribution reconstruction
 *
 * \post the local reconstruction of the original data is returned.
 */
std::vector<double>
kde::reconstruction(const std::vector<double> &distribution,
                    const std::vector<std::array<double, 3>> &one_over_bandwidth,
                    const quick_index &qindex, const double discontinuity_cutoff) const {
  Require(qindex.dim < 3 && qindex.dim > 0);
  const size_t local_size = distribution.size();
  // be sure that the quick_index matches this data size
  Require(qindex.locations.size() == local_size);
  Require(one_over_bandwidth.size() == local_size);

  // used for the zero accumulation conservation
  std::vector<double> result(local_size, 0.0);
  std::vector<double> normal(local_size, 0.0);
  if (qindex.domain_decomposed) {

    std::vector<double> ghost_distribution(qindex.local_ghost_buffer_size);
    qindex.collect_ghost_data(distribution, ghost_distribution);
    std::vector<std::array<double, 3>> ghost_one_over_bandwidth(qindex.local_ghost_buffer_size,
                                                                {0.0, 0.0, 0.0});
    qindex.collect_ghost_data(one_over_bandwidth, ghost_one_over_bandwidth);

    std::array<double, 3> win_min{0.0, 0.0, 0.0};
    std::array<double, 3> win_max{0.0, 0.0, 0.0};
    // now apply the kernel to the local ranks
    for (size_t i = 0; i < local_size; i++) {
      const std::array<double, 3> r0 = qindex.locations[i];
      const std::array<double, 3> one_over_h0 = one_over_bandwidth[i];
      calc_win_min_max(qindex, r0, one_over_h0, win_min, win_max);
      const std::vector<size_t> coarse_bins = qindex.window_coarse_index_list(win_min, win_max);
      // fetch local contribution
      for (auto &cb : coarse_bins) {
        // skip bins that aren't present in the map (for constness)
        auto mapItr = qindex.coarse_index_map.find(cb);
        if (mapItr != qindex.coarse_index_map.end()) {
          // loop over local data
          for (auto &l : mapItr->second) {
            const double weight = calc_weight(r0, one_over_h0, qindex.locations[l],
                                              one_over_bandwidth[l], qindex, discontinuity_cutoff);
            result[i] += distribution[l] * weight;
            normal[i] += weight;
          }
        }
        auto gmapItr = qindex.local_ghost_index_map.find(cb);
        if (gmapItr != qindex.local_ghost_index_map.end()) {
          // loop over ghost data
          for (auto &g : gmapItr->second) {
            const double weight =
                calc_weight(r0, one_over_h0, qindex.local_ghost_locations[g],
                            ghost_one_over_bandwidth[g], qindex, discontinuity_cutoff);
            result[i] += ghost_distribution[g] * weight;
            normal[i] += weight;
          }
        }
      }
    }
  } else { // local reconstruction only

    std::array<double, 3> win_min{0.0, 0.0, 0.0};
    std::array<double, 3> win_max{0.0, 0.0, 0.0};
    // now apply the kernel to the local ranks
    for (size_t i = 0; i < local_size; i++) {
      const std::array<double, 3> r0 = qindex.locations[i];
      const std::array<double, 3> one_over_h0 = one_over_bandwidth[i];
      calc_win_min_max(qindex, r0, one_over_h0, win_min, win_max);
      const std::vector<size_t> coarse_bins = qindex.window_coarse_index_list(win_min, win_max);
      for (auto &cb : coarse_bins) {
        // skip bins that aren't present in the map (can't use [] operator with constness)
        auto mapItr = qindex.coarse_index_map.find(cb);
        if (mapItr != qindex.coarse_index_map.end()) {
          // loop over local data
          for (auto &l : mapItr->second) {
            const double weight = calc_weight(r0, one_over_h0, qindex.locations[l],
                                              one_over_bandwidth[l], qindex, discontinuity_cutoff);
            result[i] += distribution[l] * weight;
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
  }

  return result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief KDE reconstruction done in logarithmic data space
 * 
 * \pre The local reconstruction data is passed into this function which includes the original data
 * distribution, its spatial position, and the optimal bandwidth to be used at each point. The
 * original data distribution is transformed into log space prior and post reconstruction. This is
 * helpful for strongly peaked data and should be exact for exponential distributions.
 *
 * \param[in] distribution original data to be reconstructed
 * \param[in] one_over_bandwidth inverse bandwidth size to be used at each data location
 * \param[in] qindex quick_index class to be used for data access.
 * \param[in] discontinuity_cutoff maximum size of value discrepancies to include in the
 * reconstruction
 * \return final local KDE function distribution reconstruction
 *
 * \post the local reconstruction of the original data is returned.
 */
std::vector<double>
kde::log_reconstruction(const std::vector<double> &distribution,
                        const std::vector<std::array<double, 3>> &one_over_bandwidth,
                        const quick_index &qindex, const double discontinuity_cutoff) const {
  Require(qindex.dim < 3 && qindex.dim > 0);
  const size_t local_size = distribution.size();
  Require(qindex.locations.size() == local_size);
  Require(one_over_bandwidth.size() == local_size);

  // used for the zero accumulation conservation
  std::vector<double> result(local_size, 0.0);
  std::vector<double> normal(local_size, 0.0);
  double min_value = *std::min_element(distribution.begin(), distribution.end());
  double max_value = *std::max_element(distribution.begin(), distribution.end());
  double log_bias = fabs(min_value) + (max_value - min_value);
  if (qindex.domain_decomposed) {

    rtt_c4::global_min(min_value);
    rtt_c4::global_min(max_value);

    std::vector<double> ghost_distribution(qindex.local_ghost_buffer_size);
    qindex.collect_ghost_data(distribution, ghost_distribution);
    std::vector<std::array<double, 3>> ghost_one_over_bandwidth(qindex.local_ghost_buffer_size,
                                                                {0.0, 0.0, 0.0});
    qindex.collect_ghost_data(one_over_bandwidth, ghost_one_over_bandwidth);

    log_bias = fabs(min_value) + (max_value - min_value);
    // if the log bias is zero the answer must be zero everywhere
    if (!(log_bias > 0.0))
      return result;
    // now apply the kernel to the local ranks
    std::array<double, 3> win_min{0.0, 0.0, 0.0};
    std::array<double, 3> win_max{0.0, 0.0, 0.0};
    for (size_t i = 0; i < local_size; i++) {
      const std::array<double, 3> r0 = qindex.locations[i];
      const std::array<double, 3> one_over_h0 = one_over_bandwidth[i];
      calc_win_min_max(qindex, r0, one_over_h0, win_min, win_max);
      const std::vector<size_t> coarse_bins = qindex.window_coarse_index_list(win_min, win_max);
      // fetch local contribution
      for (auto &cb : coarse_bins) {
        // skip bins that aren't present in the map (can't use [] operator with constness)
        auto mapItr = qindex.coarse_index_map.find(cb);
        if (mapItr != qindex.coarse_index_map.end()) {
          // loop over local data
          for (auto &l : mapItr->second) {
            const double weight = calc_weight(r0, one_over_h0, qindex.locations[l],
                                              one_over_bandwidth[l], qindex, discontinuity_cutoff);
            result[i] += log_transform(distribution[l], log_bias) * weight;
            normal[i] += weight;
          }
        }
        auto gmapItr = qindex.local_ghost_index_map.find(cb);
        if (gmapItr != qindex.local_ghost_index_map.end()) {
          // loop over ghost data
          for (auto &g : gmapItr->second) {
            const double weight =
                calc_weight(r0, one_over_h0, qindex.local_ghost_locations[g],
                            ghost_one_over_bandwidth[g], qindex, discontinuity_cutoff);
            result[i] += log_transform(ghost_distribution[g], log_bias) * weight;
            normal[i] += weight;
          }
        }
      }
    }
  } else { // local reconstruction only

    // if the log bias is zero the answer must be zero everywhere
    if (!(log_bias > 0.0))
      return result;

    // now apply the kernel to the local ranks
    std::array<double, 3> win_min{0.0, 0.0, 0.0};
    std::array<double, 3> win_max{0.0, 0.0, 0.0};
    for (size_t i = 0; i < local_size; i++) {
      const std::array<double, 3> r0 = qindex.locations[i];
      const std::array<double, 3> one_over_h0 = one_over_bandwidth[i];
      calc_win_min_max(qindex, r0, one_over_h0, win_min, win_max);
      const std::vector<size_t> coarse_bins = qindex.window_coarse_index_list(win_min, win_max);
      // fetch local contribution
      for (auto &cb : coarse_bins) {
        // skip bins that aren't present in the map (can't use [] operator with constness)
        auto mapItr = qindex.coarse_index_map.find(cb);
        if (mapItr != qindex.coarse_index_map.end()) {
          // loop over local data
          for (auto &l : mapItr->second) {
            const double weight = calc_weight(r0, one_over_h0, qindex.locations[l],
                                              one_over_bandwidth[l], qindex, discontinuity_cutoff);
            result[i] += log_transform(distribution[l], log_bias) * weight;
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
  }

  return result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief KDE apply conservation
 * 
 * \pre Apply conservation fix to the new distribution so sum(original_distribution) ==
 * sum(new_distribution)
 *
 * \param[in] original_distribution original data to be reconstructed
 * \param[in,out] new_distribution original data to be reconstructed
 * \param[in] domain_decomposed bool
 *
 */
void kde::apply_conservation(const std::vector<double> &original_distribution,
                             std::vector<double> &new_distribution,
                             const bool domain_decomposed) const {

  const size_t local_size = original_distribution.size();
  Insist(new_distribution.size() == local_size,
         "Original and new distributions must be the same size");

  // compute absolute solution
  std::vector<double> abs_distribution(local_size, 0.0);
  for (size_t i = 0; i < local_size; i++) {
    if (!rtt_dsxx::soft_equiv(new_distribution[i], original_distribution[i], 1e-12))
      abs_distribution[i] = fabs(new_distribution[i]);
  }

  // compute totals to be used in residual calculation
  double original_conservation =
      std::accumulate(original_distribution.begin(), original_distribution.end(), 0.0);
  double reconstruction_conservation =
      std::accumulate(new_distribution.begin(), new_distribution.end(), 0.0);
  double abs_distribution_conservation =
      std::accumulate(abs_distribution.begin(), abs_distribution.end(), 0.0);

  if (domain_decomposed) {
    // accumulate global contribution
    rtt_c4::global_sum(original_conservation);
    rtt_c4::global_sum(reconstruction_conservation);
    rtt_c4::global_sum(abs_distribution_conservation);
  }

  // Apply residual
  if (abs_distribution_conservation > 0.0) {
    const double res = original_conservation - reconstruction_conservation;
    for (size_t i = 0; i < local_size; i++)
      new_distribution[i] += res * abs_distribution[i] / abs_distribution_conservation;
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Calculate window min and max bounds.
 *
 *  Calculate the bounding window (via win_min (x_min,y_min,z_min) and win_max (x_max, y_max,
 *  z_max)) given a central location and the bandwidth size in each dimension (dx,dy) for Cartesian
 *  or (dr,arc_length) for spherical.
 * 
 * \param[in] qindex quick index class for finding bounds xy bounds of a wedge shape
 * \param[in] position is the central location of the bounds
 * \param[in] one_over_bandwidth size of the reconstruction domain in each dimension. This is
 * (dx,dy) for Caresian and (dr, arc_length) for spherical. 
 * \param[in,out] win_min is the minimum corner of the bounding box (x_min, y_min, z_min)
 * \param[in,out] win_max is the maximum corner of the bounding box (x_max, y_max, z_max)
 *
 */
void kde::calc_win_min_max(const quick_index &qindex, const std::array<double, 3> &position,
                           const std::array<double, 3> &one_over_bandwidth,
                           std::array<double, 3> &win_min, std::array<double, 3> &win_max) const {
  const size_t dim = qindex.dim;
  Require(dim > 0);
  Require(one_over_bandwidth[0] > 0.0);
  Require(dim > 1 ? one_over_bandwidth[1] > 0.0 : true);
  Require(dim > 2 ? one_over_bandwidth[2] > 0.0 : true);
  for (size_t d = 0; d < dim; d++) {
    win_min[d] = position[d] - 1.0 / one_over_bandwidth[d];
    win_max[d] = position[d] + 1.0 / one_over_bandwidth[d];
  }
}

} // end namespace rtt_kde

//------------------------------------------------------------------------------------------------//
// end of kde/kde.cc
//------------------------------------------------------------------------------------------------//
