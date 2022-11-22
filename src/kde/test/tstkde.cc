//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/test/tstkde.cc
 * \author Mathew Cleveland
 * \date   Nov. 10th 2020
 * \brief  KDE function tests
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "kde/kde.hh"
#include "kde/quick_index.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include <numeric>

using namespace rtt_dsxx;
using namespace rtt_c4;
using namespace rtt_kde;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
//
void test_replication(ParallelUnitTest &ut) { // NOLINT [hicpp-function-size]
  kde test_kde;

  // test the epan kernel
  double value = test_kde.epan_kernel(0.0);
  if (!rtt_dsxx::soft_equiv(value, 0.75))
    ITFAILS;

  // spherical reconstruction
  {
    const bool spherical = true;
    const std::array<double, 3> sphere_center{0.0, -1.0, 0.0};
    kde sphere_kde;
    // reflect on the theta boundary
    kde theta_reflected_sphere_kde({false, false, true, true, false, false});
    const std::array<double, 8> radial_edges{0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0};
    const std::array<double, 9> cosine_edges{-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0};
    const size_t data_size = radial_edges.size() * cosine_edges.size();
    std::vector<std::array<double, 3>> position_array(data_size,
                                                      std::array<double, 3>{0.0, 0.0, 0.0});

    std::vector<double> shell_data(data_size, 0.0);
    std::vector<double> spoke_data(data_size, 0.0);
    std::vector<double> const_data(data_size, 1.0);
    std::vector<double> shell_weights(data_size, 0.0);
    size_t point_i = 0;
    size_t ri = 0;
    for (auto &r : radial_edges) {
      size_t mui = 0;
      for (auto &mu : cosine_edges) {
        spoke_data[point_i] = static_cast<double>(mui) + 1.0;
        shell_data[point_i] = static_cast<double>(ri) + 1.0;
        if (ri < 4)
          shell_weights[point_i] = 1e3;
        else
          shell_weights[point_i] = 1.0;
        double rel_y = r * mu;
        position_array[point_i][0] =
            rtt_dsxx::soft_equiv(r * r, rel_y * rel_y, 1e-6) ? 0.0 : sqrt(r * r - rel_y * rel_y);
        position_array[point_i][1] = sphere_center[1] + rel_y;
        point_i++;
        mui++;
      }
      ri++;
    }

    // zero reconstruction array
    {
      std::vector<double> zero_data(data_size, 0.0);
      // R-theta space bandwidths in spherical reconstruction
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          data_size, std::array<double, 3>{1.0, 1.0e12, 0.0});
      const bool dd = false;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 1.0;
      const size_t dim = 2;
      quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);

      std::vector<int> maskids(1, 1);
      std::vector<int> reconstruction_mask(data_size, 1);
      std::vector<double> smooth_result = sphere_kde.reconstruction(
          zero_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result = sphere_kde.log_reconstruction(
          zero_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      std::vector<double> sampled_smooth_result = sphere_kde.sampled_reconstruction(
          zero_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(zero_data, maskids, reconstruction_mask, smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(spoke_data, maskids, reconstruction_mask, log_smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(spoke_data, maskids, reconstruction_mask, sampled_smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < data_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], zero_data[i]))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], zero_data[i]))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(sampled_smooth_result[i], zero_data[i]))
          ITFAILS;
      }

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(std::accumulate(zero_data.begin(), zero_data.end(), 0.0),
                                std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(
              std::accumulate(zero_data.begin(), zero_data.end(), 0.0),
              std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(
              std::accumulate(zero_data.begin(), zero_data.end(), 0.0),
              std::accumulate(sampled_smooth_result.begin(), sampled_smooth_result.end(), 0.0)))
        ITFAILS;
    }
    // spoke reconstruction array
    {
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          data_size, std::array<double, 3>{1.0, 1.0e12, 0.0});
      const bool dd = false;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 1.0;
      const size_t dim = 2;
      quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);

      std::vector<int> maskids(1, 1);
      std::vector<int> reconstruction_mask(data_size, 1);
      std::vector<double> smooth_result = sphere_kde.reconstruction(
          spoke_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result = sphere_kde.log_reconstruction(
          spoke_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(spoke_data, maskids, reconstruction_mask, smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(spoke_data, maskids, reconstruction_mask, log_smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < data_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], spoke_data[i]))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], spoke_data[i]))
          ITFAILS;
      }

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
                                std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(
              std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
              std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
        ITFAILS;
    }

    // shell reconstruction array
    {
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          data_size, std::array<double, 3>{1.0e12, 1.0, 0.0});
      const bool dd = false;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 1.0;
      const size_t dim = 2;
      quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);

      std::vector<int> maskids(1, 1);
      std::vector<int> reconstruction_mask(data_size, 1);
      std::vector<double> smooth_result = sphere_kde.reconstruction(
          shell_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result = sphere_kde.log_reconstruction(
          shell_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(shell_data, maskids, reconstruction_mask, smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(shell_data, maskids, reconstruction_mask, log_smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < data_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], shell_data[i]))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], shell_data[i]))
          ITFAILS;
      }

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(std::accumulate(shell_data.begin(), shell_data.end(), 0.0),
                                std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(
              std::accumulate(shell_data.begin(), shell_data.end(), 0.0),
              std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
        ITFAILS;
    }

    // spoke smoothing on shell array
    {

      std::vector<double> spoke_smoothed_shells{
          3.7717,  3.7717,  3.7717,  3.7717,  3.7717,  3.7717,  3.7717,  3.7717,  3.7717,
          3.83452, 3.83452, 3.83452, 3.83452, 3.83452, 3.83452, 3.83452, 3.83452, 3.83452,
          3.89388, 3.89388, 3.89388, 3.89388, 3.89388, 3.89388, 3.89388, 3.89388, 3.89388,
          3.95025, 3.95025, 3.95025, 3.95025, 3.95025, 3.95025, 3.95025, 3.95025, 3.95025,
          4.24469, 4.24469, 4.24469, 4.24469, 4.24469, 4.24469, 4.24469, 4.24469, 4.24469,
          4.67528, 4.67528, 4.67528, 4.67528, 4.67528, 4.67528, 4.67528, 4.67528, 4.67528,
          5.21495, 5.21495, 5.21495, 5.21495, 5.21495, 5.21495, 5.21495, 5.21495, 5.21495,
          6.41471, 6.41472, 6.41472, 6.41472, 6.41471, 6.41471, 6.41471, 6.41471, 6.41471};
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          data_size, std::array<double, 3>{1.0, 1.0e12, 0.0});
      const bool dd = false;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 1.0;
      const size_t dim = 2;
      quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);

      std::vector<int> maskids(1, 1);
      std::vector<int> reconstruction_mask(data_size, 1);
      std::vector<double> smooth_result = sphere_kde.reconstruction(
          shell_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result = sphere_kde.reconstruction(
          shell_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(shell_data, maskids, reconstruction_mask, smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(shell_data, maskids, reconstruction_mask, log_smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < data_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], spoke_smoothed_shells[i], 1e-3))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], spoke_smoothed_shells[i], 1e-3))
          ITFAILS;
      }

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(std::accumulate(shell_data.begin(), shell_data.end(), 0.0),
                                std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(
              std::accumulate(shell_data.begin(), shell_data.end(), 0.0),
              std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
        ITFAILS;
    }

    // spoke smoothing on shell array with weights
    {
      // relative weights are used to reduce the smearing between shells <= 4 and shells > 4

      std::vector<double> spoke_weighted_smoothed_shells{
          2.49765, 2.49765, 2.49765, 2.49765, 2.49765, 2.49765, 2.49765, 2.49765, 2.49765,
          2.49922, 2.49922, 2.49922, 2.49922, 2.49922, 2.49922, 2.49922, 2.49922, 2.49922,
          2.50078, 2.50078, 2.50079, 2.50078, 2.50078, 2.50078, 2.50078, 2.50078, 2.50078,
          2.50235, 2.50235, 2.50235, 2.50236, 2.50235, 2.50235, 2.50235, 2.50235, 2.50235,
          6.2,     6.2,     6.2,     6.2,     6.2,     6.2,     6.2,     6.2,     6.2,
          6.41379, 6.41379, 6.41379, 6.41379, 6.41379, 6.41379, 6.41379, 6.41379, 6.41379,
          6.58621, 6.58621, 6.58621, 6.58621, 6.58621, 6.58621, 6.58621, 6.58621, 6.58621,
          6.8,     6.8,     6.8,     6.8,     6.8,     6.8,     6.8,     6.8,     6.8};
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          data_size, std::array<double, 3>{1.0, 1.0e12, 0.0});
      const bool dd = false;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 1.0;
      const size_t dim = 2;
      quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);

      std::vector<int> maskids(1, 1);
      std::vector<int> reconstruction_mask(data_size, 1);
      std::vector<double> smooth_result = sphere_kde.weighted_reconstruction(
          shell_data, shell_weights, reconstruction_mask, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(shell_data, maskids, reconstruction_mask, smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < data_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], spoke_weighted_smoothed_shells[i], 1e-3))
          ITFAILS;
      }

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(std::accumulate(shell_data.begin(), shell_data.end(), 0.0),
                                std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
        ITFAILS;
    }

    // shell smoothing on spoke array
    {
      std::vector<double> shell_smoothed_spoke{
          1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
          1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
          1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
          1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
          1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
          1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
          1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
          1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675};
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          data_size, std::array<double, 3>{1.0e12, 1.0, 0.0});
      const bool dd = false;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 1.0;
      const size_t dim = 2;
      quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);

      std::vector<int> maskids(1, 1);
      std::vector<int> reconstruction_mask(data_size, 1);
      std::vector<double> smooth_result = sphere_kde.reconstruction(
          spoke_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result = sphere_kde.reconstruction(
          spoke_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(spoke_data, maskids, reconstruction_mask, smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(spoke_data, maskids, reconstruction_mask, log_smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < data_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], shell_smoothed_spoke[i], 1e-3))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], shell_smoothed_spoke[i], 1e-3))
          ITFAILS;
      }

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
                                std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(
              std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
              std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
        ITFAILS;
    }

    // shell smoothing on spoke array with mask over inner sphere points
    {
      std::vector<double> masked_shell_smoothed_spoke{
          1.0,     2.0,     3.0,     4.0,     5, 6.0,     7.0,     8.0,     9.0,
          1.0,     2.0,     3.0,     4.0,     5, 6.0,     7.0,     8.0,     9.0,
          1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
          1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
          1.0,     2.0,     3.0,     4.0,     5, 6.0,     7.0,     8.0,     9.0,
          1.0,     2.0,     3.0,     4.0,     5, 6.0,     7.0,     8.0,     9.0,
          1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
          1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675};
      std::vector<int> shell_mask{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};

      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          data_size, std::array<double, 3>{1.0e12, 1.0, 0.0});
      const bool dd = false;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 1.0;
      const size_t dim = 2;
      quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);

      std::vector<int> maskids = {1, 2};
      std::vector<double> smooth_result =
          sphere_kde.reconstruction(spoke_data, shell_mask, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result =
          sphere_kde.reconstruction(spoke_data, shell_mask, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(spoke_data, maskids, shell_mask, smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(spoke_data, maskids, shell_mask, log_smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < data_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], masked_shell_smoothed_spoke[i], 1e-3))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], masked_shell_smoothed_spoke[i], 1e-3))
          ITFAILS;
      }

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
                                std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(
              std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
              std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
        ITFAILS;
    }

    // shell smoothing on spoke array with theta reflection
    {
      std::vector<double> shell_smoothed_spoke{
          1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675,
          1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675,
          1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675,
          1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675,
          1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675,
          1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675,
          1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675,
          1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675};
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          data_size, std::array<double, 3>{1.0e12, 1.0, 0.0});
      const bool dd = false;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 1.0;
      const size_t dim = 2;
      quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);

      std::vector<int> maskids = {1};
      std::vector<int> reconstruction_mask(data_size, 1);
      std::vector<double> smooth_result = theta_reflected_sphere_kde.reconstruction(
          spoke_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result = theta_reflected_sphere_kde.reconstruction(
          spoke_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(spoke_data, maskids, reconstruction_mask, smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(spoke_data, maskids, reconstruction_mask, log_smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < data_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], shell_smoothed_spoke[i], 1e-3))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], shell_smoothed_spoke[i], 1e-3))
          ITFAILS;
      }

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
                                std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(
              std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
              std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
        ITFAILS;
    }
  }

  // No mean reconstruction because of small basis functions
  {
    std::vector<double> data{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 0.1, 0.0, 0.0});
    const bool dd = false;
    // two bins per point
    const size_t n_coarse_bins = 5;
    const double max_window_size = 0.1;
    const size_t dim = 1;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < 10; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.1))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.1))
        ITFAILS;
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(
            std::accumulate(data.begin(), data.end(), 0.0),
            std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
      ITFAILS;
  }

  // 2D No mean reconstruction because of small basis functions
  {
    std::vector<double> data{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 0.1, 1.0 / 0.1, 0.0});
    const bool dd = false;
    // two bins per point
    const size_t n_coarse_bins = 5;
    const double max_window_size = 0.1;
    const size_t dim = 2;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < 10; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.1))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.1))
        ITFAILS;
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(
            std::accumulate(data.begin(), data.end(), 0.0),
            std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
      ITFAILS;
  }

  // "Smoothed" reconstruction.
  {
    std::vector<double> data{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 4.0, 0.0, 0.0});
    const bool dd = false;
    // one bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 4.0;
    const size_t dim = 1;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);
    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < 10; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.1, 1e-1))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.1, 1e-1))
        ITFAILS;
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(
            std::accumulate(data.begin(), data.end(), 0.0),
            std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
      ITFAILS;
  }

  // 2D "Smoothed" reconstruction.
  {
    std::vector<double> data{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 4.0, 1.0 / 4.0, 0.0});
    const bool dd = false;
    // one bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 4.0;
    const size_t dim = 2;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);
    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < 10; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.1, 1e-1))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.1, 1e-1))
        ITFAILS;
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(
            std::accumulate(data.begin(), data.end(), 0.0),
            std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
      ITFAILS;
  }

  // No reconstruction because of small basis functions
  {
    std::vector<double> data{0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 0.1, 0.0, 0.0});
    const bool dd = false;
    // 2X bin per point
    const size_t n_coarse_bins = 20;
    const double max_window_size = 0.1;
    const size_t dim = 1;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);
    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < 10; i++) {
      if (!rtt_dsxx::soft_equiv(data[i], smooth_result[i]))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(data[i], log_smooth_result[i]))
        ITFAILS;
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(
            std::accumulate(data.begin(), data.end(), 0.0),
            std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
      ITFAILS;
  }

  // 2D No reconstruction because of small basis in both directions
  {
    std::vector<double> data{0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 0.1, 1.0 / 0.1, 0.0});
    const bool dd = false;
    // 2X bin per point
    const size_t n_coarse_bins = 20;
    const double max_window_size = 0.1;
    const size_t dim = 2;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);
    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < 10; i++) {
      if (!rtt_dsxx::soft_equiv(data[i], smooth_result[i]))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(data[i], log_smooth_result[i]))
        ITFAILS;
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(
            std::accumulate(data.begin(), data.end(), 0.0),
            std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
      ITFAILS;
  }

  // 2D No reconstruction because of small bias in dim=1 keeps dim=2 from
  // accumulating offset data. This test can't be achieved in the opposite
  // direction without a small bandwidth in both dirs because the rows are
  // exactly in line with one another, while the columns are offset.
  {
    std::vector<double> data{0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 0.1, 1.0 / 4.0, 0.0});
    const bool dd = false;
    // 2X bin per point
    const size_t n_coarse_bins = 20;
    const double max_window_size = 4.0;
    const size_t dim = 2;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);
    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < 10; i++) {
      if (!rtt_dsxx::soft_equiv(data[i], smooth_result[i]))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(data[i], log_smooth_result[i]))
        ITFAILS;
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(
            std::accumulate(data.begin(), data.end(), 0.0),
            std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
      ITFAILS;
  }

  // 2D reconstruct only along dim=1 for each row in dim=2
  {
    std::vector<double> data{0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 4.0, 1.0 / 0.1, 0.0});
    const bool dd = false;
    // 2X bin per point
    const size_t n_coarse_bins = 20;
    const double max_window_size = 4.0;
    const size_t dim = 2;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);
    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < 10; i++) {
      if (i < 5) {
        // 0.14 = (0.1*3+0.2*2)/5
        if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.14, 3e-2))
          ITFAILS;
        // 0.14 = (0.1*3+0.2*2)/5
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.14, 3e-2))
          ITFAILS;
      } else {
        // 0.16 = (0.1*2+0.2*3)/5
        if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.16, 3e-2))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.16, 3e-2))
          ITFAILS;
      }
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(
            std::accumulate(data.begin(), data.end(), 0.0),
            std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
      ITFAILS;
  }

  // 2D reconstruct mainly along dim=2 (rows are offset by 0.5 so we have to
  // have a larger bandwidth in dim=1 to get any smoothing in dim=2) for each
  // column in dim=1
  {
    std::vector<double> data{0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 0.6, 1.0 / 4.0, 0.0});
    const bool dd = false;
    // 2X bin per point
    const size_t n_coarse_bins = 20;
    const double max_window_size = 4.0;
    const size_t dim = 2;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);
    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    std::vector<double> bench{0.122267, 0.181788, 0.118212, 0.181788, 0.118212,
                              0.181788, 0.118212, 0.181788, 0.118212, 0.177733};

    std::vector<double> log_bench{0.121638, 0.182268, 0.117873, 0.182268, 0.117873,
                                  0.182268, 0.117873, 0.182268, 0.117873, 0.177799};

    // Check smooth result
    for (int i = 0; i < 10; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], bench[i], 1e-4))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], log_bench[i], 1e-4))
        ITFAILS;
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(
            std::accumulate(data.begin(), data.end(), 0.0),
            std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
      ITFAILS;
  }

  // 2D Smoothed reconstruction should be close to the problem mean of 0.15
  {
    std::vector<double> data{0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 4.0, 1.0 / 4.0, 0.0});
    const bool dd = false;
    // 1 bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 4.0;
    const size_t dim = 2;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);
    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < 10; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.15, 1e-1))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.15, 1e-1))
        ITFAILS;
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
    // Energy conservation
    if (!rtt_dsxx::soft_equiv(
            std::accumulate(data.begin(), data.end(), 0.0),
            std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
      ITFAILS;
  }

  // Smoothed reconstruction should be close to the problem mean of 0.15
  {
    std::vector<double> data{0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 4.0, 0.0, 0.0});
    const bool dd = false;
    // 1 bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 4.0;
    const size_t dim = 1;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);
    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < 10; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.15, 1e-1))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.15, 1e-1))
        ITFAILS;
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(
            std::accumulate(data.begin(), data.end(), 0.0),
            std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
      ITFAILS;
  }

  // No variable band width test.
  {
    std::vector<double> data{0.01, 0.02, 0.1, 0.2, 0.1, 0.02, 0.01, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0, 0.0, 0.0});

    const bool dd = false;
    // 1 bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 1.0;
    const size_t dim = 1;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);

    std::vector<double> bench{0.01446,   0.0172074, 0.10425,  0.172074, 0.131586,
                              0.0172074, 0.040488,  0.172074, 0.131586, 0.15906};

    // Check smooth result
    for (int i = 0; i < 10; i++) {
      if (!rtt_dsxx::soft_equiv(bench[i], smooth_result[i], 1e-4))
        ITFAILS;
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
  }

  // 2D No variable band width test.
  {
    std::vector<double> data{0.01, 0.02, 0.1, 0.2, 0.1, 0.02, 0.01, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0, 1.0 / 4.0, 0.0});

    const bool dd = false;
    // 1 bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 4.0;
    const size_t dim = 2;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);

    std::vector<double> bench{0.0142901, 0.0172733, 0.104099, 0.172733, 0.130699,
                              0.0172733, 0.0396694, 0.172733, 0.130699, 0.160531};

    // Check smooth result
    for (int i = 0; i < 10; i++) {
      if (!rtt_dsxx::soft_equiv(bench[i], smooth_result[i], 1e-4))
        ITFAILS;
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
  }

  // 2D  variable band width test.
  {
    std::vector<double> data{0.01, 0.02, 0.1, 0.2, 0.1, 0.02, 0.01, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0, 1.0 / 4.0, 0.0});

    // lets make the array a little bit more complicated
    one_over_bandwidth_array[9] = {1.0 / 0.5, 1.0 / 4.0, 0.0};
    one_over_bandwidth_array[3] = {1.0 / 1.0, 1.0 / 0.1, 0.0};
    one_over_bandwidth_array[4] = {1.0 / 0.5, 1.0 / 4.0, 0.0};
    one_over_bandwidth_array[2] = {1.0 / 0.1, 1.0 / 4.0, 0.0};
    const bool dd = false;
    // 1 bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 4.0;
    const size_t dim = 2;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);

    std::vector<double> bench{0.0137349, 0.0166022, 0.0972147, 0.194429, 0.0972147,
                              0.0166022, 0.0381284, 0.166022,  0.125622, 0.194429};

    // Check smooth result
    for (int i = 0; i < 10; i++) {
      if (!rtt_dsxx::soft_equiv(bench[i], smooth_result[i], 1e-4))
        ITFAILS;
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
  }

  //  step band width test.
  {
    std::vector<double> data{0.01, 0.02, 0.1, 0.2, 0.1, 0.02, 0.01, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0, 0., 0.0});

    // transition at 1.75
    // lets make the array a little bit more complicated
    one_over_bandwidth_array[0] = {1.0 / 1.75, 0.0, 0.0};
    one_over_bandwidth_array[1] = {1.0 / 0.75, 0.0, 0.0};
    one_over_bandwidth_array[2] = {1.0 / 0.25, 0.0, 0.0};
    one_over_bandwidth_array[3] = {1.0 / 1.25, 0.0, 0.0};
    one_over_bandwidth_array[4] = {1.0 / 2.25, 0.0, 0.0};
    one_over_bandwidth_array[5] = {1.0 / 1.25, 0.0, 0.0};
    one_over_bandwidth_array[6] = {1.0 / 0.25, 0.0, 0.0};
    one_over_bandwidth_array[7] = {1.0 / 0.75, 0.0, 0.0};
    one_over_bandwidth_array[8] = {1.0 / 1.75, 0.0, 0.0};
    one_over_bandwidth_array[9] = {1.0 / 2.75, 0.0, 0.0};
    const bool dd = false;
    // 1 bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 3.0;
    const size_t dim = 1;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);

    std::vector<double> bench{0.01588,   0.0177126, 0.101982, 0.157172, 0.154663,
                              0.0163707, 0.0101982, 0.177126, 0.153908, 0.154988};

    // Check smooth result
    for (int i = 0; i < 10; i++) {
      if (!rtt_dsxx::soft_equiv(bench[i], smooth_result[i], 1e-4))
        ITFAILS;
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
  }

  // what if half of it is negative and the mean is zero
  {
    std::vector<double> data{-0.2, 0.2, -0.2, 0.2, -0.2, 0.2, -0.2, 0.2, -0.2, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 4.0, 0.0, 0.0});

    // lets make the array a little bit more complicated
    const bool dd = false;
    // 1 bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 4.0;
    const size_t dim = 1;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(10, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, reconstruction_mask, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    for (int i = 0; i < 10; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.0, 1e-2))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.0, 1e-2))
        ITFAILS;
    }

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0)))
      ITFAILS;
    // Energy conservation
    if (!rtt_dsxx::soft_equiv(
            std::accumulate(data.begin(), data.end(), 0.0),
            std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0)))
      ITFAILS;
  }

  if (ut.numFails == 0) {
    PASSMSG("KDE checks pass");
  } else {
    FAILMSG("KDE checks failed");
  }
}

//------------------------------------------------------------------------------------------------//
void test_decomposition(ParallelUnitTest &ut) { // NOLINT [hicpp-function-size]
  kde test_kde;

  // test the epan kernel
  double value = test_kde.epan_kernel(0.0);
  if (!rtt_dsxx::soft_equiv(value, 0.75))
    ITFAILS;

  if (rtt_c4::nodes() != 3)
    ITFAILS;

  // spherical reconstruction
  {
    const bool spherical = true;
    const size_t local_size = 24;
    const std::array<double, 3> sphere_center{0.0, -1.0, 0.0};
    kde sphere_kde;
    // reflect on the theta boundary
    kde theta_reflected_sphere_kde({false, false, true, true, false, false});
    const std::array<double, 8> radial_edges{0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0};
    const std::array<double, 9> cosine_edges{-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0};
    const size_t data_size = radial_edges.size() * cosine_edges.size();
    std::vector<std::array<double, 3>> position_array(data_size,
                                                      std::array<double, 3>{0.0, 0.0, 0.0});

    std::vector<double> shell_data(data_size, 0.0);
    std::vector<double> spoke_data(data_size, 0.0);
    std::vector<double> shell_weights(data_size, 0.0);
    size_t point_i = 0;
    size_t ri = 0;
    for (auto &r : radial_edges) {
      size_t mui = 0;
      for (auto &mu : cosine_edges) {
        spoke_data[point_i] = static_cast<double>(mui) + 1.0;
        shell_data[point_i] = static_cast<double>(ri) + 1.0;
        if (ri < 4)
          shell_weights[point_i] = 1e3;
        else
          shell_weights[point_i] = 1.0;
        double rel_y = r * mu;
        position_array[point_i][0] =
            rtt_dsxx::soft_equiv(r * r, rel_y * rel_y, 1e-6) ? 0.0 : sqrt(r * r - rel_y * rel_y);
        position_array[point_i][1] = sphere_center[1] + rel_y;
        point_i++;
        mui++;
      }
      ri++;
    }
    std::vector<double> spoke_smoothed_shells{
        3.7717,  3.7717,  3.7717,  3.7717,  3.7717,  3.7717,  3.7717,  3.7717,  3.7717,
        3.83452, 3.83452, 3.83452, 3.83452, 3.83452, 3.83452, 3.83452, 3.83452, 3.83452,
        3.89388, 3.89388, 3.89388, 3.89388, 3.89388, 3.89388, 3.89388, 3.89388, 3.89388,
        3.95025, 3.95025, 3.95025, 3.95025, 3.95025, 3.95025, 3.95025, 3.95025, 3.95025,
        4.24469, 4.24469, 4.24469, 4.24469, 4.24469, 4.24469, 4.24469, 4.24469, 4.24469,
        4.67528, 4.67528, 4.67528, 4.67528, 4.67528, 4.67528, 4.67528, 4.67528, 4.67528,
        5.21495, 5.21495, 5.21495, 5.21495, 5.21495, 5.21495, 5.21495, 5.21495, 5.21495,
        6.41471, 6.41472, 6.41472, 6.41472, 6.41471, 6.41471, 6.41471, 6.41471, 6.41471};
    std::vector<double> shell_smoothed_spoke{
        1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
        1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
        1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
        1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
        1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
        1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
        1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
        1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675};
    std::vector<double> masked_shell_smoothed_spoke{
        1.0,     2.0,     3.0,     4.0,     5, 6.0,     7.0,     8.0,     9.0,
        1.0,     2.0,     3.0,     4.0,     5, 6.0,     7.0,     8.0,     9.0,
        1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
        1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
        1.0,     2.0,     3.0,     4.0,     5, 6.0,     7.0,     8.0,     9.0,
        1.0,     2.0,     3.0,     4.0,     5, 6.0,     7.0,     8.0,     9.0,
        1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675,
        1.32325, 2.77318, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.22682, 8.67675};
    std::vector<int> shell_mask{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    std::vector<double> shell_smoothed_spoke_theta_ref{
        1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675,
        1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675,
        1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675,
        1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675,
        1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675,
        1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675,
        1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675,
        1.32325, 2.5489, 3.67834, 4.30104, 5, 5.69896, 6.32166, 7.4511, 8.67675};
    std::vector<double> spoke_weighted_smoothed_shells{
        2.49765, 2.49765, 2.49765, 2.49765, 2.49765, 2.49765, 2.49765, 2.49765, 2.49765,
        2.49922, 2.49922, 2.49922, 2.49922, 2.49922, 2.49922, 2.49922, 2.49922, 2.49922,
        2.50078, 2.50078, 2.50079, 2.50078, 2.50078, 2.50078, 2.50078, 2.50078, 2.50078,
        2.50235, 2.50235, 2.50235, 2.50236, 2.50235, 2.50235, 2.50235, 2.50235, 2.50235,
        6.2,     6.2,     6.2,     6.2,     6.2,     6.2,     6.2,     6.2,     6.2,
        6.41379, 6.41379, 6.41379, 6.41379, 6.41379, 6.41379, 6.41379, 6.41379, 6.41379,
        6.58621, 6.58621, 6.58621, 6.58621, 6.58621, 6.58621, 6.58621, 6.58621, 6.58621,
        6.8,     6.8,     6.8,     6.8,     6.8,     6.8,     6.8,     6.8,     6.8};

    std::vector<double> dd_const_data(local_size, 1.0);
    std::vector<double> dd_spoke_data(local_size);
    std::vector<double> dd_shell_data(local_size);
    std::vector<double> dd_shell_weights(local_size);
    std::vector<double> dd_spoke_smoothed_shells(local_size);
    std::vector<double> dd_spoke_weighted_smoothed_shells(local_size);
    std::vector<double> dd_shell_smoothed_spoke(local_size);
    std::vector<double> dd_masked_shell_smoothed_spoke(local_size);
    std::vector<int> dd_shell_mask(local_size);
    std::vector<double> dd_shell_smoothed_spoke_theta_ref(local_size);
    std::vector<std::array<double, 3>> dd_position_array(local_size);
    for (size_t i = 0; i < local_size; i++) {
      dd_spoke_data[i] = spoke_data[i + rtt_c4::node() * local_size];
      dd_shell_data[i] = shell_data[i + rtt_c4::node() * local_size];
      dd_shell_weights[i] = shell_weights[i + rtt_c4::node() * local_size];
      dd_spoke_smoothed_shells[i] = spoke_smoothed_shells[i + rtt_c4::node() * local_size];
      dd_spoke_weighted_smoothed_shells[i] =
          spoke_weighted_smoothed_shells[i + rtt_c4::node() * local_size];
      dd_shell_smoothed_spoke[i] = shell_smoothed_spoke[i + rtt_c4::node() * local_size];
      dd_masked_shell_smoothed_spoke[i] =
          masked_shell_smoothed_spoke[i + rtt_c4::node() * local_size];
      dd_shell_mask[i] = shell_mask[i + rtt_c4::node() * local_size];
      dd_shell_smoothed_spoke_theta_ref[i] =
          shell_smoothed_spoke_theta_ref[i + rtt_c4::node() * local_size];
      dd_position_array[i] = position_array[i + rtt_c4::node() * local_size];
    }

    // zero reconstruction array on a shell mesh
    {
      std::vector<double> zero_data(local_size, 0.0);
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          local_size, std::array<double, 3>{1.0, 1.0e12, 0.0});
      const bool dd = true;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 2.0;
      const size_t dim = 2;
      quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);
      std::vector<int> maskids = {1};
      std::vector<int> reconstruction_mask(local_size, 1);
      std::vector<double> smooth_result = sphere_kde.reconstruction(
          zero_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result = sphere_kde.log_reconstruction(
          zero_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      std::vector<double> sampled_smooth_result = sphere_kde.sampled_reconstruction(
          zero_data, reconstruction_mask, one_over_bandwidth_array, qindex);

      // Apply Conservation
      sphere_kde.apply_conservation(zero_data, maskids, reconstruction_mask, smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(zero_data, maskids, reconstruction_mask, log_smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(zero_data, maskids, reconstruction_mask, sampled_smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < local_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], zero_data[i]))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], zero_data[i]))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(sampled_smooth_result[i], zero_data[i]))
          ITFAILS;
      }

      double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
      rtt_c4::global_sum(smooth_conservation);
      double log_smooth_conservation =
          std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
      rtt_c4::global_sum(log_smooth_conservation);
      double sampled_smooth_conservation =
          std::accumulate(sampled_smooth_result.begin(), sampled_smooth_result.end(), 0.0);
      rtt_c4::global_sum(sampled_smooth_conservation);

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(0.0, smooth_conservation))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(0.0, log_smooth_conservation))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(0.0, sampled_smooth_conservation))
        ITFAILS;
    }

    // spoke reconstruction array
    {
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          local_size, std::array<double, 3>{1.0, 1.0e12, 0.0});
      const bool dd = true;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 2.0;
      const size_t dim = 2;
      quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);

      std::vector<int> maskids = {1};
      std::vector<int> reconstruction_mask(local_size, 1);
      std::vector<double> smooth_result = sphere_kde.reconstruction(
          dd_spoke_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result = sphere_kde.log_reconstruction(
          dd_spoke_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(dd_spoke_data, maskids, reconstruction_mask, smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(dd_spoke_data, maskids, reconstruction_mask, log_smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < local_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], dd_spoke_data[i]))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], dd_spoke_data[i]))
          ITFAILS;
      }

      double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
      rtt_c4::global_sum(smooth_conservation);
      double log_smooth_conservation =
          std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
      rtt_c4::global_sum(log_smooth_conservation);

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
                                smooth_conservation))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
                                log_smooth_conservation))
        ITFAILS;
    }

    // shell reconstruction array
    {
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          local_size, std::array<double, 3>{1.0e12, 1.0, 0.0});
      const bool dd = true;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 2.0;
      const size_t dim = 2;
      quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);

      std::vector<int> maskids = {1};
      std::vector<int> reconstruction_mask(local_size, 1);
      std::vector<double> smooth_result = sphere_kde.reconstruction(
          dd_shell_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result = sphere_kde.log_reconstruction(
          dd_shell_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(dd_shell_data, maskids, reconstruction_mask, smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(dd_shell_data, maskids, reconstruction_mask, log_smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < local_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], dd_shell_data[i]))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], dd_shell_data[i]))
          ITFAILS;
      }

      double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
      rtt_c4::global_sum(smooth_conservation);
      double log_smooth_conservation =
          std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
      rtt_c4::global_sum(log_smooth_conservation);

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(std::accumulate(shell_data.begin(), shell_data.end(), 0.0),
                                smooth_conservation))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(std::accumulate(shell_data.begin(), shell_data.end(), 0.0),
                                log_smooth_conservation))
        ITFAILS;
    }

    // spoke smoothing on shell array
    {
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          local_size, std::array<double, 3>{1.0, 1.0e12, 0.0});
      const bool dd = true;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 2.0;
      const size_t dim = 2;
      quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);

      std::vector<int> maskids = {1};
      std::vector<int> reconstruction_mask(local_size, 1);
      std::vector<double> smooth_result = sphere_kde.reconstruction(
          dd_shell_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result = sphere_kde.reconstruction(
          dd_shell_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(dd_shell_data, maskids, reconstruction_mask, smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(dd_shell_data, maskids, reconstruction_mask, log_smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < local_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], dd_spoke_smoothed_shells[i], 1e-3))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], dd_spoke_smoothed_shells[i], 1e-3))
          ITFAILS;
      }

      double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
      rtt_c4::global_sum(smooth_conservation);
      double log_smooth_conservation =
          std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
      rtt_c4::global_sum(log_smooth_conservation);

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(std::accumulate(shell_data.begin(), shell_data.end(), 0.0),
                                smooth_conservation))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(std::accumulate(shell_data.begin(), shell_data.end(), 0.0),
                                log_smooth_conservation))
        ITFAILS;
    }

    // spoke smoothing on shell array with shell weights
    {
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          local_size, std::array<double, 3>{1.0, 1.0e12, 0.0});
      const bool dd = true;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 2.0;
      const size_t dim = 2;
      quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);

      std::vector<int> maskids = {1};
      std::vector<int> reconstruction_mask(local_size, 1);
      std::vector<double> smooth_result = sphere_kde.weighted_reconstruction(
          dd_shell_data, dd_shell_weights, reconstruction_mask, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(dd_shell_data, maskids, reconstruction_mask, smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < local_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], dd_spoke_weighted_smoothed_shells[i], 1e-3))
          ITFAILS;
      }

      double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
      rtt_c4::global_sum(smooth_conservation);

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(std::accumulate(shell_data.begin(), shell_data.end(), 0.0),
                                smooth_conservation))
        ITFAILS;
    }

    // shell smoothing on spoke array
    {
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          local_size, std::array<double, 3>{1.0e12, 1.0, 0.0});
      const bool dd = true;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 1.0;
      const size_t dim = 2;
      quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);

      std::vector<int> maskids = {1};
      std::vector<int> reconstruction_mask(local_size, 1);
      std::vector<double> smooth_result = sphere_kde.reconstruction(
          dd_spoke_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result = sphere_kde.reconstruction(
          dd_spoke_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(dd_spoke_data, maskids, reconstruction_mask, smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(dd_spoke_data, maskids, reconstruction_mask, log_smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < local_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], dd_shell_smoothed_spoke[i], 1e-3))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], dd_shell_smoothed_spoke[i], 1e-3))
          ITFAILS;
      }

      double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
      rtt_c4::global_sum(smooth_conservation);
      double log_smooth_conservation =
          std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
      rtt_c4::global_sum(log_smooth_conservation);

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
                                smooth_conservation))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
                                log_smooth_conservation))
        ITFAILS;
    }

    // shell smoothing on spoke array with mask
    {
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          local_size, std::array<double, 3>{1.0e12, 1.0, 0.0});
      const bool dd = true;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 1.0;
      const size_t dim = 2;
      quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);

      std::vector<double> smooth_result =
          sphere_kde.reconstruction(dd_spoke_data, dd_shell_mask, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result =
          sphere_kde.reconstruction(dd_spoke_data, dd_shell_mask, one_over_bandwidth_array, qindex);
      // Apply Conservation
      std::vector<int> maskids = {2, 1};
      sphere_kde.apply_conservation(dd_spoke_data, maskids, dd_shell_mask, smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(dd_spoke_data, maskids, dd_shell_mask, log_smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < local_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], dd_masked_shell_smoothed_spoke[i], 1e-3))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], dd_masked_shell_smoothed_spoke[i], 1e-3))
          ITFAILS;
      }

      double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
      rtt_c4::global_sum(smooth_conservation);
      double log_smooth_conservation =
          std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
      rtt_c4::global_sum(log_smooth_conservation);

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
                                smooth_conservation))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
                                log_smooth_conservation))
        ITFAILS;
    }

    // shell smoothing on spoke array with theta reflection
    {
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          local_size, std::array<double, 3>{1.0e12, 1.0, 0.0});
      const bool dd = true;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 1.0;
      const size_t dim = 2;
      quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd, spherical,
                         sphere_center);

      std::vector<int> maskids = {1};
      std::vector<int> reconstruction_mask(local_size, 1);
      std::vector<double> smooth_result = theta_reflected_sphere_kde.reconstruction(
          dd_spoke_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result = theta_reflected_sphere_kde.reconstruction(
          dd_spoke_data, reconstruction_mask, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(dd_spoke_data, maskids, reconstruction_mask, smooth_result,
                                    qindex.domain_decomposed);
      sphere_kde.apply_conservation(dd_spoke_data, maskids, reconstruction_mask, log_smooth_result,
                                    qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < local_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], dd_shell_smoothed_spoke_theta_ref[i], 1e-3))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], dd_shell_smoothed_spoke_theta_ref[i], 1e-3))
          ITFAILS;
      }

      double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
      rtt_c4::global_sum(smooth_conservation);
      double log_smooth_conservation =
          std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
      rtt_c4::global_sum(log_smooth_conservation);

      // Energy conservation
      if (!rtt_dsxx::soft_equiv(std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
                                smooth_conservation))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(std::accumulate(spoke_data.begin(), spoke_data.end(), 0.0),
                                log_smooth_conservation))
        ITFAILS;
    }
  }

  int local_size = 3;
  // give the odd size to the final rank to make striding easy
  if (rtt_c4::node() == 2)
    local_size = 4;

  // No mean reconstruction because of small basis functions
  {
    std::vector<double> data{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 0.1, 0.0, 0.0});

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 1 bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 0.1;
    const size_t dim = 1;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result = test_kde.log_reconstruction(
        dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.1))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.1))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);
    double log_smooth_conservation =
        std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
    rtt_c4::global_sum(log_smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              log_smooth_conservation))
      ITFAILS;
  }

  // 2D No mean reconstruction because of small basis functions
  {
    std::vector<double> data{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 0.1, 1.0 / 0.1, 0.0});

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 1 bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 0.1;
    const size_t dim = 2;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result = test_kde.log_reconstruction(
        dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.1))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.1))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);
    double log_smooth_conservation =
        std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
    rtt_c4::global_sum(log_smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              log_smooth_conservation))
      ITFAILS;
  }

  // "Smoothed" reconstruction.
  {
    std::vector<double> data{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 4.0, 0.0, 0.0});

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 1/2 bin per point
    const size_t n_coarse_bins = 5;
    const double max_window_size = 4.0;
    const size_t dim = 1;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result = test_kde.log_reconstruction(
        dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.1))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.1))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);
    double log_smooth_conservation =
        std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
    rtt_c4::global_sum(log_smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              log_smooth_conservation))
      ITFAILS;
  }

  // 2D "Smoothed" reconstruction.
  {
    std::vector<double> data{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 4.0, 1.0 / 4.0, 0.0});

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 1/2 bin per point
    const size_t n_coarse_bins = 5;
    const double max_window_size = 4.0;
    const size_t dim = 1;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result = test_kde.log_reconstruction(
        dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.1))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.1))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);
    double log_smooth_conservation =
        std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
    rtt_c4::global_sum(log_smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              log_smooth_conservation))
      ITFAILS;
  }

  // No reconstruction because of small basis functions
  {
    std::vector<double> data{0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 0.1, 0.0, 0.0});

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 2x bin per point
    const size_t n_coarse_bins = 20;
    const double max_window_size = 0.1;
    const size_t dim = 1;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result = test_kde.log_reconstruction(
        dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(data[i + rtt_c4::node() * 3], smooth_result[i]))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(data[i + rtt_c4::node() * 3], log_smooth_result[i]))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);
    double log_smooth_conservation =
        std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
    rtt_c4::global_sum(log_smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              log_smooth_conservation))
      ITFAILS;
  }

  // 2D No reconstruction because of small basis functions
  {
    std::vector<double> data{0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 0.1, 1.0, 0.0});

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 2x bin per point
    const size_t n_coarse_bins = 20;
    const double max_window_size = 1.0;
    const size_t dim = 2;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result = test_kde.log_reconstruction(
        dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(data[i + rtt_c4::node() * 3], smooth_result[i]))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(data[i + rtt_c4::node() * 3], log_smooth_result[i]))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);
    double log_smooth_conservation =
        std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
    rtt_c4::global_sum(log_smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              log_smooth_conservation))
      ITFAILS;
  }

  // 2D No reconstruction because of small bias in dim=1 keeps dim=2 from
  // accumulating offset data. This test can't be achieved in the opposite
  // direction without a small bandwidth in both dirs because the rows are
  // exactly in line with one another, while the columns are offset.
  {
    std::vector<double> data{0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 0.1, 1.0 / 4.0, 0.0});

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 2x bin per point
    const size_t n_coarse_bins = 20;
    const double max_window_size = 4.0;
    const size_t dim = 2;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result = test_kde.log_reconstruction(
        dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(data[i + rtt_c4::node() * 3], smooth_result[i]))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(data[i + rtt_c4::node() * 3], log_smooth_result[i]))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);
    double log_smooth_conservation =
        std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
    rtt_c4::global_sum(log_smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              log_smooth_conservation))
      ITFAILS;
  }

  // 2D reconstruct only along dim=1 for each row in dim=2
  {
    std::vector<double> data{0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 4.0, 1.0 / 0.1, 0.0});

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 2x bin per point
    const size_t n_coarse_bins = 20;
    const double max_window_size = 4.0;
    const size_t dim = 2;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result = test_kde.log_reconstruction(
        dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    // Check smooth result
    for (int i = 0; i < local_size; i++) {
      if (dd_position_array[i][1] > 0.0) {
        // 0.14 = (0.1*3+0.2*2)/5
        if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.14, 3e-2))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.14, 3e-2))
          ITFAILS;
      } else {
        // 0.16 = (0.1*2+0.2*3)/5
        if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.16, 3e-2))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.16, 3e-2))
          ITFAILS;
      }
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);
    double log_smooth_conservation =
        std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
    rtt_c4::global_sum(log_smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              log_smooth_conservation))
      ITFAILS;
  }

  // 2D reconstruct mainly along dim=2 (rows are offset by 0.5 so we have to
  // have a larger bandwidth in dim=1 to get any smoothing in dim=2) for each
  // column in dim=1
  {
    std::vector<double> data{0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 0.6, 1.0 / 4.0, 0.0});

    std::vector<double> bench{0.122267, 0.181788, 0.118212, 0.181788, 0.118212,
                              0.181788, 0.118212, 0.181788, 0.118212, 0.177733};
    std::vector<double> log_bench{0.121638, 0.182268, 0.117873, 0.182268, 0.117873,
                                  0.182268, 0.117873, 0.182268, 0.117873, 0.177799};

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});
    std::vector<double> dd_bench(local_size, 0.0);
    std::vector<double> log_dd_bench(local_size, 0.0);

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
      dd_bench[i] = bench[i + rtt_c4::node() * 3];
      log_dd_bench[i] = log_bench[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 2x bin per point
    const size_t n_coarse_bins = 20;
    const double max_window_size = 4.0;
    const size_t dim = 2;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result = test_kde.log_reconstruction(
        dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], dd_bench[i], 1e-4))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], log_dd_bench[i], 1e-4))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);
    double log_smooth_conservation =
        std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
    rtt_c4::global_sum(log_smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              log_smooth_conservation))
      ITFAILS;
  }

  // Smoothed reconstruction should be close to the problem mean of 0.15
  {
    std::vector<double> data{0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 4.0, 0.0, 0.0});

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 1x bin per point
    const size_t n_coarse_bins = 10;
    // window size must be 2X bigger then biggest bandwidth
    const double max_window_size = 9.0;
    const size_t dim = 1;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result = test_kde.log_reconstruction(
        dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.15, 1e-2))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.15, 1e-2))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);
    double log_smooth_conservation =
        std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
    rtt_c4::global_sum(log_smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              log_smooth_conservation))
      ITFAILS;
  }

  // 2D Smoothed reconstruction should be close to the problem mean of 0.15
  {
    std::vector<double> data{0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 4.0, 1.0 / 4.0, 0.0});

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 1x bin per point
    const size_t n_coarse_bins = 10;
    // window size must be 2X bigger then biggest bandwidth
    const double max_window_size = 9.0;
    const size_t dim = 2;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result = test_kde.log_reconstruction(
        dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.15, 1e-2))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.15, 1e-2))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);
    double log_smooth_conservation =
        std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(log_smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              log_smooth_conservation))
      ITFAILS;
  }

  // No  variable band width test.
  {
    std::vector<double> data{0.01, 0.02, 0.1, 0.2, 0.1, 0.02, 0.01, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0, 0.0, 0.0});
    std::vector<double> bench{0.01446,   0.0172074, 0.10425,  0.172074, 0.131586,
                              0.0172074, 0.040488,  0.172074, 0.131586, 0.15906};

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 1x bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 1.0;
    const size_t dim = 1;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(bench[i + rtt_c4::node() * 3], smooth_result[i], 1e-4))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
  }

  // 2D no  variable band width test.
  {
    std::vector<double> data{0.01, 0.02, 0.1, 0.2, 0.1, 0.02, 0.01, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0, 1.0 / 4.0, 0.0});
    std::vector<double> bench{0.0142901, 0.0172733, 0.104099, 0.172733, 0.130699,
                              0.0172733, 0.0396694, 0.172733, 0.130699, 0.160531};

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 1x bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 4.0;
    const size_t dim = 2;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(bench[i + rtt_c4::node() * 3], smooth_result[i], 1e-4))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
  }

  //  variable band width test.
  {
    std::vector<double> data{0.01, 0.02, 0.1, 0.2, 0.1, 0.02, 0.01, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0, 0.0, 0.0});

    // lets make the array a little bit more complicated
    one_over_bandwidth_array[9] = {1.0 / 0.5, 0.0, 0.0};
    one_over_bandwidth_array[3] = {1.0 / 0.1, 0.0, 0.0};
    one_over_bandwidth_array[4] = {1.0 / 0.5, 0.0, 0.0};
    one_over_bandwidth_array[2] = {1.0 / 2.0, 0.0, 0.0};

    std::vector<double> bench{0.0139053, 0.0165473, 0.0953673, 0.194674, 0.0973372,
                              0.0165473, 0.0389349, 0.165473,  0.126538, 0.194674};

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 1x bin per point
    const size_t n_coarse_bins = 10;
    // max window size must be 2x the max bandwidth size
    const double max_window_size = 4.0;
    const size_t dim = 1;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(bench[i + rtt_c4::node() * 3], smooth_result[i], 1e-4))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
  }

  // 2D variable band width test.
  {
    std::vector<double> data{0.01, 0.02, 0.1, 0.2, 0.1, 0.02, 0.01, 0.2, 0.1, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0, 1.0 / 4.0, 0.0});

    // lets make the array a little bit more complicated
    one_over_bandwidth_array[9] = {1.0 / 0.5, 1.0 / 4.0, 0.0};
    one_over_bandwidth_array[3] = {1.0 / 1.0, 1.0 / 0.1, 0.0};
    one_over_bandwidth_array[4] = {1.0 / 0.5, 1.0 / 4.0, 0.0};
    one_over_bandwidth_array[2] = {1.0 / 0.1, 1.0 / 4.0, 0.0};

    std::vector<double> bench{0.0137349, 0.0166022, 0.0972147, 0.194429, 0.0972147,
                              0.0166022, 0.0381284, 0.166022,  0.125622, 0.194429};

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 1x bin per point
    const size_t n_coarse_bins = 10;
    // max window size must be 2x the max bandwidth size
    const double max_window_size = 4.0;
    const size_t dim = 2;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);

    // Check smooth result
    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(bench[i + rtt_c4::node() * 3], smooth_result[i], 1e-4))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
  }

  // what if half of it is negative and the mean is zero for a reconstruction
  {
    std::vector<double> data{-0.2, 0.2, -0.2, 0.2, -0.2, 0.2, -0.2, 0.2, -0.2, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 4.0, 0.0, 0.0});

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 1x bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 4.0;
    const size_t dim = 1;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result = test_kde.log_reconstruction(
        dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.0, 1e-2))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.0, 1e-2))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);
    double log_smooth_conservation =
        std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
    rtt_c4::global_sum(log_smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              log_smooth_conservation))
      ITFAILS;
  }

  // what if half of it is negative and the mean is zero for a reconstruction
  // what if we also reflect the bc
  {
    kde refl_kde({true, true, true, true, true, true});
    std::vector<double> data{-0.2, 0.2, -0.2, 0.2, -0.2, 0.2, -0.2, 0.2, -0.2, 0.2};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    std::vector<std::array<double, 3>> one_over_bandwidth_array(
        10, std::array<double, 3>{1.0 / 4.0, 1.0 / 4.0, 0.0});

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});
    std::vector<std::array<double, 3>> dd_one_over_bandwidth_array(
        local_size, std::array<double, 3>{0.0, 0., 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
      dd_one_over_bandwidth_array[i] = one_over_bandwidth_array[i + rtt_c4::node() * 3];
    }

    const bool dd = true;
    // 1x bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 4.0;
    const size_t dim = 2;
    quick_index qindex(dim, dd_position_array, max_window_size, n_coarse_bins, dd);

    std::vector<int> maskids = {1};
    std::vector<int> reconstruction_mask(local_size, 1);
    std::vector<double> smooth_result =
        refl_kde.reconstruction(dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result = refl_kde.log_reconstruction(
        dd_data, reconstruction_mask, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    refl_kde.apply_conservation(dd_data, maskids, reconstruction_mask, smooth_result,
                                qindex.domain_decomposed);
    refl_kde.apply_conservation(dd_data, maskids, reconstruction_mask, log_smooth_result,
                                qindex.domain_decomposed);

    for (int i = 0; i < local_size; i++) {
      if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.0, 1e-2))
        ITFAILS;
      if (!rtt_dsxx::soft_equiv(log_smooth_result[i], 0.0, 1e-2))
        ITFAILS;
    }

    double smooth_conservation = std::accumulate(smooth_result.begin(), smooth_result.end(), 0.0);
    rtt_c4::global_sum(smooth_conservation);
    double log_smooth_conservation =
        std::accumulate(log_smooth_result.begin(), log_smooth_result.end(), 0.0);
    rtt_c4::global_sum(log_smooth_conservation);

    // Energy conservation
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0), smooth_conservation))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(std::accumulate(data.begin(), data.end(), 0.0),
                              log_smooth_conservation))
      ITFAILS;
  }

  if (ut.numFails == 0) {
    PASSMSG("KDE DD checks pass");
  } else {
    FAILMSG("KDE DD checks failed");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ParallelUnitTest ut(argc, argv, release);
  try {
    // >>> UNIT TESTS
    test_replication(ut);
    if (nodes() == 3)
      test_decomposition(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstkde.cc
//------------------------------------------------------------------------------------------------//
