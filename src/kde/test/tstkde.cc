//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/test/tstkde.cc
 * \author Mathew Cleveland
 * \date   Nov. 10th 2020
 * \brief  KDE function tests
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
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
void test_replication(ParallelUnitTest &ut) {
  kde test_kde;

  // test the epan kernel
  double value = test_kde.epan_kernel(0.0);
  if (!rtt_dsxx::soft_equiv(value, 0.75))
    ITFAILS;

  // test some public sphere calculations
  {
    const std::array<double, 3> sphere_center{0.0, 0.0, 0.0};
    const std::array<double, 3> location{sqrt(2), sqrt(2), 0.0};
    const std::array<double, 3> location2{0.0, 2.0, 0.0};
    const double radius = 2.0;
    const double small_radius = 1.0;
    const double pi_over_4 = 0.78539816;
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(test_kde.calc_radius(sphere_center, location), 2.0));
    FAIL_IF_NOT(
        rtt_dsxx::soft_equiv(test_kde.calc_arch_length(sphere_center, radius, location, location2),
                             2.0 * pi_over_4, 1e-6));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(
        test_kde.calc_arch_length(sphere_center, small_radius, location, location2), pi_over_4,
        1e-6));
  }

  // spherical reconstruction
  {
    const std::array<double, 3> sphere_center{0.0, -1.0, 0.0};
    const double max_radius = 1.0;
    const double min_radius = 0.0;
    kde sphere_kde;
    sphere_kde.set_sphere_center(sphere_center, min_radius, max_radius);
    const std::array<double, 8> radial_edges{0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0};
    const std::array<double, 9> cosine_edges{-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0};
    const size_t data_size = radial_edges.size() * cosine_edges.size();
    std::vector<std::array<double, 3>> position_array(data_size,
                                                      std::array<double, 3>{0.0, 0.0, 0.0});

    const std::array<double, 8> smoothed_shell_ref{3.7717,  3.83452, 3.89388, 3.95025,
                                                   4.24469, 4.67528, 5.21495, 6.41471};
    std::vector<double> shell_data(data_size, 0.0);
    std::vector<double> spoke_data(data_size, 0.0);
    std::vector<double> const_data(data_size, 1.0);
    std::vector<double> spoke_smoothed_shells(data_size, 1.0);
    size_t point_i = 0;
    size_t ri = 0;
    for (auto &r : radial_edges) {
      size_t mui = 0;
      for (auto &mu : cosine_edges) {
        spoke_data[point_i] = static_cast<double>(mui) + 1.0;
        shell_data[point_i] = static_cast<double>(ri) + 1.0;
        spoke_smoothed_shells[point_i] = smoothed_shell_ref[ri];
        double rel_y = r * mu;
        position_array[point_i][0] =
            rtt_dsxx::soft_equiv(r, rel_y, 1e-6) ? 0.0 : sqrt(r * r - rel_y * rel_y);
        position_array[point_i][1] = sphere_center[1] + rel_y;
        point_i++;
        mui++;
      }
      ri++;
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
      quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

      std::vector<double> smooth_result =
          sphere_kde.reconstruction(spoke_data, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result =
          sphere_kde.log_reconstruction(spoke_data, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(spoke_data, smooth_result, qindex.domain_decomposed);
      sphere_kde.apply_conservation(spoke_data, log_smooth_result, qindex.domain_decomposed);

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
      quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

      std::vector<double> smooth_result =
          sphere_kde.reconstruction(shell_data, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result =
          sphere_kde.log_reconstruction(shell_data, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(shell_data, smooth_result, qindex.domain_decomposed);
      sphere_kde.apply_conservation(shell_data, log_smooth_result, qindex.domain_decomposed);

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
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          data_size, std::array<double, 3>{1.0, 1.0e12, 0.0});
      const bool dd = false;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 1.0;
      const size_t dim = 2;
      quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

      std::vector<double> smooth_result =
          sphere_kde.reconstruction(shell_data, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result =
          sphere_kde.reconstruction(shell_data, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(shell_data, smooth_result, qindex.domain_decomposed);
      sphere_kde.apply_conservation(shell_data, log_smooth_result, qindex.domain_decomposed);

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

    // shell smoothing on spoke array
    {
      std::vector<double> shell_smooth_spoke{
          4.83303, 4.83501, 4.83589, 4.83663, 4.83732, 4.838,   4.83874, 4.83963, 4.8416,
          4.82008, 4.82805, 4.8316,  4.83456, 4.83732, 4.84007, 4.84303, 4.84658, 4.85456,
          4.79812, 4.81637, 4.82441, 4.8311,  4.83732, 4.84353, 4.85022, 4.85827, 4.87651,
          4.76662, 4.79982, 4.81427, 4.82622, 4.83732, 4.84841, 4.86036, 4.87481, 4.90802,
          4.0521,  4.2347,  4.68473, 4.7647,  5,       4.98777, 5.15372, 6.36747, 5.62253,
          2.63366, 3.76559, 4.12929, 4.48833, 5,       6.92886, 6.67554, 7.37286, 7.04097,
          1.67316, 3.02095, 3.73636, 4.4428,  5,       8.70717, 8.01969, 6.65368, 8.00148,
          1.2802,  2.68295, 3.55866, 4.16109, 5,       5.51354, 6.11598, 6.99168, 8.39443};
      std::vector<std::array<double, 3>> one_over_bandwidth_array(
          data_size, std::array<double, 3>{1.0e12, 1.0, 0.0});
      const bool dd = false;
      // two bins per point
      const size_t n_coarse_bins = 5;
      const double max_window_size = 1.0;
      const size_t dim = 2;
      quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

      std::vector<double> smooth_result =
          sphere_kde.reconstruction(spoke_data, one_over_bandwidth_array, qindex);
      std::vector<double> log_smooth_result =
          sphere_kde.reconstruction(spoke_data, one_over_bandwidth_array, qindex);
      // Apply Conservation
      sphere_kde.apply_conservation(spoke_data, smooth_result, qindex.domain_decomposed);
      sphere_kde.apply_conservation(spoke_data, log_smooth_result, qindex.domain_decomposed);

      // Check smooth result
      for (size_t i = 0; i < data_size; i++) {
        if (!rtt_dsxx::soft_equiv(smooth_result[i], shell_smooth_spoke[i], 1e-3))
          ITFAILS;
        if (!rtt_dsxx::soft_equiv(log_smooth_result[i], shell_smooth_spoke[i], 1e-3))
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
        10, std::array<double, 3>{1.0 / 0.1, 0., 0.0});
    const bool dd = false;
    // two bins per point
    const size_t n_coarse_bins = 5;
    const double max_window_size = 0.1;
    const size_t dim = 1;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(data, log_smooth_result, qindex.domain_decomposed);

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(data, log_smooth_result, qindex.domain_decomposed);

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
        10, std::array<double, 3>{1.0 / 4.0, 0., 0.0});
    const bool dd = false;
    // one bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 4.0;
    const size_t dim = 1;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(data, log_smooth_result, qindex.domain_decomposed);

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
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(data, log_smooth_result, qindex.domain_decomposed);

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
        10, std::array<double, 3>{1.0 / 0.1, 0., 0.0});
    const bool dd = false;
    // 2X bin per point
    const size_t n_coarse_bins = 20;
    const double max_window_size = 0.1;
    const size_t dim = 1;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(data, log_smooth_result, qindex.domain_decomposed);

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
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(data, log_smooth_result, qindex.domain_decomposed);

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
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(data, log_smooth_result, qindex.domain_decomposed);

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
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(data, log_smooth_result, qindex.domain_decomposed);

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
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(data, log_smooth_result, qindex.domain_decomposed);

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
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(data, log_smooth_result, qindex.domain_decomposed);

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
        10, std::array<double, 3>{1.0 / 4.0, 0., 0.0});
    const bool dd = false;
    // 1 bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 4.0;
    const size_t dim = 1;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);
    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(data, log_smooth_result, qindex.domain_decomposed);

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
        10, std::array<double, 3>{1.0, 0., 0.0});

    const bool dd = false;
    // 1 bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 1.0;
    const size_t dim = 1;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);

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
    one_over_bandwidth_array[9] = {1.0 / 0.5, 1.0 / 4.0, 0.};
    one_over_bandwidth_array[3] = {1.0 / 1.0, 1.0 / 0.1, 0.};
    one_over_bandwidth_array[4] = {1.0 / 0.5, 1.0 / 4.0, 0.};
    one_over_bandwidth_array[2] = {1.0 / 0.1, 1.0 / 4.0, 0.};
    const bool dd = false;
    // 1 bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 4.0;
    const size_t dim = 2;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);

    std::vector<double> bench{0.0131256, 0.0158657, 0.1,      0.2,      0.1,
                              0.0158657, 0.0364369, 0.158657, 0.120049, 0.2};

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
    one_over_bandwidth_array[0] = {1.0 / 1.75, 0., 0.};
    one_over_bandwidth_array[1] = {1.0 / 0.75, 0., 0.};
    one_over_bandwidth_array[2] = {1.0 / 0.25, 0., 0.};
    one_over_bandwidth_array[3] = {1.0 / 1.25, 0., 0.};
    one_over_bandwidth_array[4] = {1.0 / 2.25, 0., 0.};
    one_over_bandwidth_array[5] = {1.0 / 1.25, 0., 0.};
    one_over_bandwidth_array[6] = {1.0 / 0.25, 0., 0.};
    one_over_bandwidth_array[7] = {1.0 / 0.75, 0., 0.};
    one_over_bandwidth_array[8] = {1.0 / 1.75, 0., 0.};
    one_over_bandwidth_array[9] = {1.0 / 2.75, 0., 0.};
    const bool dd = false;
    // 1 bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 3.0;
    const size_t dim = 1;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);

    std::vector<double> bench{0.0159208, 0.0177581, 0.1,      0.157576, 0.15506,
                              0.0164128, 0.01,      0.177581, 0.154304, 0.155386};

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
        10, std::array<double, 3>{1.0 / 4.0, 0., 0.0});

    // lets make the array a little bit more complicated
    const bool dd = false;
    // 1 bin per point
    const size_t n_coarse_bins = 10;
    const double max_window_size = 4.0;
    const size_t dim = 1;
    quick_index qindex(dim, position_array, max_window_size, n_coarse_bins, dd);

    std::vector<double> smooth_result =
        test_kde.reconstruction(data, one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(data, one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(data, log_smooth_result, qindex.domain_decomposed);

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

void test_decomposition(ParallelUnitTest &ut) {
  kde test_kde;

  // test the epan kernel
  double value = test_kde.epan_kernel(0.0);
  if (!rtt_dsxx::soft_equiv(value, 0.75))
    ITFAILS;

  if (rtt_c4::nodes() != 3)
    ITFAILS;

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
        10, std::array<double, 3>{1.0 / 0.1, 0., 0.0});

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, log_smooth_result, qindex.domain_decomposed);

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, log_smooth_result, qindex.domain_decomposed);

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
        10, std::array<double, 3>{1.0 / 4.0, 0., 0.0});

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, log_smooth_result, qindex.domain_decomposed);

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, log_smooth_result, qindex.domain_decomposed);

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
        10, std::array<double, 3>{1.0 / 0.1, 0., 0.0});

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, log_smooth_result, qindex.domain_decomposed);

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, log_smooth_result, qindex.domain_decomposed);

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, log_smooth_result, qindex.domain_decomposed);

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, log_smooth_result, qindex.domain_decomposed);

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, log_smooth_result, qindex.domain_decomposed);

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
        10, std::array<double, 3>{1.0 / 4.0, 0., 0.0});

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, log_smooth_result, qindex.domain_decomposed);

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, log_smooth_result, qindex.domain_decomposed);

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
        10, std::array<double, 3>{1.0, 0., 0.0});
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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);

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
        10, std::array<double, 3>{1.0, 0., 0.0});

    // lets make the array a little bit more complicated
    one_over_bandwidth_array[9] = {1.0 / 0.5, 0., 0.};
    one_over_bandwidth_array[3] = {1.0 / 0.1, 0., 0.};
    one_over_bandwidth_array[4] = {1.0 / 0.5, 0., 0.};
    one_over_bandwidth_array[2] = {1.0 / 2.0, 0., 0.};

    std::vector<double> bench{0.0135142, 0.0160819, 0.0926847, 0.2,      0.1,
                              0.0160819, 0.0378397, 0.160819,  0.122979, 0.2};

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);

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
    one_over_bandwidth_array[9] = {1.0 / 0.5, 1.0 / 4.0, 0.};
    one_over_bandwidth_array[3] = {1.0 / 1.0, 1.0 / 0.1, 0.};
    one_over_bandwidth_array[4] = {1.0 / 0.5, 1.0 / 4.0, 0.};
    one_over_bandwidth_array[2] = {1.0 / 0.1, 1.0 / 4.0, 0.};

    std::vector<double> bench{0.0131256, 0.0158657, 0.1,      0.2,      0.1,
                              0.0158657, 0.0364369, 0.158657, 0.120049, 0.2};

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);

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
        10, std::array<double, 3>{1.0 / 4.0, 0., 0.0});

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

    std::vector<double> smooth_result =
        test_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        test_kde.log_reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    test_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);
    test_kde.apply_conservation(dd_data, log_smooth_result, qindex.domain_decomposed);

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

    std::vector<double> smooth_result =
        refl_kde.reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    std::vector<double> log_smooth_result =
        refl_kde.log_reconstruction(dd_data, dd_one_over_bandwidth_array, qindex);
    // Apply Conservation
    refl_kde.apply_conservation(dd_data, smooth_result, qindex.domain_decomposed);
    refl_kde.apply_conservation(dd_data, log_smooth_result, qindex.domain_decomposed);

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
