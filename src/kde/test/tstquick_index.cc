//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/test/tstquick_index.cc
 * \author Mathew Cleveland
 * \date   Aug. 10th 2021
 * \brief  quick_index testing function
 * \note   Copyright (C) 2021-2021 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "kde/quick_index.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/dbc.hh"
#include <numeric>

using namespace rtt_dsxx;
using namespace rtt_c4;
using namespace rtt_kde;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
//
void test_replication(ParallelUnitTest &ut) {
  {
    std::vector<double> data{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    // in rep mode the max window size does nothing so set it large
    const double max_window_size = 100.0;
    const size_t bins_per_dim = 10UL;
    const bool dd = false;
    const size_t dim = 1;
    quick_index qindex(dim, position_array, max_window_size, bins_per_dim, dd);
    // Check public data
    //------------------------
    if (qindex.domain_decomposed)
      ITFAILS;
    if (qindex.coarse_bin_resolution != bins_per_dim)
      ITFAILS;
    if (!soft_equiv(qindex.max_window_size, max_window_size))
      ITFAILS;
    // Check global bounding box
    if (!soft_equiv(qindex.bounding_box_min[0], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[1], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[2], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[0], 4.5))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[1], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[2], 0.0))
      ITFAILS;
    // Check local coarse_index map
    // build up a global gold to check the map
    std::map<size_t, std::vector<size_t>> gold_map;
    gold_map[0] = {0};
    gold_map[1] = {5};
    gold_map[2] = {1};
    gold_map[3] = {6};
    gold_map[4] = {2};
    gold_map[5] = {7};
    gold_map[6] = {3};
    gold_map[7] = {8};
    gold_map[8] = {4};
    gold_map[9] = {9};
    if (gold_map.size() != qindex.coarse_index_map.size())
      ITFAILS;
    for (auto &map : qindex.coarse_index_map)
      for (size_t i = 0; i < map.second.size(); i++)
        if (gold_map[map.first][i] != map.second[i])
          ITFAILS;

    // Check non-spherical orthogonal distance calculation
    auto distance = qindex.calc_orthogonal_distance({-1, -1, -1}, {1, 1, 1}, 10.0);
    for (auto &val : distance)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(val, 2.0));
  }

  if (ut.numFails == 0) {
    PASSMSG("quick_index checks pass");
  } else {
    FAILMSG("quick_index checks failed");
  }
}

void test_replication_sphere(ParallelUnitTest &ut) {
  {
    std::vector<double> data{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    // in rep mode the max window size does nothing so set it large
    const bool spherical = true;
    const std::array<double, 3> sphere_center{2.0, 0.0, 0.0};
    const double max_window_size = 100.0;
    const size_t bins_per_dim = 10UL;
    const bool dd = false;
    const size_t dim = 2;
    quick_index qindex(dim, position_array, max_window_size, bins_per_dim, dd, spherical,
                       sphere_center);
    // Check public data
    //------------------------
    if (qindex.domain_decomposed)
      ITFAILS;
    if (qindex.coarse_bin_resolution != bins_per_dim)
      ITFAILS;
    if (!soft_equiv(qindex.max_window_size, max_window_size))
      ITFAILS;
    // Check global bounding box
    if (!soft_equiv(qindex.bounding_box_min[0], 0.5))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[1], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[2], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[0], 2.54951, 1e-4))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[1], 5.17604, 1e-4))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[2], 0.0))
      ITFAILS;
    // Check local coarse_index map
    // build up a global gold to check the map
    std::map<size_t, std::vector<size_t>> gold_map;
    gold_map[0] = {2};
    gold_map[23] = {3};
    gold_map[27] = {4};
    gold_map[35] = {8};
    gold_map[39] = {9};
    gold_map[41] = {7};
    gold_map[71] = {6};
    gold_map[85] = {5};
    gold_map[93] = {1};
    gold_map[97] = {0};
    if (gold_map.size() != qindex.coarse_index_map.size())
      ITFAILS;
    for (auto &map : qindex.coarse_index_map)
      for (size_t i = 0; i < map.second.size(); i++)
        if (gold_map[map.first][i] != map.second[i])
          ITFAILS;

    // Check non-spherical orthogonal distance calculation
    auto distance = qindex.calc_orthogonal_distance({-1, 0.5, -1}, {1, 1, 1}, 4.0);
    for (auto &val : distance)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(val, 2.0));
  }

  if (ut.numFails == 0) {
    PASSMSG("quick_index sphere checks pass");
  } else {
    FAILMSG("quick_index sphere checks failed");
  }
}

void test_decomposition(ParallelUnitTest &ut) {
  if (rtt_c4::nodes() != 3)
    ITFAILS;

  int local_size = 3;
  // give the odd size to the final rank to make striding easy
  if (rtt_c4::node() == 2)
    local_size = 4;

  {
    std::vector<double> data{3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    // This cell spatial ordering is difficult for this setup in that every
    // rank requires a sub set of information from every other rank
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::vector<double>> dd_3x_data(3, std::vector<double>(local_size, 0.0));
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_3x_data[0][i] = data[i + rtt_c4::node() * 3];
      dd_3x_data[1][i] = data[i + rtt_c4::node() * 3] + 1;
      dd_3x_data[2][i] = -data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
    }

    // in dd mode the max window size determines the number of ghost cells
    const double max_window_size = 1.0;
    const size_t bins_per_dim = 10UL;
    const bool dd = true;
    const size_t dim = 1;
    quick_index qindex(dim, dd_position_array, max_window_size, bins_per_dim, dd);

    // Check the local state data
    if (!qindex.domain_decomposed)
      ITFAILS;
    if (qindex.coarse_bin_resolution != bins_per_dim)
      ITFAILS;
    if (!soft_equiv(qindex.max_window_size, max_window_size))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[0], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[1], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[2], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[0], 4.5))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[1], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[2], 0.0))
      ITFAILS;
    // Check local coarse_index map
    // local indexing will not match the domain replicated case (different
    // number of points per rank so different local indexing)
    std::map<size_t, std::vector<size_t>> gold_map;
    if (rtt_c4::node() == 0) {
      gold_map[0] = {0}; // 0.0
      gold_map[2] = {1}; // 1.0
      gold_map[4] = {2}; // 2.0
    } else if (rtt_c4::node() == 1) {
      gold_map[6] = {0}; // 3.0
      gold_map[8] = {1}; // 4.0
      gold_map[1] = {2}; // 0.5
    } else {
      gold_map[3] = {0}; // 1.5
      gold_map[5] = {1}; // 2.5
      gold_map[7] = {2}; // 3.5
      gold_map[9] = {3}; // 4.5
    }
    if (gold_map.size() != qindex.coarse_index_map.size())
      ITFAILS;
    for (auto &map : qindex.coarse_index_map)
      for (size_t i = 0; i < map.second.size(); i++)
        if (gold_map[map.first][i] != map.second[i])
          ITFAILS;

    // Check Domain Decomposed Data
    // local bounding box extends beyond local data based on the window size
    if (rtt_c4::node() == 0) {
      if (!soft_equiv(qindex.local_bounding_box_min[0], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[1], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[2], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[0], 2.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[1], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[2], 0.0))
        ITFAILS;
    } else if (rtt_c4::node() == 1) {
      if (!soft_equiv(qindex.local_bounding_box_min[0], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[1], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[2], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[0], 4.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[1], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[2], 0.0))
        ITFAILS;
    } else {
      if (!soft_equiv(qindex.local_bounding_box_min[0], 1.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[1], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[2], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[0], 4.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[1], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[2], 0.0))
        ITFAILS;
    }
    // global bins that span the local domains
    std::vector<size_t> gold_bins;
    if (rtt_c4::node() == 0) {
      gold_bins = {0, 1, 2, 3, 4, 5};
    } else if (rtt_c4::node() == 1) {
      gold_bins = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    } else {
      gold_bins = {2, 3, 4, 5, 6, 7, 8, 9};
    }

    if (gold_bins.size() != qindex.local_bins.size())
      ITFAILS;
    for (size_t i = 0; i < qindex.local_bins.size(); i++)
      if (gold_bins[i] != qindex.local_bins[i])
        ITFAILS;

    // local ghost index map (how to find general location of the ghost data)
    std::map<size_t, std::vector<size_t>> gold_ghost_index_map;
    if (rtt_c4::node() == 0) {
      gold_ghost_index_map[1] = {0}; // 0.5 from rank 1
      gold_ghost_index_map[3] = {1}; // 1.5 from rank 2
      gold_ghost_index_map[5] = {2}; // 2.5 from rank 2
    } else if (rtt_c4::node() == 1) {
      gold_ghost_index_map[0] = {0}; // 0.0 from rank 0
      gold_ghost_index_map[2] = {1}; // 1.0 from rank 0
      gold_ghost_index_map[4] = {2}; // 2.0 from rank 0
      gold_ghost_index_map[3] = {3}; // 1.5 from rank 2
      gold_ghost_index_map[5] = {4}; // 2.5 from rank 2
      gold_ghost_index_map[7] = {5}; // 3.5 from rank 2
      gold_ghost_index_map[9] = {6}; // 4.5 from rank 2
    } else {
      gold_ghost_index_map[2] = {0}; // 1.0 from rank 0
      gold_ghost_index_map[4] = {1}; // 2.0 from rank 0
      gold_ghost_index_map[6] = {2}; // 3.0 from rank 1
      gold_ghost_index_map[8] = {3}; // 4.0 from rank 1
    }
    if (gold_ghost_index_map.size() != qindex.local_ghost_index_map.size())
      ITFAILS;
    for (auto &map : qindex.local_ghost_index_map) {
      if (gold_ghost_index_map[map.first].size() != map.second.size())
        ITFAILS;
      for (size_t i = 0; i < map.second.size(); i++) {
        if (map.second[i] != gold_ghost_index_map[map.first][i])
          ITFAILS;
      }
    }

    // Check the local ghost locations (this tangentially checks the private
    // put_window_map which is used to build this local data).
    std::vector<std::array<double, 3>> gold_ghost_locations;
    if (rtt_c4::node() == 0) {
      gold_ghost_locations = {{0.5, 0.0, 0.0}, {1.5, 0.0, 0.0}, {2.5, 0.0, 0.0}};
    } else if (rtt_c4::node() == 1) {
      gold_ghost_locations = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {1.5, 0.0, 0.0},
                              {2.5, 0.0, 0.0}, {3.5, 0.0, 0.0}, {4.5, 0.0, 0.0}};
    } else {
      gold_ghost_locations = {{1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0}, {4.0, 0.0, 0.0}};
    }
    if (gold_ghost_locations.size() != qindex.local_ghost_locations.size())
      ITFAILS;
    for (size_t i = 0; i < qindex.local_ghost_locations.size(); i++) {
      for (size_t d = 0; d < 3; d++)
        if (!rtt_dsxx::soft_equiv(gold_ghost_locations[i][d], qindex.local_ghost_locations[i][d]))
          ITFAILS;
    }

    // Check collect_ghost_data vector call
    std::vector<double> ghost_data(qindex.local_ghost_buffer_size, 0.0);
    qindex.collect_ghost_data(dd_data, ghost_data);
    std::vector<std::vector<double>> ghost_3x_data(
        3, std::vector<double>(qindex.local_ghost_buffer_size, 0.0));
    qindex.collect_ghost_data(dd_3x_data, ghost_3x_data);

    std::vector<double> gold_ghost_data;
    std::vector<std::vector<double>> gold_3x_ghost_data(3);
    if (rtt_c4::node() == 0) {
      gold_ghost_data = {8.0, 9.0, 10.0};
      gold_3x_ghost_data[0] = {8.0, 9.0, 10.0};
      gold_3x_ghost_data[1] = {9.0, 10.0, 11.0};
      gold_3x_ghost_data[2] = {-8.0, -9.0, -10.0};
    } else if (rtt_c4::node() == 1) {
      gold_ghost_data = {3.0, 4.0, 5.0, 9.0, 10.0, 11.0, 12.0};
      gold_3x_ghost_data[0] = {3.0, 4.0, 5.0, 9.0, 10.0, 11.0, 12.0};
      gold_3x_ghost_data[1] = {4.0, 5.0, 6.0, 10.0, 11.0, 12.0, 13.0};
      gold_3x_ghost_data[2] = {-3.0, -4.0, -5.0, -9.0, -10.0, -11.0, -12.0};
    } else {
      gold_ghost_data = {4.0, 5.0, 6.0, 7.0};
      gold_3x_ghost_data[0] = {4.0, 5.0, 6.0, 7.0};
      gold_3x_ghost_data[1] = {5.0, 6.0, 7.0, 8.0};
      gold_3x_ghost_data[2] = {-4.0, -5.0, -6.0, -7.0};
    }
    for (size_t i = 0; i < ghost_data.size(); i++)
      if (!rtt_dsxx::soft_equiv(ghost_data[i], gold_ghost_data[i]))
        ITFAILS;
    for (size_t i = 0; i < ghost_3x_data[0].size(); i++)
      if (!rtt_dsxx::soft_equiv(ghost_3x_data[0][i], gold_3x_ghost_data[0][i]))
        ITFAILS;
    for (size_t i = 0; i < ghost_3x_data[1].size(); i++)
      if (!rtt_dsxx::soft_equiv(ghost_3x_data[1][i], gold_3x_ghost_data[1][i]))
        ITFAILS;
    for (size_t i = 0; i < ghost_3x_data[2].size(); i++)
      if (!rtt_dsxx::soft_equiv(ghost_3x_data[2][i], gold_3x_ghost_data[2][i]))
        ITFAILS;

    // check max window mapping (more bins then data) functions
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, 0.0, 0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, 0.0, 0.0};
      const std::array<size_t, 3> bin_sizes{5, 0, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "max";
      std::vector<double> window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);

      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 3.0, 0.0, 8.0};
        gold_window_3x_data[0] = {0.0, 0.0, 3.0, 0.0, 8.0};
        gold_window_3x_data[1] = {0.0, 0.0, 4.0, 0.0, 9.0};
        gold_window_3x_data[2] = {0.0, 0.0, -3.0, 0.0, -8.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {10.0, 0.0, 6.0, 0.0, 11.0};
        gold_window_3x_data[0] = {10.0, 0.0, 6.0, 0.0, 11.0};
        gold_window_3x_data[1] = {11.0, 0.0, 7.0, 0.0, 12.0};
        gold_window_3x_data[2] = {-10.0, 0.0, -6.0, 0.0, -11.0};
      } else {
        gold_window_data = {4.0, 0.0, 9.0, 0.0, 5.0};
        gold_window_3x_data[0] = {4.0, 0.0, 9.0, 0.0, 5.0};
        gold_window_3x_data[1] = {5.0, 0.0, 10.0, 0.0, 6.0};
        gold_window_3x_data[2] = {-4.0, 0.0, -9.0, 0.0, -5.0};
      }

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }

    // check min window mapping (more bins then data) functions
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, 0.0, 0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, 0.0, 0.0};
      const std::array<size_t, 3> bin_sizes{5, 0, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "min";
      std::vector<double> window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 3.0, 0.0, 8.0};
        gold_window_3x_data[0] = {0.0, 0.0, 3.0, 0.0, 8.0};
        gold_window_3x_data[1] = {0.0, 0.0, 4.0, 0.0, 9.0};
        gold_window_3x_data[2] = {0.0, 0.0, -3.0, 0.0, -8.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {10.0, 0.0, 6.0, 0.0, 11.0};
        gold_window_3x_data[0] = {10.0, 0.0, 6.0, 0.0, 11.0};
        gold_window_3x_data[1] = {11.0, 0.0, 7.0, 0.0, 12.0};
        gold_window_3x_data[2] = {-10.0, 0.0, -6.0, 0.0, -11.0};
      } else {
        gold_window_data = {4.0, 0.0, 9.0, 0.0, 5.0};
        gold_window_3x_data[0] = {4.0, 0.0, 9.0, 0.0, 5.0};
        gold_window_3x_data[1] = {5.0, 0.0, 10.0, 0.0, 6.0};
        gold_window_3x_data[2] = {-4.0, 0.0, -9.0, 0.0, -5.0};
      }

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }

    // check min_fill window mapping (more bins then data) functions
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, 0.0, 0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, 0.0, 0.0};
      const std::array<size_t, 3> bin_sizes{5, 0, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "min_fill";
      std::vector<double> window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 3.0, 3.0, 8.0};
        gold_window_3x_data[0] = {0.0, 0.0, 3.0, 3.0, 8.0};
        gold_window_3x_data[1] = {0.0, 0.0, 4.0, 4.0, 9.0};
        gold_window_3x_data[2] = {0.0, 0.0, -3.0, -3.0, -8.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {10.0, 10.0, 6.0, 6.0, 11.0};
        gold_window_3x_data[0] = {10.0, 10.0, 6.0, 6.0, 11.0};
        gold_window_3x_data[1] = {11.0, 11.0, 7.0, 7.0, 12.0};
        gold_window_3x_data[2] = {-10.0, -10.0, -6.0, -6.0, -11.0};
      } else {
        gold_window_data = {4.0, 4.0, 9.0, 9.0, 5.0};
        gold_window_3x_data[0] = {4.0, 4.0, 9.0, 9.0, 5.0};
        gold_window_3x_data[1] = {5.0, 5.0, 10.0, 10.0, 6.0};
        gold_window_3x_data[2] = {-4.0, -4.0, -9.0, -9.0, -5.0};
      }

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }

    // check ave window mapping (more bins then data) functions
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, 0.0, 0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, 0.0, 0.0};
      const std::array<size_t, 3> bin_sizes{5, 0, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "ave";
      std::vector<double> window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 3.0, 0.0, 8.0};
        gold_window_3x_data[0] = {0.0, 0.0, 3.0, 0.0, 8.0};
        gold_window_3x_data[1] = {0.0, 0.0, 4.0, 0.0, 9.0};
        gold_window_3x_data[2] = {0.0, 0.0, -3.0, 0.0, -8.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {10.0, 0.0, 6.0, 0.0, 11.0};
        gold_window_3x_data[0] = {10.0, 0.0, 6.0, 0.0, 11.0};
        gold_window_3x_data[1] = {11.0, 0.0, 7.0, 0.0, 12.0};
        gold_window_3x_data[2] = {-10.0, 0.0, -6.0, 0.0, -11.0};
      } else {
        gold_window_data = {4.0, 0.0, 9.0, 0.0, 5.0};
        gold_window_3x_data[0] = {4.0, 0.0, 9.0, 0.0, 5.0};
        gold_window_3x_data[1] = {5.0, 0.0, 10.0, 0.0, 6.0};
        gold_window_3x_data[2] = {-4.0, 0.0, -9.0, 0.0, -5.0};
      }

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }

    // check ave_fill window mapping (more bins then data) functions
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, 0.0, 0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, 0.0, 0.0};
      const std::array<size_t, 3> bin_sizes{5, 0, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "ave_fill";
      std::vector<double> window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 3.0, 3.0, 8.0};
        gold_window_3x_data[0] = {0.0, 0.0, 3.0, 3.0, 8.0};
        gold_window_3x_data[1] = {0.0, 0.0, 4.0, 4.0, 9.0};
        gold_window_3x_data[2] = {0.0, 0.0, -3.0, -3.0, -8.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {10.0, 10.0, 6.0, 6.0, 11.0};
        gold_window_3x_data[0] = {10.0, 10.0, 6.0, 6.0, 11.0};
        gold_window_3x_data[1] = {11.0, 11.0, 7.0, 7.0, 12.0};
        gold_window_3x_data[2] = {-10.0, -10.0, -6.0, -6.0, -11.0};
      } else {
        gold_window_data = {4.0, 4.0, 9.0, 9.0, 5.0};
        gold_window_3x_data[0] = {4.0, 4.0, 9.0, 9.0, 5.0};
        gold_window_3x_data[1] = {5.0, 5.0, 10.0, 10.0, 6.0};
        gold_window_3x_data[2] = {-4.0, -4.0, -9.0, -9.0, -5.0};
      }

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }

    // check normalize window mapping (more bins then data) functions
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, 0.0, 0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, 0.0, 0.0};
      const std::array<size_t, 3> bin_sizes{5, 0, 0};
      const bool normalize = true;
      const bool bias = false;
      const std::string map_type = "ave";
      std::vector<double> window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 3.0 / 11.0, 0.0, 8.0 / 11.0};
        gold_window_3x_data[0] = {0.0, 0.0, 3.0 / 11.0, 0.0, 8.0 / 11.0};
        gold_window_3x_data[1] = {0.0, 0.0, 4.0 / 13.0, 0.0, 9.0 / 13.0};
        gold_window_3x_data[2] = {0.0, 0.0, 3.0 / 11.0, 0.0, 8.0 / 11.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {10.0 / 27.0, 0.0, 6.0 / 27.0, 0.0, 11.0 / 27.0};
        gold_window_3x_data[0] = {10.0 / 27.0, 0.0, 6.0 / 27.0, 0.0, 11.0 / 27.0};
        gold_window_3x_data[1] = {11.0 / 30.0, 0.0, 7.0 / 30.0, 0.0, 12.0 / 30.0};
        gold_window_3x_data[2] = {10.0 / 27.0, 0.0, 6.0 / 27.0, 0.0, 11.0 / 27.0};
      } else {
        gold_window_data = {4.0 / 18.0, 0.0, 9.0 / 18.0, 0.0, 5.0 / 18.0};
        gold_window_3x_data[0] = {4.0 / 18.0, 0.0, 9.0 / 18.0, 0.0, 5.0 / 18.0};
        gold_window_3x_data[1] = {5.0 / 21.0, 0.0, 10.0 / 21.0, 0.0, 6.0 / 21.0};
        gold_window_3x_data[2] = {4.0 / 18.0, 0.0, 9.0 / 18.0, 0.0, 5.0 / 18.0};
      }

      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
    }

    // check ave_fill + normalize window mapping (more bins then data) functions
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, 0.0, 0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, 0.0, 0.0};
      const std::array<size_t, 3> bin_sizes{5, 0, 0};
      const bool normalize = true;
      const bool bias = false;
      const std::string map_type = "ave_fill";
      std::vector<double> window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 3.0 / 14.0, 3.0 / 14.0, 8.0 / 14.0};
        gold_window_3x_data[0] = {0.0, 0.0, 3.0 / 14.0, 3.0 / 14.0, 8.0 / 14.0};
        gold_window_3x_data[1] = {0.0, 0.0, 4.0 / 17.0, 4.0 / 17, 9.0 / 17.0};
        gold_window_3x_data[2] = {0.0, 0.0, 3.0 / 14.0, 3.0 / 14.0, 8.0 / 14.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {10.0 / 43.0, 10.0 / 43.0, 6.0 / 43.0, 6.0 / 43.0, 11.0 / 43.0};
        gold_window_3x_data[0] = {10.0 / 43.0, 10.0 / 43.0, 6.0 / 43.0, 6.0 / 43.0, 11.0 / 43.0};
        gold_window_3x_data[1] = {11.0 / 48.0, 11.0 / 48.0, 7.0 / 48.0, 7.0 / 48.0, 12.0 / 48.0};
        gold_window_3x_data[2] = {10.0 / 43.0, 10.0 / 43.0, 6.0 / 43.0, 6.0 / 43.0, 11.0 / 43.0};
      } else {
        gold_window_data = {4.0 / 31.0, 4.0 / 31.0, 9.0 / 31.0, 9.0 / 31.0, 5.0 / 31.0};
        gold_window_3x_data[0] = {4.0 / 31.0, 4.0 / 31.0, 9.0 / 31.0, 9.0 / 31.0, 5.0 / 31.0};
        gold_window_3x_data[1] = {5.0 / 36.0, 5.0 / 36.0, 10.0 / 36.0, 10.0 / 36.0, 6.0 / 36.0};
        gold_window_3x_data[2] = {4.0 / 31.0, 4.0 / 31.0, 9.0 / 31.0, 9.0 / 31.0, 5.0 / 31.0};
      }

      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
    }

    // check bias window mapping (more bins then data) functions
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, 0.0, 0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, 0.0, 0.0};
      const std::array<size_t, 3> bin_sizes{5, 0, 0};
      const bool normalize = false;
      const bool bias = true;
      const std::string map_type = "ave";
      std::vector<double> window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 3.0, 0.0, 8.0};
        gold_window_3x_data[0] = {0.0, 0.0, 3.0, 0.0, 8.0};
        gold_window_3x_data[1] = {0.0, 0.0, 4.0, 0.0, 9.0};
        gold_window_3x_data[2] = {0.0, 0.0, 8.0 - 3.0, 0.0, 8.0 - 8.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {10.0, 0.0, 6.0, 0.0, 11.0};
        gold_window_3x_data[0] = {10.0, 0.0, 6.0, 0.0, 11.0};
        gold_window_3x_data[1] = {11.0, 0.0, 7.0, 0.0, 12.0};
        gold_window_3x_data[2] = {11.0 - 10.0, 0.0, 11.0 - 6.0, 0.0, 11.0 - 11.0};
      } else {
        gold_window_data = {4.0, 0.0, 9.0, 0.0, 5.0};
        gold_window_3x_data[0] = {4.0, 0.0, 9.0, 0.0, 5.0};
        gold_window_3x_data[1] = {5.0, 0.0, 10.0, 0.0, 6.0};
        gold_window_3x_data[2] = {9.0 - 4.0, 0.0, 9.0 - 9.0, 0.0, 9.0 - 5.0};
      }

      for (size_t v = 0; v < 3; v++) {
        for (size_t i = 0; i < bin_sizes[0]; i++) {
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
        }
      }
      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
    }

    // check bias and normal window mapping (more bins then data) functions
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, 0.0, 0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, 0.0, 0.0};
      const std::array<size_t, 3> bin_sizes{5, 0, 0};
      const bool normalize = true;
      const bool bias = true;
      const std::string map_type = "ave";
      std::vector<double> window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 3.0 / 11.0, 0.0, 8.0 / 11.0};
        gold_window_3x_data[0] = {0.0, 0.0, 3.0 / 11.0, 0.0, 8.0 / 11.0};
        gold_window_3x_data[1] = {0.0, 0.0, 4.0 / 13.0, 0.0, 9.0 / 13.0};
        gold_window_3x_data[2] = {0.0, 0.0, 1.0, 0.0, 0.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {10.0 / 27.0, 0.0, 6.0 / 27.0, 0.0, 11.0 / 27.0};
        gold_window_3x_data[0] = {10.0 / 27.0, 0.0, 6.0 / 27.0, 0.0, 11.0 / 27.0};
        gold_window_3x_data[1] = {11.0 / 30.0, 0.0, 7.0 / 30.0, 0.0, 12.0 / 30.0};
        gold_window_3x_data[2] = {1.0 / 6.0, 0.0, 5.0 / 6.0, 0.0, 0.0};
      } else {
        gold_window_data = {4.0 / 18.0, 0.0, 9.0 / 18.0, 0.0, 5.0 / 18.0};
        gold_window_3x_data[0] = {4.0 / 18.0, 0.0, 9.0 / 18.0, 0.0, 5.0 / 18.0};
        gold_window_3x_data[1] = {5.0 / 21.0, 0.0, 10.0 / 21.0, 0.0, 6.0 / 21.0};
        gold_window_3x_data[2] = {5.0 / 9.0, 0.0, 0.0, 0.0, 4.0 / 9.0};
      }

      for (size_t v = 0; v < 3; v++) {
        for (size_t i = 0; i < bin_sizes[0]; i++) {
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
        }
      }
      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
    }

    // check max mapping (fewer bins then data) functions
    {
      // put two particles in the topmost bin
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, 0.0, 0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, 0.0, 0.0};
      const std::array<size_t, 3> bin_sizes{2, 0, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "max";
      std::vector<double> window_data(2, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(2, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 8.0};
        gold_window_3x_data[0] = {0.0, 8.0};
        gold_window_3x_data[1] = {0.0, 9.0};
        gold_window_3x_data[2] = {0.0, -3.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {10.0, 11.0};
        gold_window_3x_data[0] = {10.0, 11.0};
        gold_window_3x_data[1] = {11.0, 12.0};
        gold_window_3x_data[2] = {-10.0, -6.0};
      } else {
        gold_window_data = {4.0, 9.0};
        gold_window_3x_data[0] = {4.0, 9.0};
        gold_window_3x_data[1] = {5.0, 10.0};
        gold_window_3x_data[2] = {-4.0, -5.0};
      }

      for (size_t v = 0; v < 3; v++) {
        for (size_t i = 0; i < bin_sizes[0]; i++) {
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
        }
      }
      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
    }

    // check min mapping (fewer bins then data) functions
    {
      // put two particles in the topmost bin
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, 0.0, 0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, 0.0, 0.0};
      const std::array<size_t, 3> bin_sizes{2, 0, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "min";
      std::vector<double> window_data(2, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(2, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 3.0};
        gold_window_3x_data[0] = {0.0, 3.0};
        gold_window_3x_data[1] = {0.0, 4.0};
        gold_window_3x_data[2] = {0.0, -8.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {10.0, 6.0};
        gold_window_3x_data[0] = {10.0, 6.0};
        gold_window_3x_data[1] = {11.0, 7.0};
        gold_window_3x_data[2] = {-10.0, -11.0};
      } else {
        gold_window_data = {4.0, 5.0};
        gold_window_3x_data[0] = {4.0, 5.0};
        gold_window_3x_data[1] = {5.0, 6.0};
        gold_window_3x_data[2] = {-4.0, -9.0};
      }

      for (size_t v = 0; v < 3; v++) {
        for (size_t i = 0; i < bin_sizes[0]; i++) {
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
        }
      }
      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
    }

    // check ave mapping (fewer bins then data) functions
    {
      // put two particles in the topmost bin
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, 0.0, 0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, 0.0, 0.0};
      const std::array<size_t, 3> bin_sizes{2, 0, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "ave";
      // use the negative data for the single array operations this time
      std::vector<double> window_data(2, 0.0);
      qindex.map_data_to_grid_window(dd_3x_data[2], ghost_3x_data[2], window_data, min, max,
                                     bin_sizes, map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(2, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);

      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, -5.5};
        gold_window_3x_data[0] = {0.0, 5.5};
        gold_window_3x_data[1] = {0.0, 6.5};
        gold_window_3x_data[2] = {0.0, -5.5};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {-10.0, -8.5};
        gold_window_3x_data[0] = {10.0, 8.5};
        gold_window_3x_data[1] = {11.0, 9.5};
        gold_window_3x_data[2] = {-10.0, -8.5};
      } else {
        gold_window_data = {-4.0, -7.0};
        gold_window_3x_data[0] = {4.0, 7.0};
        gold_window_3x_data[1] = {5.0, 8.0};
        gold_window_3x_data[2] = {-4.0, -7.0};
      }

      for (size_t v = 0; v < 3; v++) {
        for (size_t i = 0; i < bin_sizes[0]; i++) {
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
        }
      }
      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
    }

    // check nearest mapping (fewer bins then data) functions
    // this is the same as average because the simple spacing
    {
      // put two particles in the topmost bin
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, 0.0, 0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, 0.0, 0.0};
      const std::array<size_t, 3> bin_sizes{2, 0, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "nearest";
      // use the negative data for the single array operations this time
      std::vector<double> window_data(2, 0.0);
      qindex.map_data_to_grid_window(dd_3x_data[2], ghost_3x_data[2], window_data, min, max,
                                     bin_sizes, map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(2, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);

      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, -5.5};
        gold_window_3x_data[0] = {0.0, 5.5};
        gold_window_3x_data[1] = {0.0, 6.5};
        gold_window_3x_data[2] = {0.0, -5.5};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {-10.0, -8.5};
        gold_window_3x_data[0] = {10.0, 8.5};
        gold_window_3x_data[1] = {11.0, 9.5};
        gold_window_3x_data[2] = {-10.0, -8.5};
      } else {
        gold_window_data = {-4.0, -7.0};
        gold_window_3x_data[0] = {4.0, 7.0};
        gold_window_3x_data[1] = {5.0, 8.0};
        gold_window_3x_data[2] = {-4.0, -7.0};
      }

      for (size_t v = 0; v < 3; v++) {
        for (size_t i = 0; i < bin_sizes[0]; i++) {
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
        }
      }
      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
    }

    // check nearest mapping (fewer bins then data) functions
    // not the same as average because the window is on the center point
    {
      // put two particles in the topmost bin
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, 0.0, 0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, 0.0, 0.0};
      const std::array<size_t, 3> bin_sizes{1, 0, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "nearest";
      // use the negative data for the single array operations this time
      std::vector<double> window_data(1, 0.0);
      qindex.map_data_to_grid_window(dd_3x_data[2], ghost_3x_data[2], window_data, min, max,
                                     bin_sizes, map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(1, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      if (rtt_c4::node() == 0) {
        gold_window_data = {-3.0};
        gold_window_3x_data[0] = {3.0};
        gold_window_3x_data[1] = {4.0};
        gold_window_3x_data[2] = {-3.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {-6};
        gold_window_3x_data[0] = {6};
        gold_window_3x_data[1] = {7};
        gold_window_3x_data[2] = {-6};
      } else {
        gold_window_data = {-9.0};
        gold_window_3x_data[0] = {9.0};
        gold_window_3x_data[1] = {10.0};
        gold_window_3x_data[2] = {-9.0};
      }

      for (size_t v = 0; v < 3; v++) {
        for (size_t i = 0; i < bin_sizes[0]; i++) {
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
        }
      }
      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
    }

    // check bias and normal window mapping (fewer bins then data) functions
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, 0.0, 0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, 0.0, 0.0};
      const std::array<size_t, 3> bin_sizes{2, 0, 0};
      const bool normalize = true;
      const bool bias = true;
      const std::string map_type = "max";
      std::vector<double> window_data(2, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(2, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);

      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 1.0};
        gold_window_3x_data[0] = {0.0, 1.0};
        gold_window_3x_data[1] = {0.0, 1.0};
        // this one is a bit tricky the bias causes the negative value to go to
        // zero so there is nothing to rescale to one.
        gold_window_3x_data[2] = {0.0, 0.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {10.0 / 21.0, 11.0 / 21.0};
        gold_window_3x_data[0] = {10.0 / 21.0, 11.0 / 21.0};
        gold_window_3x_data[1] = {11.0 / 23.0, 12.0 / 23.0};
        gold_window_3x_data[2] = {0.0, 1.0};
      } else {
        gold_window_data = {4.0 / 13.0, 9.0 / 13.0};
        gold_window_3x_data[0] = {4.0 / 13.0, 9.0 / 13.0};
        gold_window_3x_data[1] = {5.0 / 15.0, 10.0 / 15.0};
        gold_window_3x_data[2] = {1.0, 0.0};
      }

      for (size_t v = 0; v < 3; v++) {
        for (size_t i = 0; i < bin_sizes[0]; i++) {
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
        }
      }
      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
    }
  }

  //2d
  {
    std::vector<double> data{3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    // This cell spatial ordering is difficult for this setup in that every
    // rank requires a sub set of information from every other rank
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::vector<double>> dd_3x_data(3, std::vector<double>(local_size, 0.0));
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_3x_data[0][i] = data[i + rtt_c4::node() * 3];
      dd_3x_data[1][i] = data[i + rtt_c4::node() * 3] + 1;
      dd_3x_data[2][i] = -data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
    }

    // in dd mode the max window size determines the number of ghost cells
    const double max_window_size = 1.0;
    const size_t bins_per_dim = 10UL;
    const bool dd = true;
    const size_t dim = 2;
    quick_index qindex(dim, dd_position_array, max_window_size, bins_per_dim, dd);

    // Check the local state data
    if (!qindex.domain_decomposed)
      ITFAILS;
    if (qindex.coarse_bin_resolution != bins_per_dim)
      ITFAILS;
    if (!soft_equiv(qindex.max_window_size, max_window_size))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[0], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[1], -0.5))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[2], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[0], 4.5))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[1], 0.5))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[2], 0.0))
      ITFAILS;
    // Check local coarse_index map
    // local indexing will not match the domain replicated case (different
    // number of points per rank so different local indexing)
    std::map<size_t, std::vector<size_t>> gold_map;
    if (rtt_c4::node() == 0) {
      gold_map[90] = {0}; // 0.0 0.5
      gold_map[92] = {1}; // 1.0 0.5
      gold_map[94] = {2}; // 2.0 0.5
    } else if (rtt_c4::node() == 1) {
      gold_map[96] = {0}; // 3.0 0.5
      gold_map[98] = {1}; // 4.0 0.5
      gold_map[1] = {2};  // 0.5 -0.5
    } else {
      gold_map[3] = {0}; // 1.5 -0.5
      gold_map[5] = {1}; // 2.5 -0.5
      gold_map[7] = {2}; // 3.5 -0.5
      gold_map[9] = {3}; // 4.5 -0.5
    }
    if (gold_map.size() != qindex.coarse_index_map.size())
      ITFAILS;
    for (auto &map : qindex.coarse_index_map)
      for (size_t i = 0; i < map.second.size(); i++)
        if (gold_map[map.first][i] != map.second[i])
          ITFAILS;

    // Check Domain Decomposed Data
    // local bounding box extends beyond local data based on the window size
    if (rtt_c4::node() == 0) {
      if (!soft_equiv(qindex.local_bounding_box_min[0], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[1], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[2], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[0], 2.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[1], 0.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[2], 0.0))
        ITFAILS;
    } else if (rtt_c4::node() == 1) {
      if (!soft_equiv(qindex.local_bounding_box_min[0], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[1], -0.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[2], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[0], 4.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[1], 0.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[2], 0.0))
        ITFAILS;
    } else {
      if (!soft_equiv(qindex.local_bounding_box_min[0], 1.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[1], -0.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[2], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[0], 4.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[1], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[2], 0.0))
        ITFAILS;
    }
    // global bins that span the local domains
    std::vector<size_t> gold_bins;
    if (rtt_c4::node() == 0) {
      gold_bins = {50, 51, 52, 53, 54, 55, 60, 61, 62, 63, 64, 65, 70, 71, 72,
                   73, 74, 75, 80, 81, 82, 83, 84, 85, 90, 91, 92, 93, 94, 95};
    } else if (rtt_c4::node() == 1) {
      gold_bins = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                   20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                   40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                   60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                   80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99};
    } else {
      gold_bins = {2,  3,  4,  5,  6,  7,  8,  9,  12, 13, 14, 15, 16, 17, 18, 19,
                   22, 23, 24, 25, 26, 27, 28, 29, 32, 33, 34, 35, 36, 37, 38, 39,
                   42, 43, 44, 45, 46, 47, 48, 49, 52, 53, 54, 55, 56, 57, 58, 59};
    }

    if (gold_bins.size() != qindex.local_bins.size())
      ITFAILS;
    for (size_t i = 0; i < qindex.local_bins.size(); i++)
      if (gold_bins[i] != qindex.local_bins[i])
        ITFAILS;

    // local ghost index map (how to find general location of the ghost data)
    std::map<size_t, std::vector<size_t>> gold_ghost_index_map;
    if (rtt_c4::node() == 1) {
      gold_ghost_index_map[90] = {0}; // 0.0, 0.5 from rank 0
      gold_ghost_index_map[92] = {1}; // 1.0, 0.5 from rank 0
      gold_ghost_index_map[94] = {2}; // 2.0, 0.5 from rank 0
      gold_ghost_index_map[3] = {3};  // 1.5, -0.5 from rank 2
      gold_ghost_index_map[5] = {4};  // 2.5, -0.5 from rank 2
      gold_ghost_index_map[7] = {5};  // 3.5, -0.5 from rank 2
      gold_ghost_index_map[9] = {6};  // 4.5, -0.5 from rank 2
    }

    if (gold_ghost_index_map.size() != qindex.local_ghost_index_map.size())
      ITFAILS;
    for (auto &map : qindex.local_ghost_index_map) {
      if (gold_ghost_index_map[map.first].size() != map.second.size())
        ITFAILS;
      for (size_t i = 0; i < map.second.size(); i++) {
        if (map.second[i] != gold_ghost_index_map[map.first][i])
          ITFAILS;
      }
    }

    // Check the local ghost locations (this tangentially checks the private
    // put_window_map which is used to build this local data).
    std::vector<std::array<double, 3>> gold_ghost_locations;
    if (rtt_c4::node() == 1) {
      gold_ghost_locations = {{0.0, 0.5, 0.0},  {1.0, 0.5, 0.0},  {2.0, 0.5, 0.0}, {1.5, -0.5, 0.0},
                              {2.5, -0.5, 0.0}, {3.5, -0.5, 0.0}, {4.5, -0.5, 0.0}};
    }
    if (gold_ghost_locations.size() != qindex.local_ghost_locations.size())
      ITFAILS;
    for (size_t i = 0; i < qindex.local_ghost_locations.size(); i++) {
      for (size_t d = 0; d < 3; d++)
        if (!rtt_dsxx::soft_equiv(gold_ghost_locations[i][d], qindex.local_ghost_locations[i][d]))
          ITFAILS;
    }

    // Check collect_ghost_data vector call
    std::vector<double> ghost_data(qindex.local_ghost_buffer_size, 0.0);
    qindex.collect_ghost_data(dd_data, ghost_data);
    std::vector<std::vector<double>> ghost_3x_data(
        3, std::vector<double>(qindex.local_ghost_buffer_size, 0.0));
    qindex.collect_ghost_data(dd_3x_data, ghost_3x_data);

    std::vector<double> gold_ghost_data;
    std::vector<std::vector<double>> gold_3x_ghost_data(3);
    if (rtt_c4::node() == 1) {
      gold_ghost_data = {3.0, 4.0, 5.0, 9.0, 10.0, 11.0, 12.0};
      gold_3x_ghost_data[0] = {3.0, 4.0, 5.0, 9.0, 10.0, 11.0, 12.0};
      gold_3x_ghost_data[1] = {4.0, 5.0, 6.0, 10.0, 11.0, 12.0, 13.0};
      gold_3x_ghost_data[2] = {-3.0, -4.0, -5.0, -9.0, -10.0, -11.0, -12.0};
    }

    for (size_t i = 0; i < ghost_data.size(); i++)
      if (!rtt_dsxx::soft_equiv(ghost_data[i], gold_ghost_data[i]))
        ITFAILS;
    for (size_t i = 0; i < ghost_3x_data[0].size(); i++)
      if (!rtt_dsxx::soft_equiv(ghost_3x_data[0][i], gold_3x_ghost_data[0][i]))
        ITFAILS;
    for (size_t i = 0; i < ghost_3x_data[1].size(); i++)
      if (!rtt_dsxx::soft_equiv(ghost_3x_data[1][i], gold_3x_ghost_data[1][i]))
        ITFAILS;
    for (size_t i = 0; i < ghost_3x_data[2].size(); i++)
      if (!rtt_dsxx::soft_equiv(ghost_3x_data[2][i], gold_3x_ghost_data[2][i]))
        ITFAILS;

    // check max window mapping (more bins then data) functions
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, dd_position_array[0][1] - 0.5,
                                      0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, dd_position_array[0][1] + 0.5,
                                      0.0};
      const std::array<size_t, 3> bin_sizes{5, 1, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "max";
      std::vector<double> window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      // different result then 1D because the 1.0 y offset of the data
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 3.0, 0.0, 0.0};
        gold_window_3x_data[0] = {0.0, 0.0, 3.0, 0.0, 0.0};
        gold_window_3x_data[1] = {0.0, 0.0, 4.0, 0.0, 0.0};
        gold_window_3x_data[2] = {0.0, 0.0, -3.0, 0.0, 0.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {0.0, 0.0, 6.0, 0.0, 0.0};
        gold_window_3x_data[0] = {0.0, 0.0, 6.0, 0.0, 0.0};
        gold_window_3x_data[1] = {0.0, 0.0, 7.0, 0.0, 0.0};
        gold_window_3x_data[2] = {0.0, 0.0, -6.0, 0.0, 0.0};
      } else {
        gold_window_data = {0.0, 0.0, 9.0, 0.0, 0.0};
        gold_window_3x_data[0] = {0.0, 0.0, 9.0, 0.0, 0.0};
        gold_window_3x_data[1] = {0.0, 0.0, 10.0, 0.0, 0.0};
        gold_window_3x_data[2] = {0.0, 0.0, -9.0, 0.0, 0.0};
      }

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }

    // check max_fill window mapping (more bins then data) functions
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{dd_position_array[0][0] - 0.5, dd_position_array[0][1] - 0.5,
                                      0.0};
      const std::array<double, 3> max{dd_position_array[0][0] + 0.5, dd_position_array[0][1] + 0.5,
                                      0.0};
      const std::array<size_t, 3> bin_sizes{5, 1, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "max_fill";
      std::vector<double> window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      // different result then 1D because the 1.0 y offset of the data
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 3.0, 3.0, 3.0};
        gold_window_3x_data[0] = {0.0, 0.0, 3.0, 3.0, 3.0};
        gold_window_3x_data[1] = {0.0, 0.0, 4.0, 4.0, 4.0};
        gold_window_3x_data[2] = {0.0, 0.0, -3.0, -3.0, -3.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {0.0, 0.0, 6.0, 6.0, 6.0};
        gold_window_3x_data[0] = {0.0, 0.0, 6.0, 6.0, 6.0};
        gold_window_3x_data[1] = {0.0, 0.0, 7.0, 7.0, 7.0};
        gold_window_3x_data[2] = {0.0, 0.0, -6.0, -6.0, -6.0};
      } else {
        gold_window_data = {0.0, 0.0, 9.0, 9.0, 9.0};
        gold_window_3x_data[0] = {0.0, 0.0, 9.0, 9.0, 9.0};
        gold_window_3x_data[1] = {0.0, 0.0, 10.0, 10.0, 10.0};
        gold_window_3x_data[2] = {0.0, 0.0, -9.0, -9.0, -9.0};
      }

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }

    // check y dim mapping only on rank 1 because the ghost data extends far enough to touch all domain space
    if (rtt_c4::node() == 1) {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{0.5, -0.5, 0.0};
      const std::array<double, 3> max{1.5, 0.5, 0.0};
      const std::array<size_t, 3> bin_sizes{1, 2, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "min";
      std::vector<double> window_data(2, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(2, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      gold_window_data = {8.0, 4.0};
      gold_window_3x_data[0] = {8.0, 4.0};
      gold_window_3x_data[1] = {9.0, 5.0};
      gold_window_3x_data[2] = {-9.0, -4.0};

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }

    // check y dim mapping only on rank 1 because the ghost data extends far enough to touch all domain space
    if (rtt_c4::node() == 1) {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{0.5, -0.5, 0.0};
      const std::array<double, 3> max{1.5, 0.5, 0.0};
      const std::array<size_t, 3> bin_sizes{1, 2, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "ave";
      std::vector<double> window_data(2, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(2, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      gold_window_data = {8.5, 4.0};
      gold_window_3x_data[0] = {8.5, 4.0};
      gold_window_3x_data[1] = {9.5, 5.0};
      gold_window_3x_data[2] = {-8.5, -4.0};

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }

    // check y dim mapping only on rank 1 because the ghost data extends far enough to touch all domain space
    if (rtt_c4::node() == 1) {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{0.45, -0.5, 0.0};
      const std::array<double, 3> max{0.55, 0.5, 0.0};
      const std::array<size_t, 3> bin_sizes{1, 2, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "ave";
      std::vector<double> window_data(2, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(2, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);

      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      gold_window_data = {8.0, 0.0};
      gold_window_3x_data[0] = {8.0, 0.0};
      gold_window_3x_data[1] = {9.0, 0.0};
      gold_window_3x_data[2] = {-8.0, 0.0};

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }

    // check y dim mapping only on rank 1 because the ghost data extends far enough to touch all domain space
    if (rtt_c4::node() == 1) {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> min{0.5, -0.5, 0.0};
      const std::array<double, 3> max{1.5, 0.5, 0.0};
      const std::array<size_t, 3> bin_sizes{1, 2, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "max";
      std::vector<double> window_data(2, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, window_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<std::vector<double>> window_3x_data(3, std::vector<double>(2, 0.0));
      qindex.map_data_to_grid_window(dd_3x_data, ghost_3x_data, window_3x_data, min, max, bin_sizes,
                                     map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      gold_window_data = {9.0, 4.0};
      gold_window_3x_data[0] = {9.0, 4.0};
      gold_window_3x_data[1] = {10.0, 5.0};
      gold_window_3x_data[2] = {-8.0, -4.0};

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }
  }

  if (ut.numFails == 0) {
    PASSMSG("quick_index DD checks pass");
  } else {
    FAILMSG("quick_index DD checks failed");
  }
}

void test_decomposition_sphere(ParallelUnitTest &ut) {
  if (rtt_c4::nodes() != 3)
    ITFAILS;

  //2d sphere
  {

    const size_t local_size = 4;
    const std::array<double, 3> sphere_center{0.0, 0.0, 0.0};
    const std::array<double, 2> radial_edges{0.5, 1.0};
    const std::array<double, 6> cosine_edges{-.99, 0, .99, -.99, 0, .99};
    const size_t data_size = radial_edges.size() * cosine_edges.size();
    std::vector<std::array<double, 3>> position_array(data_size,
                                                      std::array<double, 3>{0.0, 0.0, 0.0});

    std::vector<double> shell_data(data_size, 0.0);
    std::vector<double> spoke_data(data_size, 0.0);
    size_t point_i = 0;
    size_t mui = 0;
    for (auto &mu : cosine_edges) {
      size_t ri = 0;
      for (auto &r : radial_edges) {
        double sign = mui < 3 ? 1.0 : -1.0;
        spoke_data[point_i] = static_cast<double>(mui) + 1.0;
        shell_data[point_i] = static_cast<double>(ri) + 1.0;
        double rel_y = r * mu;
        position_array[point_i][0] =
            rtt_dsxx::soft_equiv(r * r, rel_y * rel_y, 1e-6) ? 0.0 : sqrt(r * r - rel_y * rel_y);
        position_array[point_i][0] *= sign;
        position_array[point_i][1] = sphere_center[1] + rel_y;
        point_i++;
        ri++;
      }
      mui++;
    }

    std::vector<double> dd_data(local_size, 1.0);
    std::vector<std::vector<double>> dd_3x_data(3, std::vector<double>(local_size, 1.0));
    std::vector<double> dd_spoke_data(local_size);
    std::vector<double> dd_shell_data(local_size);
    std::vector<std::array<double, 3>> dd_position_array(local_size);
    for (size_t i = 0; i < local_size; i++) {
      dd_spoke_data[i] = spoke_data[i + rtt_c4::node() * local_size];
      dd_shell_data[i] = shell_data[i + rtt_c4::node() * local_size];
      dd_position_array[i] = position_array[i + rtt_c4::node() * local_size];
    }

    // in dd mode the max window size does nothing so set it large
    const bool spherical = true;
    const double max_window_size = 1.0;
    const size_t bins_per_dim = 10UL;
    const bool dd = true;
    const size_t dim = 2;
    quick_index qindex(dim, dd_position_array, max_window_size, bins_per_dim, dd, spherical,
                       sphere_center);
    // Check the local state data
    if (!qindex.domain_decomposed)
      ITFAILS;
    if (qindex.coarse_bin_resolution != bins_per_dim)
      ITFAILS;
    if (!soft_equiv(qindex.max_window_size, max_window_size))
      ITFAILS;

    // Check global bounding box
    if (!soft_equiv(qindex.bounding_box_min[0], 0.5))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[1], 0.141539, 1e-4))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[2], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[0], 1.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[1], 6.14165, 1e-4))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[2], 0.0))
      ITFAILS;
    // Check local coarse_index map
    // local indexing will not match the domain replicated case (different
    // number of points per rank so different local indexing)
    std::map<size_t, std::vector<size_t>> gold_map;
    if (rtt_c4::node() == 0) {
      gold_map[20] = {2};
      gold_map[29] = {3};
      gold_map[40] = {0};
      gold_map[49] = {1};
    } else if (rtt_c4::node() == 1) {
      gold_map[0] = {0};
      gold_map[9] = {1};
      gold_map[50] = {2};
      gold_map[59] = {3};
    } else {
      gold_map[70] = {0};
      gold_map[79] = {1};
      gold_map[90] = {2};
      gold_map[99] = {3};
    }
    if (gold_map.size() != qindex.coarse_index_map.size())
      ITFAILS;

    for (auto &map : qindex.coarse_index_map)
      for (size_t i = 0; i < map.second.size(); i++)
        if (gold_map[map.first][i] != map.second[i])
          ITFAILS;

    // Check Domain Decomposed Data
    // local bounding box extends beyond local data based on the window size
    if (rtt_c4::node() == 0) {
      if (!soft_equiv(qindex.local_bounding_box_min[0], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[1], 1.23746, 1e-4))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[2], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[0], 1.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[1], 3.33339, 1e-4))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[2], 0.0))
        ITFAILS;
    } else if (rtt_c4::node() == 1) {
      if (!soft_equiv(qindex.local_bounding_box_min[0], 0.0))
        ITFAILS;
      // overlaps the theta=0=pi*2 boundary
      if (!soft_equiv(qindex.local_bounding_box_min[1], -0.191794, 1e-4))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[2], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[0], 1.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[1], 3.61647, 1e-4))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[2], 0.0))
        ITFAILS;
    } else {
      if (!soft_equiv(qindex.local_bounding_box_min[0], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[1], 4.37906, 1e-4))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[2], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[0], 1.5))
        ITFAILS;
      // overlaps the theta=0=pi*2 boundary
      if (!soft_equiv(qindex.local_bounding_box_max[1], 6.47498, 1e-4))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[2], 0.0))
        ITFAILS;
    }
    // global bins that span the local domains
    std::vector<size_t> gold_bins;
    if (rtt_c4::node() == 0) {
      gold_bins = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                   27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
                   44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59};
    } else if (rtt_c4::node() == 1) {
      gold_bins = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17,
                   18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
                   36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,
                   54, 55, 56, 57, 58, 59, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99};
    } else {
      gold_bins = {70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
                   90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 0,  1,  2,  3,  4,  5,  6,  7,  8,  9};
    }

    if (gold_bins.size() != qindex.local_bins.size())
      ITFAILS;
    for (size_t i = 0; i < qindex.local_bins.size(); i++)
      if (gold_bins[i] != qindex.local_bins[i])
        ITFAILS;

    // local ghost index map (how to find general location of the ghost data)
    std::map<size_t, std::vector<size_t>> gold_ghost_index_map;
    if (rtt_c4::node() == 0) {
      gold_ghost_index_map[50] = {0};
      gold_ghost_index_map[59] = {1};
    } else if (rtt_c4::node() == 1) {
      gold_ghost_index_map[20] = {0};
      gold_ghost_index_map[29] = {1};
      gold_ghost_index_map[40] = {2};
      gold_ghost_index_map[49] = {3};
      gold_ghost_index_map[90] = {4};
      gold_ghost_index_map[99] = {5};
    } else {
      gold_ghost_index_map[0] = {0};
      gold_ghost_index_map[9] = {1};
    }

    if (gold_ghost_index_map.size() != qindex.local_ghost_index_map.size())
      ITFAILS;
    for (auto &map : qindex.local_ghost_index_map) {
      if (gold_ghost_index_map[map.first].size() != map.second.size())
        ITFAILS;
      for (size_t i = 0; i < map.second.size(); i++) {
        if (map.second[i] != gold_ghost_index_map[map.first][i])
          ITFAILS;
      }
    }

    // Check the local ghost locations (this tangentially checks the private
    // put_window_map which is used to build this local data).
    std::vector<std::array<double, 3>> gold_ghost_locations;
    if (rtt_c4::node() == 0) {
      gold_ghost_locations = {{0.5, 3.28313, 0.0}, {1, 3.28313, 0.0}};
    } else if (rtt_c4::node() == 1) {
      gold_ghost_locations = {{0.5, 1.5708, 0.0}, {1, 1.5708, 0.0},    {0.5, 3.00005, 0.0},
                              {1, 3.00005, 0.0},  {0.5, 6.14165, 0.0}, {1, 6.14165, 0.0}};
    } else {
      gold_ghost_locations = {{0.5, 0.141539, 0.0}, {1, 0.141539, 0.0}};
    }

    if (gold_ghost_locations.size() != qindex.local_ghost_locations.size())
      ITFAILS;

    for (size_t i = 0; i < qindex.local_ghost_locations.size(); i++) {
      for (size_t d = 0; d < 2; d++)
        if (!rtt_dsxx::soft_equiv(gold_ghost_locations[i][d], qindex.local_ghost_locations[i][d],
                                  1e-4))
          ITFAILS;
    }
  }

  //2d half sphere
  {

    const size_t local_size = 24;
    const std::array<double, 3> sphere_center{0.0, -1.0, 0.0};
    const std::array<double, 8> radial_edges{0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0};
    const std::array<double, 9> cosine_edges{-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0};
    const size_t data_size = radial_edges.size() * cosine_edges.size();
    std::vector<std::array<double, 3>> position_array(data_size,
                                                      std::array<double, 3>{0.0, 0.0, 0.0});

    std::vector<double> shell_data(data_size, 0.0);
    std::vector<double> spoke_data(data_size, 0.0);
    size_t point_i = 0;
    size_t ri = 0;
    for (auto &r : radial_edges) {
      size_t mui = 0;
      for (auto &mu : cosine_edges) {
        spoke_data[point_i] = static_cast<double>(mui) + 1.0;
        shell_data[point_i] = static_cast<double>(ri) + 1.0;
        double rel_y = r * mu;
        position_array[point_i][0] =
            rtt_dsxx::soft_equiv(r * r, rel_y * rel_y, 1e-6) ? 0.0 : sqrt(r * r - rel_y * rel_y);
        position_array[point_i][1] = sphere_center[1] + rel_y;
        point_i++;
        mui++;
      }
      ri++;
    }

    std::vector<double> dd_data(local_size, 1.0);
    std::vector<std::vector<double>> dd_2x_data(2, std::vector<double>(local_size, 1.0));
    std::vector<std::array<double, 3>> dd_position_array(local_size);
    for (size_t i = 0; i < local_size; i++) {
      dd_data[i] = shell_data[i + rtt_c4::node() * local_size];
      dd_2x_data[0][i] = shell_data[i + rtt_c4::node() * local_size];
      dd_2x_data[1][i] = spoke_data[i + rtt_c4::node() * local_size];
      dd_position_array[i] = position_array[i + rtt_c4::node() * local_size];
    }

    // in dd mode the max window size does nothing so set it large
    const bool spherical = true;
    const double max_window_size = 1.0;
    const size_t bins_per_dim = 10UL;
    const bool dd = true;
    const size_t dim = 2;
    quick_index qindex(dim, dd_position_array, max_window_size, bins_per_dim, dd, spherical,
                       sphere_center);
    // Check the local state data
    if (!qindex.domain_decomposed)
      ITFAILS;
    if (qindex.coarse_bin_resolution != bins_per_dim)
      ITFAILS;
    if (!soft_equiv(qindex.max_window_size, max_window_size))
      ITFAILS;

    // Check global bounding box
    if (!soft_equiv(qindex.bounding_box_min[0], 0.025))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[1], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[2], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[0], 1.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[1], rtt_units::PI, 1e-6))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[2], 0.0))
      ITFAILS;
    // Check local coarse_index map
    // local indexing will not match the domain replicated case (different
    // number of points per rank so different local indexing)
    std::map<size_t, std::vector<size_t>> gold_map;
    if (rtt_c4::node() == 0) {
      gold_map[0] = {8, 17};
      gold_map[20] = {7, 16};
      gold_map[30] = {6, 15};
      gold_map[40] = {5, 14, 23};
      gold_map[50] = {3, 4, 12, 13, 21, 22};
      gold_map[60] = {2, 11, 20};
      gold_map[70] = {1, 10, 19};
      gold_map[90] = {0, 9, 18};
    } else if (rtt_c4::node() == 1) {
      gold_map[0] = {2, 11};
      gold_map[2] = {20};
      gold_map[20] = {1, 10};
      gold_map[22] = {19};
      gold_map[30] = {0, 9};
      gold_map[32] = {18};
      gold_map[40] = {8};
      gold_map[42] = {17};
      gold_map[50] = {6, 7};
      gold_map[52] = {15, 16};
      gold_map[60] = {5};
      gold_map[62] = {14};
      gold_map[64] = {23};
      gold_map[70] = {4};
      gold_map[72] = {13};
      gold_map[74] = {22};
      gold_map[90] = {3};
      gold_map[92] = {12};
      gold_map[94] = {21};
    } else {
      gold_map[4] = {5};
      gold_map[7] = {14};
      gold_map[9] = {23};
      gold_map[24] = {4};
      gold_map[27] = {13};
      gold_map[29] = {22};
      gold_map[34] = {3};
      gold_map[37] = {12};
      gold_map[39] = {21};
      gold_map[44] = {2};
      gold_map[47] = {11};
      gold_map[49] = {20};
      gold_map[54] = {0, 1};
      gold_map[57] = {9, 10};
      gold_map[59] = {18, 19};
      gold_map[67] = {8};
      gold_map[69] = {17};
      gold_map[77] = {7};
      gold_map[79] = {16};
      gold_map[97] = {6};
      gold_map[99] = {15};
    }
    if (gold_map.size() != qindex.coarse_index_map.size())
      ITFAILS;
    for (auto &map : qindex.coarse_index_map)
      for (size_t i = 0; i < map.second.size(); i++)
        if (gold_map[map.first][i] != map.second[i])
          ITFAILS;

    // Check Domain Decomposed Data
    // local bounding box extends beyond local data based on the window size
    if (rtt_c4::node() == 0) {
      if (!soft_equiv(qindex.local_bounding_box_min[0], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[1], -0.5 / 0.575))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[2], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[0], 0.575))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[1], rtt_units::PI + 0.5 / 0.575))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[2], 0.0))
        ITFAILS;
    } else if (rtt_c4::node() == 1) {
      if (!soft_equiv(qindex.local_bounding_box_min[0], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[1], -0.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[2], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[0], 1.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[1], rtt_units::PI + 0.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[2], 0.0))
        ITFAILS;
    } else {
      if (!soft_equiv(qindex.local_bounding_box_min[0], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[1], -0.5 / 1.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[2], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[0], 1.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[1], rtt_units::PI + 0.5 / 1.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[2], 0.0))
        ITFAILS;
    }
    // global bins that span the local domains
    std::vector<size_t> gold_bins;
    if (rtt_c4::node() == 0) {
      gold_bins = {0,  1,  2,  3,  4,  5,  10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, 25, 30, 31,
                   32, 33, 34, 35, 40, 41, 42, 43, 44, 45, 50, 51, 52, 53, 54, 55, 60, 61, 62, 63,
                   64, 65, 70, 71, 72, 73, 74, 75, 80, 81, 82, 83, 84, 85, 90, 91, 92, 93, 94, 95};
    } else if (rtt_c4::node() == 1) {
      gold_bins = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                   20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                   40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                   60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                   80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99};
    } else {
      gold_bins = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                   20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                   40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                   60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                   80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99};
    }

    if (gold_bins.size() != qindex.local_bins.size())
      ITFAILS;
    for (size_t i = 0; i < qindex.local_bins.size(); i++)
      if (gold_bins[i] != qindex.local_bins[i])
        ITFAILS;

    // local ghost index map (how to find general location of the ghost data)
    std::map<size_t, std::vector<size_t>> gold_ghost_index_map;
    if (rtt_c4::node() == 0) {
      gold_ghost_index_map[0] = {0, 1};
      gold_ghost_index_map[2] = {2};
      gold_ghost_index_map[4] = {24};
      gold_ghost_index_map[20] = {3, 4};
      gold_ghost_index_map[22] = {5};
      gold_ghost_index_map[24] = {25};
      gold_ghost_index_map[30] = {6, 7};
      gold_ghost_index_map[32] = {8};
      gold_ghost_index_map[34] = {26};
      gold_ghost_index_map[40] = {9};
      gold_ghost_index_map[42] = {10};
      gold_ghost_index_map[44] = {27};
      gold_ghost_index_map[50] = {11, 12};
      gold_ghost_index_map[52] = {13, 14};
      gold_ghost_index_map[54] = {28, 29};
      gold_ghost_index_map[60] = {15};
      gold_ghost_index_map[62] = {16};
      gold_ghost_index_map[64] = {17};
      gold_ghost_index_map[70] = {18};
      gold_ghost_index_map[72] = {19};
      gold_ghost_index_map[74] = {20};
      gold_ghost_index_map[90] = {21};
      gold_ghost_index_map[92] = {22};
      gold_ghost_index_map[94] = {23};
    } else if (rtt_c4::node() == 1) {
      gold_ghost_index_map[0] = {0, 1};
      gold_ghost_index_map[4] = {24};
      gold_ghost_index_map[7] = {25};
      gold_ghost_index_map[9] = {26};
      gold_ghost_index_map[20] = {2, 3};
      gold_ghost_index_map[24] = {27};
      gold_ghost_index_map[27] = {28};
      gold_ghost_index_map[29] = {29};
      gold_ghost_index_map[30] = {4, 5};
      gold_ghost_index_map[34] = {30};
      gold_ghost_index_map[37] = {31};
      gold_ghost_index_map[39] = {32};
      gold_ghost_index_map[40] = {6, 7, 8};
      gold_ghost_index_map[44] = {33};
      gold_ghost_index_map[47] = {34};
      gold_ghost_index_map[49] = {35};
      gold_ghost_index_map[50] = {9, 10, 11, 12, 13, 14};
      gold_ghost_index_map[54] = {36, 37};
      gold_ghost_index_map[57] = {38, 39};
      gold_ghost_index_map[59] = {40, 41};
      gold_ghost_index_map[60] = {15, 16, 17};
      gold_ghost_index_map[67] = {42};
      gold_ghost_index_map[69] = {43};
      gold_ghost_index_map[70] = {18, 19, 20};
      gold_ghost_index_map[77] = {44};
      gold_ghost_index_map[79] = {45};
      gold_ghost_index_map[90] = {21, 22, 23};
      gold_ghost_index_map[97] = {46};
      gold_ghost_index_map[99] = {47};
    } else {
      gold_ghost_index_map[0] = {0, 1, 24, 25};
      gold_ghost_index_map[2] = {26};
      gold_ghost_index_map[20] = {2, 3, 27, 28};
      gold_ghost_index_map[22] = {29};
      gold_ghost_index_map[30] = {4, 5, 30, 31};
      gold_ghost_index_map[32] = {32};
      gold_ghost_index_map[40] = {6, 7, 8, 33};
      gold_ghost_index_map[42] = {34};
      gold_ghost_index_map[50] = {9, 10, 11, 12, 13, 14, 35, 36};
      gold_ghost_index_map[52] = {37, 38};
      gold_ghost_index_map[60] = {15, 16, 17, 39};
      gold_ghost_index_map[62] = {40};
      gold_ghost_index_map[64] = {41};
      gold_ghost_index_map[70] = {18, 19, 20, 42};
      gold_ghost_index_map[72] = {43};
      gold_ghost_index_map[74] = {44};
      gold_ghost_index_map[90] = {21, 22, 23, 45};
      gold_ghost_index_map[92] = {46};
      gold_ghost_index_map[94] = {47};
    }

    if (gold_ghost_index_map.size() != qindex.local_ghost_index_map.size())
      ITFAILS;
    for (auto &map : qindex.local_ghost_index_map) {
      if (gold_ghost_index_map[map.first].size() != map.second.size())
        ITFAILS;
      for (size_t i = 0; i < map.second.size(); i++) {
        if (map.second[i] != gold_ghost_index_map[map.first][i])
          ITFAILS;
      }
    }

    // Check the local ghost locations (this tangentially checks the private
    // put_window_map which is used to build this local data).
    std::vector<std::array<double, 3>> gold_ghost_locations;
    if (rtt_c4::node() == 0) {
      gold_ghost_locations = {
          {0.075, 0, 0.0},      {0.1, 0, 0.0},         {0.25, 0, 0.0},       {0.075, 0.722734, 0.0},
          {0.1, 0.722734, 0.0}, {0.25, 0.722734, 0.0}, {0.075, 1.0472, 0.0}, {0.1, 1.0472, 0.0},
          {0.25, 1.0472, 0.0},  {0.1, 1.31812, 0.0},   {0.25, 1.31812, 0.0}, {0.1, 1.82348, 0.0},
          {0.1, 1.5708, 0.0},   {0.25, 1.82348, 0.0},  {0.25, 1.5708, 0.0},  {0.1, 2.0944, 0.0},
          {0.25, 2.0944, 0.0},  {0.5, 2.0944, 0.0},    {0.1, 2.41886, 0.0},  {0.25, 2.41886, 0.0},
          {0.5, 2.41886, 0.0},  {0.1, 3.14159, 0.0},   {0.25, 3.14159, 0.0}, {0.5, 3.14159, 0.0},
          {0.5, 0, 0.0},        {0.5, 0.722734, 0.0},  {0.5, 1.0472, 0.0},   {0.5, 1.31812, 0.0},
          {0.5, 1.82348, 0.0},  {0.5, 1.5708, 0.0}};
    } else if (rtt_c4::node() == 1) {
      gold_ghost_locations = {{0.025, 0, 0.0},       {0.05, 0, 0.0},        {0.025, 0.722734, 0.0},
                              {0.05, 0.722734, 0.0}, {0.025, 1.0472, 0.0},  {0.05, 1.0472, 0.0},
                              {0.025, 1.31812, 0.0}, {0.05, 1.31812, 0.0},  {0.075, 1.31812, 0.0},
                              {0.025, 1.82348, 0.0}, {0.025, 1.5708, 0.0},  {0.05, 1.82348, 0.0},
                              {0.05, 1.5708, 0.0},   {0.075, 1.82348, 0.0}, {0.075, 1.5708, 0.0},
                              {0.025, 2.0944, 0.0},  {0.05, 2.0944, 0.0},   {0.075, 2.0944, 0.0},
                              {0.025, 2.41886, 0.0}, {0.05, 2.41886, 0.0},  {0.075, 2.41886, 0.0},
                              {0.025, 3.14159, 0.0}, {0.05, 3.14159, 0.0},  {0.075, 3.14159, 0.0},
                              {0.5, 0, 0.0},         {0.75, 0, 0.0},        {1, 0, 0.0},
                              {0.5, 0.722734, 0.0},  {0.75, 0.722734, 0.0}, {1, 0.722734, 0.0},
                              {0.5, 1.0472, 0.0},    {0.75, 1.0472, 0.0},   {1, 1.0472, 0.0},
                              {0.5, 1.31812, 0.0},   {0.75, 1.31812, 0.0},  {1, 1.31812, 0.0},
                              {0.5, 1.82348, 0.0},   {0.5, 1.5708, 0.0},    {0.75, 1.82348, 0.0},
                              {0.75, 1.5708, 0.0},   {1, 1.82348, 0.0},     {1, 1.5708, 0.0},
                              {0.75, 2.0944, 0.0},   {1, 2.0944, 0.0},      {0.75, 2.41886, 0.0},
                              {1, 2.41886, 0.0},     {0.75, 3.14159, 0.0},  {1, 3.14159, 0.0}};
    } else {
      gold_ghost_locations = {{0.025, 0, 0.0},        {0.05, 0, 0.0},        {0.025, 0.722734, 0.0},
                              {0.05, 0.722734, 0.0},  {0.025, 1.0472, 0.0},  {0.05, 1.0472, 0.0},
                              {0.025, 1.31812, 0.0},  {0.05, 1.31812, 0.0},  {0.075, 1.31812, 0.0},
                              {0.025, 1.82348, 0.0},  {0.025, 1.5708, 0.0},  {0.05, 1.82348, 0.0},
                              {0.05, 1.5708, 0.0},    {0.075, 1.82348, 0.0}, {0.075, 1.5708, 0.0},
                              {0.025, 2.0944, 0.0},   {0.05, 2.0944, 0.0},   {0.075, 2.0944, 0.0},
                              {0.025, 2.41886, 0.0},  {0.05, 2.41886, 0.0},  {0.075, 2.41886, 0.0},
                              {0.025, 3.14159, 0.0},  {0.05, 3.14159, 0.0},  {0.075, 3.14159, 0.0},
                              {0.075, 0, 0.0},        {0.1, 0, 0.0},         {0.25, 0, 0.0},
                              {0.075, 0.722734, 0.0}, {0.1, 0.722734, 0.0},  {0.25, 0.722734, 0.0},
                              {0.075, 1.0472, 0.0},   {0.1, 1.0472, 0.0},    {0.25, 1.0472, 0.0},
                              {0.1, 1.31812, 0.0},    {0.25, 1.31812, 0.0},  {0.1, 1.82348, 0.0},
                              {0.1, 1.5708, 0.0},     {0.25, 1.82348, 0.0},  {0.25, 1.5708, 0.0},
                              {0.1, 2.0944, 0.0},     {0.25, 2.0944, 0.0},   {0.5, 2.0944, 0.0},
                              {0.1, 2.41886, 0.0},    {0.25, 2.41886, 0.0},  {0.5, 2.41886, 0.0},
                              {0.1, 3.14159, 0.0},    {0.25, 3.14159, 0.0},  {0.5, 3.14159, 0.0}};
    }
    if (gold_ghost_locations.size() != qindex.local_ghost_locations.size())
      ITFAILS;
    for (size_t i = 0; i < qindex.local_ghost_locations.size(); i++) {
      for (size_t d = 0; d < 2; d++)
        if (!rtt_dsxx::soft_equiv(gold_ghost_locations[i][d], qindex.local_ghost_locations[i][d],
                                  1e-4))
          ITFAILS;
    }

    // Check collect_ghost_data vector call
    std::vector<double> ghost_data(qindex.local_ghost_buffer_size, 0.0);
    qindex.collect_ghost_data(dd_data, ghost_data);
    std::vector<std::vector<double>> ghost_2x_data(
        2, std::vector<double>(qindex.local_ghost_buffer_size, 0.0));
    qindex.collect_ghost_data(dd_2x_data, ghost_2x_data);

    std::vector<double> gold_ghost_data;
    std::vector<std::vector<double>> gold_2x_ghost_data(2);
    if (rtt_c4::node() == 0) {
      gold_ghost_data = {3, 4, 5, 3, 4, 5, 3, 4, 5, 4, 5, 4, 4, 5, 5,
                         4, 5, 6, 4, 5, 6, 4, 5, 6, 6, 6, 6, 6, 6, 6};
      gold_2x_ghost_data[0] = {3, 4, 5, 3, 4, 5, 3, 4, 5, 4, 5, 4, 4, 5, 5,
                               4, 5, 6, 4, 5, 6, 4, 5, 6, 6, 6, 6, 6, 6, 6};
      gold_2x_ghost_data[1] = {9, 9, 9, 8, 8, 8, 7, 7, 7, 6, 6, 4, 5, 4, 5,
                               3, 3, 3, 2, 2, 2, 1, 1, 1, 9, 8, 7, 6, 4, 5};
    } else if (rtt_c4::node() == 1) {
      gold_ghost_data = {1, 2, 1, 2, 1, 2, 1, 2, 3, 1, 1, 2, 2, 3, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
                         6, 7, 8, 6, 7, 8, 6, 7, 8, 6, 7, 8, 6, 6, 7, 7, 8, 8, 7, 8, 7, 8, 7, 8};
      gold_2x_ghost_data[0] = {1, 2, 1, 2, 1, 2, 1, 2, 3, 1, 1, 2, 2, 3, 3, 1,
                               2, 3, 1, 2, 3, 1, 2, 3, 6, 7, 8, 6, 7, 8, 6, 7,
                               8, 6, 7, 8, 6, 6, 7, 7, 8, 8, 7, 8, 7, 8, 7, 8};
      gold_2x_ghost_data[1] = {9, 9, 8, 8, 7, 7, 6, 6, 6, 4, 5, 4, 5, 4, 5, 3,
                               3, 3, 2, 2, 2, 1, 1, 1, 9, 9, 9, 8, 8, 8, 7, 7,
                               7, 6, 6, 6, 4, 5, 4, 5, 4, 5, 3, 3, 2, 2, 1, 1};
    } else {
      gold_ghost_data = {1, 2, 1, 2, 1, 2, 1, 2, 3, 1, 1, 2, 2, 3, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
                         3, 4, 5, 3, 4, 5, 3, 4, 5, 4, 5, 4, 4, 5, 5, 4, 5, 6, 4, 5, 6, 4, 5, 6};
      gold_2x_ghost_data[0] = {1, 2, 1, 2, 1, 2, 1, 2, 3, 1, 1, 2, 2, 3, 3, 1,
                               2, 3, 1, 2, 3, 1, 2, 3, 3, 4, 5, 3, 4, 5, 3, 4,
                               5, 4, 5, 4, 4, 5, 5, 4, 5, 6, 4, 5, 6, 4, 5, 6};
      gold_2x_ghost_data[1] = {9, 9, 8, 8, 7, 7, 6, 6, 6, 4, 5, 4, 5, 4, 5, 3,
                               3, 3, 2, 2, 2, 1, 1, 1, 9, 9, 9, 8, 8, 8, 7, 7,
                               7, 6, 6, 4, 5, 4, 5, 3, 3, 3, 2, 2, 2, 1, 1, 1};
    }

    for (size_t i = 0; i < ghost_data.size(); i++)
      if (!rtt_dsxx::soft_equiv(ghost_data[i], gold_ghost_data[i]))
        ITFAILS;
    for (size_t i = 0; i < ghost_2x_data[0].size(); i++)
      if (!rtt_dsxx::soft_equiv(ghost_2x_data[0][i], gold_2x_ghost_data[0][i]))
        ITFAILS;
    for (size_t i = 0; i < ghost_2x_data[1].size(); i++)
      if (!rtt_dsxx::soft_equiv(ghost_2x_data[1][i], gold_2x_ghost_data[1][i]))
        ITFAILS;

    // check max sphere r window mapping (spoke shape more bins then data) functions
    {
      // wedge location +- 1 degree theta and 0.4 radius
      const std::array<double, 3> dr_dtheta{0.4, 0.0174533, 0.0};
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> center{qindex.locations[0]};
      const std::array<double, 3> win_min{center[0] - dr_dtheta[0], center[1] - dr_dtheta[1], 0.0};
      const std::array<double, 3> win_max{center[0] + dr_dtheta[0], center[1] + dr_dtheta[1], 0.0};
      const std::array<size_t, 3> bin_sizes{5, 1, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "max";
      std::vector<double> sphere_window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, sphere_window_data, win_min, win_max,
                                     bin_sizes, map_type, normalize, bias);
      std::vector<std::vector<double>> sphere_window_2x_data(2, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_2x_data, ghost_2x_data, sphere_window_2x_data, win_min,
                                     win_max, bin_sizes, map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_2x_data(2);
      // different result then 1D because the 1.0 y offset of the data
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 4.0, 5.0, 0.0};
        gold_window_2x_data[0] = {0.0, 0.0, 4.0, 5.0, 0.0};
        gold_window_2x_data[1] = {0.0, 0.0, 1.0, 1.0, 0.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {0.0, 0.0, 4.0, 5.0, 0.0};
        gold_window_2x_data[0] = {0.0, 0.0, 4.0, 5.0, 0.0};
        gold_window_2x_data[1] = {0.0, 0.0, 7.0, 7.0, 0.0};
      } else {
        gold_window_data = {5.0, 0.0, 6.0, 0.0, 7.0};
        gold_window_2x_data[0] = {5.0, 0.0, 6.0, 0.0, 7.0};
        gold_window_2x_data[1] = {4.0, 0.0, 4.0, 0.0, 4.0};
      }

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(sphere_window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 2; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(sphere_window_2x_data[v][i], gold_window_2x_data[v][i]))
            ITFAILS;
    }

    // check nearest sphere r window mapping (spoke shape more bins then data) functions
    {
      // wedge location +- 30 degree theta and 0.4 radius
      const std::array<double, 3> dr_dtheta{0.4, rtt_units::PI / 6.0, 0.0};
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> center{qindex.locations[0]};
      const std::array<double, 3> win_min{center[0] - dr_dtheta[0], center[1] - dr_dtheta[1], 0.0};
      const std::array<double, 3> win_max{center[0] + dr_dtheta[0], center[1] + dr_dtheta[1], 0.0};
      const std::array<size_t, 3> bin_sizes{5, 1, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "nearest";
      std::vector<double> sphere_window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, sphere_window_data, win_min, win_max,
                                     bin_sizes, map_type, normalize, bias);
      std::vector<std::vector<double>> sphere_window_2x_data(2, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_2x_data, ghost_2x_data, sphere_window_2x_data, win_min,
                                     win_max, bin_sizes, map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_2x_data(2);
      // different result then 1D because the 1.0 y offset of the data
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 1.0, 5.0, 0.0};
        gold_window_2x_data[0] = {0.0, 0.0, 1.0, 5.0, 0.0};
        gold_window_2x_data[1] = {0.0, 0.0, 1.0, 1.0, 0.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {0.0, 0.0, 3.0, 5.0, 0.0};
        gold_window_2x_data[0] = {0.0, 0.0, 3.0, 5.0, 0.0};
        gold_window_2x_data[1] = {0.0, 0.0, 7.0, 7.0, 0.0};
      } else {
        gold_window_data = {5.0, 0.0, 6.0, 0.0, 7.0};
        gold_window_2x_data[0] = {5.0, 0.0, 6.0, 0.0, 7.0};
        gold_window_2x_data[1] = {4.0, 0.0, 4.0, 0.0, 4.0};
      }

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(sphere_window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 2; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(sphere_window_2x_data[v][i], gold_window_2x_data[v][i]))
            ITFAILS;
    }

    // check max sphere r window mapping (shell shape more bins then data) functions
    {
      // wedge location +- 45 degree theta and 0.001 radius
      const std::array<double, 3> dr_dtheta{0.001, rtt_units::PI / 4.0, 0.0};
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> center{qindex.locations[0]};
      const std::array<double, 3> win_min{center[0] - dr_dtheta[0], center[1] - dr_dtheta[1], 0.0};
      const std::array<double, 3> win_max{center[0] + dr_dtheta[0], center[1] + dr_dtheta[1], 0.0};
      const std::array<size_t, 3> bin_sizes{1, 5, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "max";
      std::vector<double> sphere_window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, sphere_window_data, win_min, win_max,
                                     bin_sizes, map_type, normalize, bias);
      std::vector<std::vector<double>> sphere_window_2x_data(2, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_2x_data, ghost_2x_data, sphere_window_2x_data, win_min,
                                     win_max, bin_sizes, map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_2x_data(2);
      // different result then 1D because the 1.0 y offset of the data
      if (rtt_c4::node() == 0) {
        gold_window_data = {1.0, 0.0, 1.0, 0.0, 0.0};
        gold_window_2x_data[0] = {1.0, 0.0, 1.0, 0.0, 0.0};
        gold_window_2x_data[1] = {2.0, 0.0, 1.0, 0.0, 0.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {0.0, 3.0, 3.0, 3.0, 3.0};
        gold_window_2x_data[0] = {0.0, 3.0, 3.0, 3.0, 3.0};
        gold_window_2x_data[1] = {0.0, 8.0, 7.0, 6.0, 5.0};
      } else {
        gold_window_data = {6.0, 6.0, 6.0, 6.0, 6.0};
        gold_window_2x_data[0] = {6.0, 6.0, 6.0, 6.0, 6.0};
        gold_window_2x_data[1] = {7.0, 5.0, 4.0, 3.0, 2.0};
      }

      for (size_t i = 0; i < bin_sizes[1]; i++)
        if (!rtt_dsxx::soft_equiv(sphere_window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 2; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(sphere_window_2x_data[v][i], gold_window_2x_data[v][i]))
            ITFAILS;
    }

    // check max sphere r window mapping (shell shape more bins then data) functions with fill
    {
      // wedge location +- 22.5 degree theta and 0.001 radius
      const std::array<double, 3> dr_dtheta{0.001, rtt_units::PI / 4.0, 0.0};
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> center{qindex.locations[0]};
      const std::array<double, 3> win_min{center[0] - dr_dtheta[0], center[1] - dr_dtheta[1], 0.0};
      const std::array<double, 3> win_max{center[0] + dr_dtheta[0], center[1] + dr_dtheta[1], 0.0};
      const std::array<size_t, 3> bin_sizes{1, 5, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "max_fill";
      std::vector<double> sphere_window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, sphere_window_data, win_min, win_max,
                                     bin_sizes, map_type, normalize, bias);
      std::vector<std::vector<double>> sphere_window_2x_data(2, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_2x_data, ghost_2x_data, sphere_window_2x_data, win_min,
                                     win_max, bin_sizes, map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_2x_data(2);
      // different result then 1D because the 1.0 y offset of the data
      if (rtt_c4::node() == 0) {
        gold_window_data = {1.0, 1.0, 1.0, 1.0, 1.0};
        gold_window_2x_data[0] = {1.0, 1.0, 1.0, 1.0, 1.0};
        gold_window_2x_data[1] = {2.0, 2.0, 1.0, 1.0, 1.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {0.0, 3.0, 3.0, 3.0, 3.0};
        gold_window_2x_data[0] = {0.0, 3.0, 3.0, 3.0, 3.0};
        gold_window_2x_data[1] = {0.0, 8.0, 7.0, 6.0, 5.0};
      } else {
        gold_window_data = {6.0, 6.0, 6.0, 6.0, 6.0};
        gold_window_2x_data[0] = {6.0, 6.0, 6.0, 6.0, 6.0};
        gold_window_2x_data[1] = {7.0, 5.0, 4.0, 3.0, 2.0};
      }

      for (size_t i = 0; i < bin_sizes[1]; i++)
        if (!rtt_dsxx::soft_equiv(sphere_window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 2; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(sphere_window_2x_data[v][i], gold_window_2x_data[v][i]))
            ITFAILS;
    }

    // check min sphere r window mapping (shell shape more bins then data) functions with fill
    {
      // wedge location +- 22.5 degree theta and 0.001 radius
      const std::array<double, 3> dr_dtheta{0.001, rtt_units::PI / 4.0, 0.0};
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> center{qindex.locations[0]};
      const std::array<double, 3> win_min{center[0] - dr_dtheta[0], center[1] - dr_dtheta[1], 0.0};
      const std::array<double, 3> win_max{center[0] + dr_dtheta[0], center[1] + dr_dtheta[1], 0.0};
      const std::array<size_t, 3> bin_sizes{1, 5, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "min_fill";
      std::vector<double> sphere_window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, sphere_window_data, win_min, win_max,
                                     bin_sizes, map_type, normalize, bias);
      std::vector<std::vector<double>> sphere_window_2x_data(2, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_2x_data, ghost_2x_data, sphere_window_2x_data, win_min,
                                     win_max, bin_sizes, map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_2x_data(2);
      // different result then 1D because the 1.0 y offset of the data
      if (rtt_c4::node() == 0) {
        gold_window_data = {1.0, 1.0, 1.0, 1.0, 1.0};
        gold_window_2x_data[0] = {1.0, 1.0, 1.0, 1.0, 1.0};
        gold_window_2x_data[1] = {2.0, 2.0, 1.0, 1.0, 1.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {0.0, 3.0, 3.0, 3.0, 3.0};
        gold_window_2x_data[0] = {0.0, 3.0, 3.0, 3.0, 3.0};
        gold_window_2x_data[1] = {0.0, 8.0, 7.0, 6.0, 5.0};
      } else {
        gold_window_data = {6.0, 6.0, 6.0, 6.0, 6.0};
        gold_window_2x_data[0] = {6.0, 6.0, 6.0, 6.0, 6.0};
        gold_window_2x_data[1] = {6.0, 5.0, 4.0, 3.0, 2.0};
      }

      for (size_t i = 0; i < bin_sizes[1]; i++)
        if (!rtt_dsxx::soft_equiv(sphere_window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 2; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(sphere_window_2x_data[v][i], gold_window_2x_data[v][i]))
            ITFAILS;
    }

    // check nearest sphere r window mapping (shell shape more bins then data) functions with fill
    {
      // wedge location +- 45 degree theta and 0.001 radius
      const std::array<double, 3> dr_dtheta{0.001, rtt_units::PI / 4.0, 0.0};
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> center{qindex.locations[0]};
      const std::array<double, 3> win_min{center[0] - dr_dtheta[0], center[1] - dr_dtheta[1], 0.0};
      const std::array<double, 3> win_max{center[0] + dr_dtheta[0], center[1] + dr_dtheta[1], 0.0};
      const std::array<size_t, 3> bin_sizes{1, 5, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "nearest_fill";
      std::vector<double> sphere_window_data(5, 0.0);
      qindex.map_data_to_grid_window(dd_data, ghost_data, sphere_window_data, win_min, win_max,
                                     bin_sizes, map_type, normalize, bias);
      std::vector<std::vector<double>> sphere_window_2x_data(2, std::vector<double>(5, 0.0));
      qindex.map_data_to_grid_window(dd_2x_data, ghost_2x_data, sphere_window_2x_data, win_min,
                                     win_max, bin_sizes, map_type, normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_2x_data(2);
      // different result then 1D because the 1.0 y offset of the data
      if (rtt_c4::node() == 0) {
        gold_window_data = {1.0, 1.0, 1.0, 1.0, 1.0};
        gold_window_2x_data[0] = {1.0, 1.0, 1.0, 1.0, 1.0};
        gold_window_2x_data[1] = {2.0, 2.0, 1.0, 1.0, 1.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {0.0, 3.0, 3.0, 3.0, 3.0};
        gold_window_2x_data[0] = {0.0, 3.0, 3.0, 3.0, 3.0};
        gold_window_2x_data[1] = {0.0, 8.0, 7.0, 6.0, 5.0};
      } else {
        gold_window_data = {6.0, 6.0, 6.0, 6.0, 6.0};
        gold_window_2x_data[0] = {6.0, 6.0, 6.0, 6.0, 6.0};
        gold_window_2x_data[1] = {6.0, 5.0, 4.0, 3.0, 2.0};
      }

      for (size_t i = 0; i < bin_sizes[1]; i++)
        if (!rtt_dsxx::soft_equiv(sphere_window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 2; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(sphere_window_2x_data[v][i], gold_window_2x_data[v][i]))
            ITFAILS;
    }
  }

  if (ut.numFails == 0) {
    PASSMSG("quick_index sphere DD checks pass");
  } else {
    FAILMSG("quick_index sphere DD checks failed");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ParallelUnitTest ut(argc, argv, release);
  try {
    // >>> UNIT TESTS
    test_replication(ut);
    test_replication_sphere(ut);
    if (nodes() == 3) {
      test_decomposition(ut);
      test_decomposition_sphere(ut);
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstquick_index.cc
//------------------------------------------------------------------------------------------------//
