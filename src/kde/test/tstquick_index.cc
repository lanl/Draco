//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/test/tstquick_index.cc
 * \author Mathew Cleveland
 * \date   Aug. 10th 2021
 * \brief  quick_index testing function
 * \note   Copyright (C) 2021 Triad National Security, LLC., All rights reserved. 
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
  }

  if (ut.numFails == 0) {
    PASSMSG("quick_index checks pass");
  } else {
    FAILMSG("quick_index checks failed");
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

    // check max sphere r window mapping (more bins then data) functions
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> center{dd_position_array[0][0], dd_position_array[0][1], 0.0};
      // wedge location +- 1 degree theta and 0.5 radius
      const std::array<double, 3> dr_dtheta{0.4, 0.0174533, 0.0};
      const std::array<size_t, 3> bin_sizes{5, 1, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "max";
      const std::array<double, 3> sphere_center{-1.0, 0.0, 0.0};
      std::vector<double> sphere_window_data(5, 0.0);
      qindex.map_data_to_sphere_grid_window(dd_data, ghost_data, sphere_window_data, sphere_center,
                                            center, dr_dtheta, bin_sizes, map_type, normalize,
                                            bias);
      std::vector<std::vector<double>> sphere_window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_sphere_grid_window(dd_3x_data, ghost_3x_data, sphere_window_3x_data,
                                            sphere_center, center, dr_dtheta, bin_sizes, map_type,
                                            normalize, bias);
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
        if (!rtt_dsxx::soft_equiv(sphere_window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(sphere_window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }

    // check max sphere r window mapping (more bins then data) functions with bias
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> center{dd_position_array[0][0], dd_position_array[0][1], 0.0};
      // wedge location +- 1 degree theta and 0.5 radius
      const std::array<double, 3> dr_dtheta{0.4, 0.0174533, 0.0};
      const std::array<size_t, 3> bin_sizes{5, 1, 0};
      const bool normalize = false;
      const bool bias = true;
      const std::string map_type = "max";
      const std::array<double, 3> sphere_center{-1.0, 0.0, 0.0};
      std::vector<double> sphere_window_data(5, 0.0);
      qindex.map_data_to_sphere_grid_window(dd_data, ghost_data, sphere_window_data, sphere_center,
                                            center, dr_dtheta, bin_sizes, map_type, normalize,
                                            bias);
      std::vector<std::vector<double>> sphere_window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_sphere_grid_window(dd_3x_data, ghost_3x_data, sphere_window_3x_data,
                                            sphere_center, center, dr_dtheta, bin_sizes, map_type,
                                            normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      // different result then 1D because the 1.0 y offset of the data
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 3.0, 0.0, 0.0};
        gold_window_3x_data[0] = {0.0, 0.0, 3.0, 0.0, 0.0};
        gold_window_3x_data[1] = {0.0, 0.0, 4.0, 0.0, 0.0};
        gold_window_3x_data[2] = {0.0, 0.0, 0.0, 0.0, 0.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {0.0, 0.0, 6.0, 0.0, 0.0};
        gold_window_3x_data[0] = {0.0, 0.0, 6.0, 0.0, 0.0};
        gold_window_3x_data[1] = {0.0, 0.0, 7.0, 0.0, 0.0};
        gold_window_3x_data[2] = {0.0, 0.0, 0.0, 0.0, 0.0};
      } else {
        gold_window_data = {0.0, 0.0, 9.0, 0.0, 0.0};
        gold_window_3x_data[0] = {0.0, 0.0, 9.0, 0.0, 0.0};
        gold_window_3x_data[1] = {0.0, 0.0, 10.0, 0.0, 0.0};
        gold_window_3x_data[2] = {0.0, 0.0, 0.0, 0.0, 0.0};
      }

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(sphere_window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(sphere_window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }

    // check max sphere r window mapping (more bins then data) normalized
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> center{dd_position_array[0][0], dd_position_array[0][1], 0.0};
      // wedge location +- 1 degree theta and 0.5 radius
      const std::array<double, 3> dr_dtheta{0.4, 0.0174533, 0.0};
      const std::array<size_t, 3> bin_sizes{5, 1, 0};
      const bool normalize = true;
      const bool bias = false;
      const std::string map_type = "max";
      const std::array<double, 3> sphere_center{-1.0, 0.0, 0.0};
      std::vector<double> sphere_window_data(5, 0.0);
      qindex.map_data_to_sphere_grid_window(dd_data, ghost_data, sphere_window_data, sphere_center,
                                            center, dr_dtheta, bin_sizes, map_type, normalize,
                                            bias);
      std::vector<std::vector<double>> sphere_window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_sphere_grid_window(dd_3x_data, ghost_3x_data, sphere_window_3x_data,
                                            sphere_center, center, dr_dtheta, bin_sizes, map_type,
                                            normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      // different result then 1D because the 1.0 y offset of the data
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 1.0, 0.0, 0.0};
        gold_window_3x_data[0] = {0.0, 0.0, 1.0, 0.0, 0.0};
        gold_window_3x_data[1] = {0.0, 0.0, 1.0, 0.0, 0.0};
        gold_window_3x_data[2] = {0.0, 0.0, 1.0, 0.0, 0.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {0.0, 0.0, 1.0, 0.0, 0.0};
        gold_window_3x_data[0] = {0.0, 0.0, 1.0, 0.0, 0.0};
        gold_window_3x_data[1] = {0.0, 0.0, 1.0, 0.0, 0.0};
        gold_window_3x_data[2] = {0.0, 0.0, 1.0, 0.0, 0.0};
      } else {
        gold_window_data = {0.0, 0.0, 1.0, 0.0, 0.0};
        gold_window_3x_data[0] = {0.0, 0.0, 1.0, 0.0, 0.0};
        gold_window_3x_data[1] = {0.0, 0.0, 1.0, 0.0, 0.0};
        gold_window_3x_data[2] = {0.0, 0.0, 1.0, 0.0, 0.0};
      }

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(sphere_window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(sphere_window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }

    // check max sphere r window mapping (more bins then data) normalized bias fill
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> center{dd_position_array[0][0], dd_position_array[0][1], 0.0};
      // wedge location +- 1 degree theta and 0.5 radius
      const std::array<double, 3> dr_dtheta{0.4, 0.0174533, 0.0};
      const std::array<size_t, 3> bin_sizes{5, 1, 0};
      const bool normalize = true;
      const bool bias = true;
      const std::string map_type = "min_fill";
      const std::array<double, 3> sphere_center{-1.0, 0.0, 0.0};
      std::vector<double> sphere_window_data(5, 0.0);
      qindex.map_data_to_sphere_grid_window(dd_data, ghost_data, sphere_window_data, sphere_center,
                                            center, dr_dtheta, bin_sizes, map_type, normalize,
                                            bias);
      std::vector<std::vector<double>> sphere_window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_sphere_grid_window(dd_3x_data, ghost_3x_data, sphere_window_3x_data,
                                            sphere_center, center, dr_dtheta, bin_sizes, map_type,
                                            normalize, bias);
      std::vector<double> gold_window_data;
      std::vector<std::vector<double>> gold_window_3x_data(3);
      // different result then 1D because the 1.0 y offset of the data
      if (rtt_c4::node() == 0) {
        gold_window_data = {0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        gold_window_3x_data[0] = {0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        gold_window_3x_data[1] = {0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        gold_window_3x_data[2] = {0.0, 0.0, 0.0, 0.0, 0.0};
      } else if (rtt_c4::node() == 1) {
        gold_window_data = {0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        gold_window_3x_data[0] = {0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        gold_window_3x_data[1] = {0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        gold_window_3x_data[2] = {0.0, 0.0, 0.0, 0.0, 0.0};
      } else {
        gold_window_data = {0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        gold_window_3x_data[0] = {0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        gold_window_3x_data[1] = {0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        gold_window_3x_data[2] = {0.0, 0.0, 0.0, 0.0, 0.0};
      }

      for (size_t i = 0; i < bin_sizes[0]; i++)
        if (!rtt_dsxx::soft_equiv(sphere_window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(sphere_window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }

    // check max sphere r window mapping (more bins then data) functions
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> center{dd_position_array[0][0], dd_position_array[0][1], 0.0};
      // wedge location +- 1 degree theta and 0.5 radius
      const std::array<double, 3> dr_dtheta{0.4, 0.0174533, 0.0};
      const std::array<size_t, 3> bin_sizes{1, 5, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "ave";
      const std::array<double, 3> sphere_center{-1.0, 0.0, 0.0};
      std::vector<double> sphere_window_data(5, 0.0);
      qindex.map_data_to_sphere_grid_window(dd_data, ghost_data, sphere_window_data, sphere_center,
                                            center, dr_dtheta, bin_sizes, map_type, normalize,
                                            bias);
      std::vector<std::vector<double>> sphere_window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_sphere_grid_window(dd_3x_data, ghost_3x_data, sphere_window_3x_data,
                                            sphere_center, center, dr_dtheta, bin_sizes, map_type,
                                            normalize, bias);
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
        if (!rtt_dsxx::soft_equiv(sphere_window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(sphere_window_3x_data[v][i], gold_window_3x_data[v][i]))
            ITFAILS;
    }

    // check max sphere r window mapping (more bins then data) functions
    {
      // build a length=1.0 window around the first point on each node
      const std::array<double, 3> center{dd_position_array[0][0], dd_position_array[0][1], 0.0};
      // wedge location +- 1 degree theta and 0.5 radius
      const std::array<double, 3> dr_dtheta{0.4, 0.0174533, 0.0};
      const std::array<size_t, 3> bin_sizes{1, 5, 0};
      const bool normalize = false;
      const bool bias = false;
      const std::string map_type = "nearest";
      const std::array<double, 3> sphere_center{-1.0, 0.0, 0.0};
      std::vector<double> sphere_window_data(5, 0.0);
      qindex.map_data_to_sphere_grid_window(dd_data, ghost_data, sphere_window_data, sphere_center,
                                            center, dr_dtheta, bin_sizes, map_type, normalize,
                                            bias);
      std::vector<std::vector<double>> sphere_window_3x_data(3, std::vector<double>(5, 0.0));
      qindex.map_data_to_sphere_grid_window(dd_3x_data, ghost_3x_data, sphere_window_3x_data,
                                            sphere_center, center, dr_dtheta, bin_sizes, map_type,
                                            normalize, bias);
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
        if (!rtt_dsxx::soft_equiv(sphere_window_data[i], gold_window_data[i]))
          ITFAILS;
      for (size_t v = 0; v < 3; v++)
        for (size_t i = 0; i < bin_sizes[0]; i++)
          if (!rtt_dsxx::soft_equiv(sphere_window_3x_data[v][i], gold_window_3x_data[v][i]))
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
// end of tstquick_index.cc
//------------------------------------------------------------------------------------------------//
