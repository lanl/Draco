//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/quick_index.cc
 * \author Mathew Cleveland
 * \brief  Explicitly defined quick_index functions.
 * \note   Copyright (C) 2021-2021 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "quick_index.hh"
#include "ds++/dbc.hh"
#include <numeric>
#include <tuple>

namespace rtt_kde {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief quick_index constructor. 
 *
 * This function builds up a global indexing table to quickly access data that is spatial located
 * near each other. It breaks up the data into equally spaced bins in each dimension. For domain
 * decomposed data it builds a one sided communication map to place local data that is need on other
 * processors for ghost cells. The ghost cell extents is determined by the max_data_window spatial
 * size such that any cell on the local domain will have access to all points that should fall into
 * the spatial window centered on any given local point.
 *
 * \param[in] dim_ specifying the data dimensionality 
 * \param[in] locations_ data locations.
 * \param[in] max_window_size_ maximum supported window size
 * \param[in] bins_per_dimension_ number of bins in each dimension
 * \param[in] domain_decomposed_
 * \param[in] spherical_ bool operator to enable spherical transform
 * \param[in] sphere_center_ origin of spherical transform
 */
quick_index::quick_index(const size_t dim_, const std::vector<std::array<double, 3>> &locations_,
                         const double max_window_size_, const size_t bins_per_dimension_,
                         const bool domain_decomposed_, const bool spherical_,
                         const std::array<double, 3> &sphere_center_)
    : dim(dim_), domain_decomposed(domain_decomposed_), spherical(spherical_),
      sphere_center(sphere_center_), coarse_bin_resolution(bins_per_dimension_),
      max_window_size(max_window_size_),
      locations(spherical ? transform_spherical(dim_, sphere_center_, locations_) : locations_),
      n_locations(locations_.size()) {
  Require(dim > 0);
  Require(coarse_bin_resolution > 0);

  // Build local bounding box
  bounding_box_min = {0, 0, 0};
  bounding_box_max = {0, 0, 0};
  // only set initial values for working dimensions
  for (size_t d = 0; d < dim; d++) {
    bounding_box_min[d] = 1e20;
    bounding_box_max[d] = -1e20;
  }

  for (auto &loc : locations) {
    for (size_t d = 0; d < dim; d++) {
      if (loc[d] < bounding_box_min[d])
        bounding_box_min[d] = loc[d];
      if (loc[d] > bounding_box_max[d])
        bounding_box_max[d] = loc[d];
    }
  }

  if (domain_decomposed) {
    // Store the local bounding box and extend to maximum non-local data size
    local_bounding_box_min = bounding_box_min;
    local_bounding_box_max = bounding_box_max;
    for (size_t d = 0; d < dim; d++) {
      double wsize = max_window_size * 0.5;
      if (spherical && d == 1) {
        // Transform to dtheta via arch_lenght=r*dtheta
        // enforce a 90 degree maximum angle
        wsize = std::min(rtt_units::PI / 2, 0.5 * max_window_size / local_bounding_box_max[0]);
      }
      local_bounding_box_min[d] -= wsize;
      local_bounding_box_max[d] += wsize;
      // No negative radius values
      if (spherical && d == 0)
        local_bounding_box_min[d] = std::max(0.0, local_bounding_box_min[d]);
    }
    // Global reduce to get the global min and max
    rtt_c4::global_min(&bounding_box_min[0], 3);
    rtt_c4::global_max(&bounding_box_max[0], 3);
    if (!spherical) {
      // spherical theta bounds can exceed global bounds because the window wraps around theta=0.
      for (size_t d = 0; d < dim; d++) {
        local_bounding_box_min[d] = std::max(local_bounding_box_min[d], bounding_box_min[d]);
        local_bounding_box_max[d] = std::min(local_bounding_box_max[d], bounding_box_max[d]);
      }
    }
  }

  // temp cast corse_bin_resolution to double for interpolation
  const auto crd = static_cast<double>(coarse_bin_resolution);

  // build up the local hash table of into global bins
  size_t locIndex = 0;
  for (auto &loc : locations) {
    std::array<size_t, 3> index{0UL, 0UL, 0UL};
    for (size_t d = 0; d < dim; d++) {
      Check(bounding_box_min[d] < bounding_box_max[d]);
      index[d] = static_cast<size_t>(std::floor(crd * (loc[d] - bounding_box_min[d]) /
                                                (bounding_box_max[d] - bounding_box_min[d])));
      index[d] = std::min(index[d], coarse_bin_resolution - 1);
    }
    // build up the local index hash
    const size_t global_index = index[0] + index[1] * coarse_bin_resolution +
                                index[2] * coarse_bin_resolution * coarse_bin_resolution;
    coarse_index_map[global_index].push_back(locIndex);
    locIndex++;
  }

  // Now we need to build up ghost location map data for domain decomposed mode
  if (domain_decomposed) {
    // temporary cast of the nodes to prevent conversion warnings
    const auto nodes = static_cast<size_t>(rtt_c4::nodes());
    const auto node = static_cast<size_t>(rtt_c4::node());

    // build list of local bins based on the local bounds
    local_bins = window_coarse_index_list(local_bounding_box_min, local_bounding_box_max);

    // build a global map for number of entries into the global bins on each processor
    // creates a (nbins**dim)*nranks sized array
    // NOTE: If this gets to big we could stride over a subset of coarse bins
    // and do multiple iterations of mpi communication to build up the map
    size_t nbins = coarse_bin_resolution;
    for (size_t d = 1; d < dim; d++)
      nbins *= coarse_bin_resolution;

    std::vector<int> global_index_per_bin_per_proc(nbins * nodes, 0UL);
    for (auto &map : coarse_index_map) {
      size_t gipbpp_index = map.first + nbins * node;
      // must cast to an int to accomidate mpi int types.
      global_index_per_bin_per_proc[gipbpp_index] = static_cast<int>(map.second.size());
    }
    rtt_c4::global_sum(&global_index_per_bin_per_proc[0], nbins * nodes);

    // calculate local ghost buffer size
    local_ghost_buffer_size = 0;
    for (size_t proc = 0; proc < nodes; proc++) {
      for (auto &bin : local_bins) {
        if (node != proc) {
          size_t gipbpp_index = bin + nbins * proc;
          // build up the local ghost index map
          for (int i = 0; i < global_index_per_bin_per_proc[gipbpp_index]; i++)
            local_ghost_index_map[bin].push_back(local_ghost_buffer_size + i);
          // accumulate the total ghost points
          local_ghost_buffer_size += global_index_per_bin_per_proc[gipbpp_index];
        }
      }
    }

    std::vector<int> global_need_bins_per_proc(nbins * nodes, 0UL);
    // global need bins
    for (auto &bin : local_bins) {
      global_need_bins_per_proc[bin + nbins * node] += 1;
    }
    rtt_c4::global_sum(&global_need_bins_per_proc[0], nbins * nodes);

    // calculate the put map so each node knows which processor to send data
    // and where to index that data
    // PERFORMANCE NOTE: This would be more efficient to use a MPI_SCAN and
    // std::partial_sum but I need to think how this would actually look.
    max_put_buffer_size = 0;
    for (int rec_proc = 0; rec_proc < rtt_c4::nodes(); rec_proc++) {
      // calculating the offset SUCKS!!! If anyone can find a better way please help.
      int offset = 0;
      for (int send_proc = 0; send_proc < rtt_c4::node(); send_proc++) {
        if (rec_proc == send_proc)
          continue;
        for (size_t bin = 0; bin < nbins; bin++) {
          if (global_need_bins_per_proc[bin + nbins * rec_proc] > 0) {
            offset += global_index_per_bin_per_proc[bin + nbins * send_proc];
          }
        }
      }
      for (auto &map : coarse_index_map) {
        if (rtt_c4::node() != rec_proc) {
          size_t gipbpp_index = map.first + nbins * rec_proc;
          if (global_need_bins_per_proc[gipbpp_index] > 0) {
            // capture the largest put buffer on this rank
            if (map.second.size() > max_put_buffer_size)
              max_put_buffer_size = map.second.size();

            // build up map data
            put_window_map[map.first].push_back(std::array<int, 2>{rec_proc, offset});
            offset += static_cast<int>(map.second.size());
          }
        }
      }
    }

    // allocate ghost locations
    local_ghost_locations =
        std::vector<std::array<double, 3>>(local_ghost_buffer_size, {0.0, 0.0, 0.0});
    // collect the local ghost locations
    collect_ghost_data(locations, local_ghost_locations);

  } // End domain decomposed data construction
}

#ifdef C4_MPI
//------------------------------------------------------------------------------------------------//
// call MPI_put using a chunk style write to avoid error in MPI_put with large local buffers.
auto put_lambda = [](auto &put, auto &put_buffer, auto &put_size, auto &win) {
  // temporary work around until RMA is available in c4
  // loop over all ranks we need to send this buffer too.
  for (auto &putv : put.second) {
    const int put_rank = putv[0];
    const int put_offset = putv[1];
    // This is dumb, but we need to write in chunks because MPI_Put writes
    // junk with large (>10,000) buffer sizes.
    int chunk_size = 1000;
    const auto nchunks = static_cast<int>(
        std::ceil(static_cast<double>(put_size) / static_cast<double>(chunk_size)));
    int nput = 0;
    for (int c = 0; c < nchunks; c++) {
      chunk_size = std::min(chunk_size, static_cast<int>(put_size) - nput);
      Check(chunk_size > 0);
      MPI_Put(&put_buffer[nput], chunk_size, MPI_DOUBLE, put_rank, put_offset, chunk_size,
              MPI_DOUBLE, win);
      nput += chunk_size;
    }
  }
};
#endif

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Collect ghost data for a vector<std::array<double, 3>>
 * 
 * Collect ghost data for vector of 3 dimensional arrays. This function uses RMA and the local
 * put_window_map to allow each rank to independently fill in its data to ghost cells of other
 * ranks.
 *
 * \param[in] local_data the local 3 dimensional data that is required to be available as ghost cell
 * data on other processors.
 * \param[in] local_ghost_data the resulting 3 dimensional ghost data data. 
 */
void quick_index::collect_ghost_data(const std::vector<std::array<double, 3>> &local_data,
                                     std::vector<std::array<double, 3>> &local_ghost_data) const {
  Require(local_data.size() == n_locations);
  Insist(domain_decomposed, "Calling collect_ghost_data with a quick_index object that specified "
                            "domain_decomposed=.false.");

  Insist(local_ghost_data.size() == local_ghost_buffer_size,
         "ghost_data input must be sized via quick_index.local_ghost_buffer_size");
#ifdef C4_MPI // temporary work around until RMA is available in c4
  // Use one sided MPI Put commands to fill up the ghost cell location data
  std::vector<double> local_ghost_buffer(local_ghost_buffer_size, 0.0);
  std::vector<double> put_buffer(max_put_buffer_size, 0.0);
  MPI_Win win;
  MPI_Win_create(local_ghost_buffer.data(), local_ghost_buffer_size * sizeof(double),
                 sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

  // working from my local data put the ghost data on the other ranks
  for (size_t d = 0; d < dim; d++) {
    Remember(int errorcode =) MPI_Win_fence(MPI_MODE_NOSTORE, win);
    Check(errorcode == MPI_SUCCESS);
    for (auto &put : put_window_map) {
      // use map.at() to allow const access
      Check((coarse_index_map.at(put.first)).size() <= max_put_buffer_size);
      // fill up the current ghost cell data for this dimension
      int putIndex = 0;
      for (auto &l : coarse_index_map.at(put.first)) {
        put_buffer[putIndex] = local_data[l][d];
        putIndex++;
      }
      put_lambda(put, put_buffer, putIndex, win);
    }
    Remember(errorcode =) MPI_Win_fence((MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED), win);
    Check(errorcode == MPI_SUCCESS);

    // alright move the position buffer to the final correct array positions
    int posIndex = 0;
    for (auto &pos : local_ghost_buffer) {
      local_ghost_data[posIndex][d] = pos;
      posIndex++;
    }
  }
  MPI_Win_free(&win);
#endif
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Collect ghost data for a vector<vector<double>> 
 * 
 * Collect ghost data for vector<vector<double>> arrays. This function uses RMA and the local
 * put_window_map to allow each rank to independently fill in its data to ghost cells of other
 * ranks.
 *
 * \param[in] local_data the local multi-dimensional data that is required to be available as ghost
 * cell data on other processors.
 * \param[in,out] local_ghost_data the resulting multi-dimensional ghost data
 */
void quick_index::collect_ghost_data(const std::vector<std::vector<double>> &local_data,
                                     std::vector<std::vector<double>> &local_ghost_data) const {
  Insist(domain_decomposed, "Calling collect_ghost_data with a quick_index object that specified "
                            "domain_decomposed=.false.");
  size_t data_dim = local_data.size();
  size_t ghost_data_dim = local_ghost_data.size();
  Insist(data_dim == ghost_data_dim,
         "The local_data.size() and the local_ghost_data.size() vectors much match");
  // Check ghost data
  for (size_t d = 0; d < ghost_data_dim; d++) {
    Insist(local_ghost_data[d].size() == local_ghost_buffer_size,
           "ghost_data[" + std::to_string(d) +
               "] input must be sized via quick_index.local_ghost_buffer_size");
  }
#ifdef C4_MPI // temporary work around until RMA is available in c4
  // Use one sided MPI Put commands to fill up the ghost cell location data
  std::vector<double> local_ghost_buffer(local_ghost_buffer_size, 0.0);
  std::vector<double> put_buffer(max_put_buffer_size, 0.0);
  MPI_Win win;
  MPI_Win_create(local_ghost_buffer.data(), local_ghost_buffer_size * sizeof(double),
                 sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

  // working from my local data put the ghost data on the other ranks
  for (size_t d = 0; d < data_dim; d++) {
    Check(local_data[d].size() == n_locations);
    Remember(int errorcode =) MPI_Win_fence(MPI_MODE_NOSTORE, win);
    Check(errorcode == MPI_SUCCESS);
    for (auto &put : put_window_map) {
      // use map.at() to allow const access
      Check((coarse_index_map.at(put.first)).size() <= max_put_buffer_size);
      // fill up the current ghost cell data for this dimension
      int putIndex = 0;
      for (auto &l : coarse_index_map.at(put.first)) {
        put_buffer[putIndex] = local_data[d][l];
        putIndex++;
      }
      put_lambda(put, put_buffer, putIndex, win);
    }
    Remember(errorcode =) MPI_Win_fence((MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED), win);
    Check(errorcode == MPI_SUCCESS);
    // alright move the position buffer to the final correct vector positions
    int posIndex = 0;
    for (auto &pos : local_ghost_buffer) {
      local_ghost_data[d][posIndex] = pos;
      posIndex++;
    }
  }
  MPI_Win_free(&win);
#endif
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Collect ghost data for a vector<double>
 * 
 * Collect ghost data for a single vector. This function uses RMA and the local put_window_map to
 * allow each rank to independently fill in its data to ghost cells of other ranks.
 *
 * \param[in] local_data the local vector data that is required to be available as ghost cell data
 * on other processors.
 * \param[in,out] local_ghost_data the resulting ghost data
 */
void quick_index::collect_ghost_data(const std::vector<double> &local_data,
                                     std::vector<double> &local_ghost_data) const {
  Require(local_data.size() == n_locations);
  Insist(domain_decomposed, "Calling collect_ghost_data with a quick_index object that specified "
                            "domain_decomposed=.false.");
  Insist(local_ghost_data.size() == local_ghost_buffer_size,
         "ghost_data input must be sized via quick_index.local_ghost_buffer_size");
#ifdef C4_MPI // temporary work around until RMA is available in c4
  std::vector<double> local_ghost_buffer(local_ghost_buffer_size, 0.0);
  std::vector<double> put_buffer(max_put_buffer_size, 0.0);
  MPI_Win win;
  MPI_Win_create(local_ghost_data.data(), local_ghost_buffer_size * sizeof(double), sizeof(double),
                 MPI_INFO_NULL, MPI_COMM_WORLD, &win);

  // working from my local data put the ghost data on the other ranks
  Remember(int errorcode =) MPI_Win_fence(MPI_MODE_NOSTORE, win);
  Check(errorcode == MPI_SUCCESS);
  for (auto put : put_window_map) {
    // use map.at() to allow const access
    Check((coarse_index_map.at(put.first)).size() <= max_put_buffer_size);
    // fill up the current ghost cell data for this dimension
    int putIndex = 0;
    for (auto &l : coarse_index_map.at(put.first)) {
      put_buffer[putIndex] = local_data[l];
      putIndex++;
    }
    put_lambda(put, put_buffer, putIndex, win);
  }
  Remember(errorcode =) MPI_Win_fence((MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED), win);
  Check(errorcode == MPI_SUCCESS);
  MPI_Win_free(&win);
#endif
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Generate a coarse index list for a window.
 * 
 *  Provides a list of global indices that are required by any given window range.
 *
 * \param[in] window_min the smallest corner point for every dimension
 * \param[in] window_max the largest corner point for every dimension
 * \return bin_list list of global bins requested for the current window.
 */
std::vector<size_t>
quick_index::window_coarse_index_list(const std::array<double, 3> &window_min,
                                      const std::array<double, 3> &window_max) const {
  Require(window_min[0] <= window_max[0]);
  Require(window_min[1] <= window_max[1]);
  Require(window_min[2] <= window_max[2]);

  // temp cast corse_bin_resolution to double for interpolation
  const auto crd = static_cast<double>(coarse_bin_resolution);

  // calculate the global index range that each processor needs to
  // accommodate the specified data window size
  std::array<size_t, 3> index_min = {0UL, 0UL, 0UL};
  std::array<size_t, 3> index_max = {0UL, 0UL, 0UL};
  size_t nbins = 1;
  for (size_t d = 0; d < dim; d++) {
    // because local bounds can extend beyond the mesh we need to force a
    // positive index if necessary
    double wmin = window_min[d];
    if (spherical && d == 1 && window_min[d] < bounding_box_min[d])
      wmin = bounding_box_min[d]; // truncate to standard theta space
    double wmax = window_max[d];
    if (spherical && d == 1 && window_max[d] > bounding_box_max[d])
      wmax = bounding_box_max[d]; // truncate to standard theta space
    index_min[d] = static_cast<size_t>(std::floor(std::max(
        crd * (wmin - bounding_box_min[d]) / (bounding_box_max[d] - bounding_box_min[d]), 0.0)));
    index_max[d] = static_cast<size_t>(std::floor(crd * (wmax - bounding_box_min[d]) /
                                                  (bounding_box_max[d] - bounding_box_min[d])));
    // because local bounds can extend beyond the mesh we need to floor to
    // the max bin size
    index_max[d] = std::min(index_max[d], coarse_bin_resolution - 1);

    // Use multiplicity to accumulate total bins;
    if ((index_max[d] - index_min[d]) > 0)
      nbins *= index_max[d] - index_min[d] + 1;
  }

  // Fill up bin list
  size_t count = 0;
  std::vector<size_t> bin_list(nbins);
  for (size_t k = index_min[2]; k <= index_max[2]; k++) {
    for (size_t j = index_min[1]; j <= index_max[1]; j++) {
      for (size_t i = index_min[0]; i <= index_max[0]; i++) {
        size_t bin_index =
            i + j * coarse_bin_resolution + k * coarse_bin_resolution * coarse_bin_resolution;
        bin_list[count] = bin_index;
        count++;
      }
    }
  }
  // Fill in the overflow around theta=0.0
  if (spherical && (window_min[1] < 0.0 || window_max[1] > 2.0 * rtt_units::PI)) {
    // Only one bound of the window should every overshoot zero
    Check(!(window_min[1] < 0.0 && window_max[1] > 2.0 * rtt_units::PI));
    size_t overlap_nbins = 1;
    for (size_t d = 0; d < dim; d++) {
      // because local bounds can extend beyond the mesh we need to force a
      // positive index if necessary
      double wmin = window_min[d];
      double wmax = window_max[d];
      if (spherical && d == 1 && window_min[d] < 0.0) {
        // Capture the overshoot theta space
        wmin = std::min(2.0 * rtt_units::PI + window_min[d], bounding_box_max[d]);
        wmax = 2.0 * rtt_units::PI;
      }
      if (spherical && d == 1 && window_max[d] > 2.0 * rtt_units::PI) {
        // Capture the overshoot theta space
        wmin = 0.0;
        wmax = std::max(window_max[d] - 2.0 * rtt_units::PI, bounding_box_min[d]);
      }
      // Truncate based on global space
      if (spherical && d == 1 && wmin < bounding_box_min[d])
        wmin = bounding_box_min[d]; // truncate to standard theta space
      if (spherical && d == 1 && wmax > bounding_box_max[d])
        wmax = bounding_box_max[d]; //trunacte to standard theta space
      index_min[d] = static_cast<size_t>(std::floor(std::max(
          crd * (wmin - bounding_box_min[d]) / (bounding_box_max[d] - bounding_box_min[d]), 0.0)));
      index_max[d] = static_cast<size_t>(std::floor(std::max(
          crd * (wmax - bounding_box_min[d]) / (bounding_box_max[d] - bounding_box_min[d]), 0.0)));
      // because local bounds can extend beyond the mesh we need to floor to
      // the max bin size
      index_min[d] = std::min(index_min[d], coarse_bin_resolution - 1);
      index_max[d] = std::min(index_max[d], coarse_bin_resolution - 1);

      // Use multiplicity to accumulate total bins;
      if ((index_max[d] - index_min[d]) > 0)
        overlap_nbins *= index_max[d] - index_min[d] + 1;
    }
    for (size_t k = index_min[2]; k <= index_max[2]; k++) {
      for (size_t j = index_min[1]; j <= index_max[1]; j++) {
        for (size_t i = index_min[0]; i <= index_max[0]; i++) {
          size_t bin_index =
              i + j * coarse_bin_resolution + k * coarse_bin_resolution * coarse_bin_resolution;
          // make sure we don't duplicate a bin here
          if (std::find(bin_list.begin(), bin_list.end(), bin_index) == bin_list.end()) {
            bin_list.push_back(bin_index);
          }
        }
      }
    }
  }

  return bin_list;
}

//------------------------------------------------------------------------------------------------//
// Lambda for getting the mapped window bin
auto get_window_bin = [](auto spherical, const auto dim, const auto &grid_bins,
                         const auto &location, const auto &window_min, const auto &window_max,
                         const auto &Remember(n_map_bins)) {
  // calculate local bin index
  bool valid = true;
  std::array<size_t, 3> bin_id{0, 0, 0};
  double distance_to_bin_center = 0.0;
  for (size_t d = 0; d < dim; d++) {
    Check((window_max[d] - window_min[d]) > 0.0);
    double loc = location[d];
    // transform location for zero theta overshoot
    if (spherical && d == 1 && window_max[d] > 2.0 * rtt_units::PI &&
        location[d] < (window_max[d] - 2.0 * rtt_units::PI))
      loc += 2.0 * rtt_units::PI;
    // transform location for zero theta overshoot
    if (spherical && d == 1 && window_min[d] < 0 &&
        location[d] > (2.0 * rtt_units::PI + window_min[d]))
      loc -= 2.0 * rtt_units::PI;
    const double bin_value =
        static_cast<double>(grid_bins[d]) * (loc - window_min[d]) / (window_max[d] - window_min[d]);
    if (bin_value < 0.0 || bin_value > static_cast<double>(grid_bins[d])) {
      valid = false;
      break;
    } else {
      bin_id[d] = static_cast<size_t>(bin_value);
      // catch any values exactly on the edge of the top bin
      bin_id[d] = std::min(grid_bins[d] - 1, bin_id[d]);
      const double bin_center =
          window_min[d] + (static_cast<double>(bin_id[d]) / static_cast<double>(grid_bins[d]) +
                           0.5 / static_cast<double>(grid_bins[d])) *
                              (window_max[d] - window_min[d]);
      // approximate in spherical geometry;
      distance_to_bin_center += (bin_center - loc) * (bin_center - loc);
    }
  }
  distance_to_bin_center =
      rtt_dsxx::soft_equiv(distance_to_bin_center, 0.0) ? 0.0 : sqrt(distance_to_bin_center);
  const size_t local_window_bin =
      bin_id[0] + bin_id[1] * grid_bins[0] + bin_id[2] * grid_bins[0] * grid_bins[1];

  Check(valid ? local_window_bin < n_map_bins : true);

  return std::tuple<bool, size_t, double>{valid, local_window_bin, distance_to_bin_center};
};

//------------------------------------------------------------------------------------------------//
// Lambda for mapping the data
auto map_data = [](auto &bias_cell_count, auto &data_count, auto &grid_data, auto &min_distance,
                   const auto &map_type, const auto &data, const auto &distance_to_bin_center,
                   const auto &local_window_bin, const auto &data_bin) {
  // regardless of map type if it is the first value to enter the bin it
  // gets set to that value
  if (data_count[local_window_bin] == 0) {
    bias_cell_count += 1.0;
    data_count[local_window_bin]++;
    min_distance[local_window_bin] = distance_to_bin_center;
    grid_data[local_window_bin] = data[data_bin];
  } else if (map_type == "max") {
    if (data[data_bin] > grid_data[local_window_bin])
      grid_data[local_window_bin] = data[data_bin];
  } else if (map_type == "min") {
    if (data[data_bin] < grid_data[local_window_bin])
      grid_data[local_window_bin] = data[data_bin];
  } else if (map_type == "ave") {
    data_count[local_window_bin] += 1;
    grid_data[local_window_bin] += data[data_bin];
  } else if (map_type == "nearest") {
    if (rtt_dsxx::soft_equiv(distance_to_bin_center, min_distance[local_window_bin])) {
      data_count[local_window_bin] += 1;
      grid_data[local_window_bin] += data[data_bin];
    } else if (distance_to_bin_center < min_distance[local_window_bin]) {
      min_distance[local_window_bin] = distance_to_bin_center;
      data_count[local_window_bin] = 1;
      grid_data[local_window_bin] = data[data_bin];
    } // else exclude the far points.
  } else {
    Insist(false, "Error: map_type=" + map_type + " is invalid. Must be max, min, or ave.");
  }
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Map data to a grid window for vector<double> data
 * 
 * Maps local+ghost data to a fixed mesh grid based on a specified weighting type. This data can
 * additionally be normalized and positively biased on the grid.
 * 
 *
 * \param[in] local_data the local data on the processor to be mapped to the window
 * \param[in] ghost_data the ghost data on the processor to be mapped to the window
 * \param[in,out] grid_data the resulting data map
 * \param[in] window_min the smallest corner point for every dimension
 * \param[in] window_max the largest corner point for every dimension
 * \param[in] grid_bins number of equally spaced bins in each dir
 * \param[in] map_type_in string indicating the mapping (max, min, ave)
 * \param[in] normalize bool operator to specify if the data should be normalized to a pdf
 * \param[in] bias bool operator to specify if the data should be moved to the
 * positive domain space
 */
void quick_index::map_data_to_grid_window(
    const std::vector<double> &local_data, const std::vector<double> &ghost_data,
    std::vector<double> &grid_data, const std::array<double, 3> &window_min,
    const std::array<double, 3> &window_max, const std::array<size_t, 3> &grid_bins,
    const std::string &map_type_in, const bool normalize, const bool bias) const {
  Require(local_data.size() == n_locations);
  Require(!(window_max[0] < window_min[0]));
  Require(!(window_max[1] < window_min[1]));
  Require(!(window_max[2] < window_min[2]));
  Require(domain_decomposed ? ghost_data.size() == local_ghost_buffer_size : true);
  Require(domain_decomposed
              ? (fabs(window_max[0] - window_min[0]) - max_window_size) / max_window_size < 1e-6
              : true);
  Remember(double ymax = spherical ? std::min(rtt_units::PI / 2.0, max_window_size / window_max[0])
                                   : max_window_size);
  Require(domain_decomposed ? (fabs(window_max[1] - window_min[1]) - ymax) / ymax < 1e-6 : true);
  Require(domain_decomposed
              ? (fabs(window_max[2] - window_min[2]) - max_window_size) / max_window_size < 1e-6
              : true);

  bool fill = false;
  std::string map_type = map_type_in;
  if (map_type_in == "max_fill") {
    Insist((grid_bins[0] > 1 && grid_bins[1] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[1] > 1 && grid_bins[0] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[2] > 1 && grid_bins[0] <= 1 && grid_bins[1] <= 1),
           "one of grid bins must be == 1, Grid must be 1D to use max_fill option");
    fill = true;
    map_type = "max";
  } else if (map_type_in == "min_fill") {
    Insist((grid_bins[0] > 1 && grid_bins[1] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[1] > 1 && grid_bins[0] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[2] > 1 && grid_bins[0] <= 1 && grid_bins[1] <= 1),
           "one of grid bins must be == 1, Grid must be 1D to use min_fill option");
    fill = true;
    map_type = "min";
  } else if (map_type_in == "ave_fill") {
    Insist((grid_bins[0] > 1 && grid_bins[1] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[1] > 1 && grid_bins[0] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[2] > 1 && grid_bins[0] <= 1 && grid_bins[1] <= 1),
           "one of grid bins must be == 1, Grid must be 1D to use ave_fill option");
    fill = true;
    map_type = "ave";
  } else if (map_type_in == "nearest_fill") {
    Insist((grid_bins[0] > 1 && grid_bins[1] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[1] > 1 && grid_bins[0] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[2] > 1 && grid_bins[0] <= 1 && grid_bins[1] <= 1),
           "one of grid bins must be == 1, Grid must be 1D to use ave_fill option");
    fill = true;
    map_type = "nearest";
  }

  for (size_t d = 0; d < dim; d++)
    Insist(grid_bins[d] > 0, "Bin size must be greater then zero for each active dimension");

  size_t n_map_bins = 1;
  for (size_t d = 0; d < dim; d++)
    n_map_bins *= grid_bins[d];

  Insist(grid_data.size() == n_map_bins,
         "grid_data must match the flatten grid_bin size for the active dimensions (in 3d "
         "grid_data.size()==grib_bins[0]*grid_bins[1]*grid_bins[2])");

  std::fill(grid_data.begin(), grid_data.end(), 0.0);

  // Grab the global bins that lie in this window
  std::vector<size_t> global_bins = window_coarse_index_list(window_min, window_max);

  std::vector<int> data_count(n_map_bins, 0);
  std::vector<double> min_distance(n_map_bins, 0);
  double bias_cell_count = 0.0;
  // Loop over all possible bins
  for (auto &cb : global_bins) {
    // skip bins that aren't present in the map (can't use [] operator with constness)
    // loop over the local data
    auto mapItr = coarse_index_map.find(cb);
    if (mapItr != coarse_index_map.end()) {
      for (auto &l : mapItr->second) {
        bool valid;
        size_t local_window_bin;
        double distance_to_bin_center;
        std::tie(valid, local_window_bin, distance_to_bin_center) = get_window_bin(
            spherical, dim, grid_bins, locations[l], window_min, window_max, n_map_bins);

        // If the bin is outside the window continue to the next poin
        if (!valid)
          continue;

        // lambda for mapping the data
        map_data(bias_cell_count, data_count, grid_data, min_distance, map_type, local_data,
                 distance_to_bin_center, local_window_bin, l);

      } // end local point loop
    }   // if valid local bin loop
    if (domain_decomposed) {
      // loop over the ghost data
      auto gmapItr = local_ghost_index_map.find(cb);
      if (gmapItr != local_ghost_index_map.end()) {
        // loop over ghost data
        for (auto &g : gmapItr->second) {
          bool valid;
          size_t local_window_bin;
          double distance_to_bin_center;
          std::tie(valid, local_window_bin, distance_to_bin_center) =
              get_window_bin(spherical, dim, grid_bins, local_ghost_locations[g], window_min,
                             window_max, n_map_bins);

          // If the bin is outside the window continue to the next poin
          if (!valid)
            continue;

          // lambda for mapping the data
          map_data(bias_cell_count, data_count, grid_data, min_distance, map_type, ghost_data,
                   distance_to_bin_center, local_window_bin, g);
        } // end ghost point loop
      }   // if valid ghost bin
    }     // if dd
  }       // end coarse bin loop

  if (map_type == "ave" || map_type == "nearest") {
    for (size_t i = 0; i < n_map_bins; i++) {
      if (data_count[i] > 0) {
        grid_data[i] /= data_count[i];
      }
    }
  }
  if (fill) {
    double last_val = 0.0;
    int last_data_count = 0;
    for (size_t i = 0; i < n_map_bins; i++) {
      if (data_count[i] > 0) {
        last_val = grid_data[i];
        last_data_count = data_count[i];
      } else {
        grid_data[i] = last_val;
        data_count[i] = last_data_count;
      }
    }
  }

  if (bias && normalize) {
    // return a positive normalized distribution
    const double bias_value =
        fabs(std::min(0.0, *std::min_element(grid_data.begin(), grid_data.end())));
    const double sum =
        std::accumulate(grid_data.begin(), grid_data.end(), 0.0) + bias_value * bias_cell_count;
    // catch zero instance
    const double scale = !rtt_dsxx::soft_equiv(sum, 0.0) ? 1.0 / sum : 1.0;
    for (size_t i = 0; i < n_map_bins; i++)
      grid_data[i] = (grid_data[i] + bias_value) * scale;
  } else if (bias) {
    // return a positive distribution
    const double bias_value =
        fabs(std::min(0.0, *std::min_element(grid_data.begin(), grid_data.end())));
    for (size_t i = 0; i < n_map_bins; i++)
      grid_data[i] += bias_value;
  } else if (normalize) {
    // return a normalized distribution
    const double sum = std::accumulate(grid_data.begin(), grid_data.end(), 0.0);
    // catch zero instance
    const double scale = !rtt_dsxx::soft_equiv(sum, 0.0) ? 1.0 / sum : 1.0;
    for (size_t i = 0; i < n_map_bins; i++)
      grid_data[i] *= scale;
  }
}

//------------------------------------------------------------------------------------------------//
// Lambda for mapping the vector data
auto map_vector_data = [](auto &bias_cell_count, auto &data_count, auto &grid_data,
                          auto &min_distance, const auto &map_type, const auto &data,
                          const auto &distance_to_bin_center, const auto &local_window_bin,
                          const auto &data_bin, const auto &vsize) {
  // regardless of map type if it is the first value to enter the bin it gets set to that value
  if (data_count[local_window_bin] == 0) {
    bias_cell_count += 1.0;
    data_count[local_window_bin]++;
    min_distance[local_window_bin] = distance_to_bin_center;
    for (size_t v = 0; v < vsize; v++)
      grid_data[v][local_window_bin] = data[v][data_bin];
  } else if (map_type == "max") {
    for (size_t v = 0; v < vsize; v++)
      if (data[v][data_bin] > grid_data[v][local_window_bin])
        grid_data[v][local_window_bin] = data[v][data_bin];
  } else if (map_type == "min") {
    for (size_t v = 0; v < vsize; v++)
      if (data[v][data_bin] < grid_data[v][local_window_bin])
        grid_data[v][local_window_bin] = data[v][data_bin];
  } else if (map_type == "ave") {
    data_count[local_window_bin] += 1;
    for (size_t v = 0; v < vsize; v++)
      grid_data[v][local_window_bin] += data[v][data_bin];
  } else if (map_type == "nearest") {
    if (rtt_dsxx::soft_equiv(distance_to_bin_center, min_distance[local_window_bin])) {
      data_count[local_window_bin] += 1;
      for (size_t v = 0; v < vsize; v++)
        grid_data[v][local_window_bin] += data[v][data_bin];
    } else if (distance_to_bin_center < min_distance[local_window_bin]) {
      min_distance[local_window_bin] = distance_to_bin_center;
      data_count[local_window_bin] = 1;
      for (size_t v = 0; v < vsize; v++)
        grid_data[v][local_window_bin] = data[v][data_bin];
    } // else exclude the far points.
  } else {
    Insist(false, "Error: map_type=" + map_type + " is invalid. Must be max, min, or ave.");
  }
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Map data to a grid window for vector<vector<double>>
 * 
 * Maps multiple local+ghost data vectors to a fixed mesh grid based on a specified weighting type.
 * This data can additionally be normalized and positively biased on the grid.
 * 
 *
 * \param[in] local_data the local data on the processor to be mapped to the window
 * \param[in] ghost_data the ghost data on the processor to be mapped to the window
 * \param[in,out] grid_data the resulting data map
 * \param[in] window_min the smallest corner point for every dimension
 * \param[in] window_max the largest corner point for every dimension
 * \param[in] grid_bins number of equally spaced bins in each dir
 * \param[in] map_type_in string indicating the mapping (max, min, ave)
 * \param[in] normalize bool operator to specify if the data should be normalized to a pdf
 * (independent of each data vector)
 * \param[in] bias bool operator to specify if the data should be moved to the positive domain space
 * (independent of each data vector)
 * \return bin_list list of global bins requested for the current window.
 */
void quick_index::map_data_to_grid_window(const std::vector<std::vector<double>> &local_data,
                                          const std::vector<std::vector<double>> &ghost_data,
                                          std::vector<std::vector<double>> &grid_data,
                                          const std::array<double, 3> &window_min,
                                          const std::array<double, 3> &window_max,
                                          const std::array<size_t, 3> &grid_bins,
                                          const std::string &map_type_in, const bool normalize,
                                          const bool bias) const {
  Require(domain_decomposed ? local_data.size() == ghost_data.size() : true);
  Require(!(window_max[0] < window_min[0]));
  Require(!(window_max[1] < window_min[1]));
  Require(!(window_max[2] < window_min[2]));
  Require(domain_decomposed
              ? (fabs(window_max[0] - window_min[0]) - max_window_size) / max_window_size < 1e-6
              : true);
  Remember(double ymax = spherical ? std::min(rtt_units::PI / 2.0, max_window_size / window_max[0])
                                   : max_window_size;)
      Require(domain_decomposed ? (fabs(window_max[1] - window_min[1]) - ymax) / ymax < 1e-6
                                : true);
  Require(domain_decomposed
              ? (fabs(window_max[2] - window_min[2]) - max_window_size) / max_window_size < 1e-6
              : true);

  bool fill = false;
  std::string map_type = map_type_in;
  if (map_type_in == "max_fill") {
    Insist((grid_bins[0] > 1 && grid_bins[1] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[1] > 1 && grid_bins[0] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[2] > 1 && grid_bins[0] <= 1 && grid_bins[1] <= 1),
           "one of grid bins must be == 1, Grid must be 1D to use max_fill option");
    fill = true;
    map_type = "max";
  } else if (map_type_in == "min_fill") {
    Insist((grid_bins[0] > 1 && grid_bins[1] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[1] > 1 && grid_bins[0] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[2] > 1 && grid_bins[0] <= 1 && grid_bins[1] <= 1),
           "one of grid bins must be == 1, Grid must be 1D to use min_fill option");
    fill = true;
    map_type = "min";
  } else if (map_type_in == "ave_fill") {
    Insist((grid_bins[0] > 1 && grid_bins[1] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[1] > 1 && grid_bins[0] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[2] > 1 && grid_bins[0] <= 1 && grid_bins[1] <= 1),
           "one of grid bins must be == 1, Grid must be 1D to use ave_fill option");
    fill = true;
    map_type = "ave";
  } else if (map_type_in == "nearest_fill") {
    Insist((grid_bins[0] > 1 && grid_bins[1] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[1] > 1 && grid_bins[0] <= 1 && grid_bins[2] <= 1) ||
               (grid_bins[2] > 1 && grid_bins[0] <= 1 && grid_bins[1] <= 1),
           "one of grid bins must be == 1, Grid must be 1D to use ave_fill option");
    fill = true;
    map_type = "nearest";
  }

  for (size_t d = 0; d < dim; d++)
    Insist(grid_bins[d] > 0, "Bin size must be greater then zero for each active dimension");

  const size_t vsize = local_data.size();
  // Grab the global bins that lie in this window
  std::vector<size_t> global_bins = window_coarse_index_list(window_min, window_max);
  size_t n_map_bins = 1;
  for (size_t d = 0; d < dim; d++) {
    n_map_bins *= grid_bins[d];
  }

  for (size_t v = 0; v < vsize; v++) {
    Insist(grid_data[v].size() == n_map_bins,
           "grid_data[" + std::to_string(v) +
               "] must match the flatten grid_bin size for the active dimensions (in 3d "
               "grid_data.size()==grib_bins[0]*grid_bins[1]*grid_bins[2])");
    std::fill(grid_data[v].begin(), grid_data[v].end(), 0.0);
  }

  // initialize grid data
  std::vector<int> data_count(n_map_bins, 0);
  std::vector<double> min_distance(n_map_bins, 0);
  double bias_cell_count = 0.0;
  // Loop over all possible bins
  for (auto &cb : global_bins) {
    // skip bins that aren't present in the map (can't use [] operator with constness)
    // loop over the local data
    auto mapItr = coarse_index_map.find(cb);
    if (mapItr != coarse_index_map.end()) {
      for (auto &l : mapItr->second) {
        bool valid;
        size_t local_window_bin;
        double distance_to_bin_center;
        std::tie(valid, local_window_bin, distance_to_bin_center) = get_window_bin(
            spherical, dim, grid_bins, locations[l], window_min, window_max, n_map_bins);
        // If the bin is outside the window continue to the next poin
        if (!valid)
          continue;
        Check(local_window_bin < n_map_bins);
        map_vector_data(bias_cell_count, data_count, grid_data, min_distance, map_type, local_data,
                        distance_to_bin_center, local_window_bin, l, vsize);
      } // end local point loop
    }   // if valid local bin loop
    if (domain_decomposed) {
      // loop over the ghost data
      auto gmapItr = local_ghost_index_map.find(cb);
      if (gmapItr != local_ghost_index_map.end()) {
        // loop over ghost data
        for (auto &g : gmapItr->second) {
          bool valid;
          size_t local_window_bin;
          double distance_to_bin_center;
          std::tie(valid, local_window_bin, distance_to_bin_center) =
              get_window_bin(spherical, dim, grid_bins, local_ghost_locations[g], window_min,
                             window_max, n_map_bins);

          // If the bin is outside the window continue to the next poin
          if (!valid)
            continue;
          map_vector_data(bias_cell_count, data_count, grid_data, min_distance, map_type,
                          ghost_data, distance_to_bin_center, local_window_bin, g, vsize);
        } // end ghost point loop
      }   // if valid ghost bin
    }     // if dd
  }       // end coarse bin loop

  if (map_type == "ave" || map_type == "nearest") {
    for (size_t i = 0; i < n_map_bins; i++) {
      for (size_t v = 0; v < vsize; v++) {
        if (data_count[i] > 0) {
          grid_data[v][i] /= data_count[i];
        }
      }
    }
  }
  if (fill) {
    std::vector<double> last_val(vsize, 0.0);
    int last_data_count = 0;
    for (size_t i = 0; i < n_map_bins; i++) {
      for (size_t v = 0; v < vsize; v++) {
        if (data_count[i] > 0) {
          last_val[v] = grid_data[v][i];
          last_data_count = data_count[i];
        } else {
          grid_data[v][i] = last_val[v];
          if (v == vsize - 1)
            data_count[i] = last_data_count;
        }
      }
    }
  }

  if (bias && normalize) {
    // return a positive normalized distribution
    for (size_t v = 0; v < vsize; v++) {
      const double bias_value =
          fabs(std::min(0.0, *std::min_element(grid_data[v].begin(), grid_data[v].end())));
      const double sum = std::accumulate(grid_data[v].begin(), grid_data[v].end(), 0.0) +
                         bias_value * bias_cell_count;
      // catch zero instance
      const double scale = !rtt_dsxx::soft_equiv(sum, 0.0) ? 1.0 / sum : 1.0;
      for (size_t i = 0; i < n_map_bins; i++)
        if (data_count[i] > 0)
          grid_data[v][i] = (grid_data[v][i] + bias_value) * scale;
    }
  } else if (bias) {
    // return a positive distribution
    for (size_t v = 0; v < vsize; v++) {
      const double bias_value =
          fabs(std::min(0.0, *std::min_element(grid_data[v].begin(), grid_data[v].end())));
      for (size_t i = 0; i < n_map_bins; i++)
        if (data_count[i] > 0)
          grid_data[v][i] += bias_value;
    }
  } else if (normalize) {
    // return a normalized distribution
    for (size_t v = 0; v < vsize; v++) {
      const double sum = std::accumulate(grid_data[v].begin(), grid_data[v].end(), 0.0);
      // catch zero instance
      const double scale = !rtt_dsxx::soft_equiv(sum, 0.0) ? 1.0 / sum : 1.0;
      for (size_t i = 0; i < n_map_bins; i++)
        if (data_count[i] > 0)
          grid_data[v][i] *= scale;
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Calculate the orthogonal distance between two points
 * 
 * Maps multiple local+ghost data vectors to a fixed mesh grid based on a specified weighting type.
 * This data can additionally be normalized and positively biased on the grid.
 * 
 *
 * \param[in] r0 initial position
 * \param[in] r final position
 * \param[in] arch_radius this is used for spherical geometry to determine at what radial point the
 * orthognal archlength is measured. This point must be bound by the initial and final point radius.
 * \return orthognal distance between the initial and final point in each direction
 */
std::array<double, 3> quick_index::calc_orthogonal_distance(const std::array<double, 3> &r0,
                                                            const std::array<double, 3> &r,
                                                            const double arch_radius) const {
  Require(spherical ? dim == 2 : true);
  Require(spherical ? !(arch_radius < 0.0) : true);
  return {r[0] - r0[0], spherical ? arch_radius * (r[1] - r0[1]) : r[1] - r0[1], r[2] - r0[2]};
}

} // namespace rtt_kde

//------------------------------------------------------------------------------------------------//
// end of quick_index.cc
//------------------------------------------------------------------------------------------------//
