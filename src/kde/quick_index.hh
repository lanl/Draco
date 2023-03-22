//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/quick_index.hh
 * \author Mathew Cleveland
 * \brief  This class generates coarse spatial indexing to quickly access near-neighbor data. This
 *         additionally provides simple interpolation schemes to map data to simple structured
 *         meshes.
 * \note   Copyright (C) 2021-2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_kde_quick_index_hh
#define rtt_kde_quick_index_hh

#include "c4/global.hh"
#include "units/MathConstants.hh"
#include <array>
#include <cmath>
#include <map>
#include <vector>

namespace rtt_kde {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Transform delta array (x, y, z) positions to (r, theta, phi) grid
 *
 * Calculate a relative r theta and phi coordinate relative to a sphere center location from a
 * standard (x,y,z) or (r,z) coordinates
 *
 * This is a bit hacky, we are defining simple min/max bounds in spherical coordinates based on the
 * defining cell corner nodes. This will create overlapping integration bounds, but we assume this
 * simple transformation is a better approximation then using point based integration. 
 *
 * \param[in] dim used to ensure it is only used in valid dimension ranges
 * \param[in] sphere_center center of sphere in (x,y,z) or (r,z) coordinates
 * \param[in] locations (x,y,z) or (r,z) locations to transform to relative (r, theta, phi) space.
 * \param[in] deltas defining the cell width in each dimension of a structured mesh
 */
inline std::vector<std::array<double, 3>>
transform_spherical_delta(const size_t dim, const std::array<double, 3> &sphere_center,
                          const std::vector<std::array<double, 3>> &locations,
                          const std::vector<std::array<double, 3>> &deltas) {
  Insist(dim == 2, "Transform_r_theta Only implemented in 2d");
  std::vector<std::array<double, 3>> r_theta_deltas(deltas);
  size_t i = 0;
  for (auto &delta : r_theta_deltas) {
    // define cell corners
    std::array<std::array<double, 3>, 4> cell_locations{};
    cell_locations[0] = std::array<double, 3>{locations[i][0] + delta[0] / 2.0,
                                              locations[i][1] + delta[1] / 2.0, 0.0};
    cell_locations[1] = std::array<double, 3>{locations[i][0] - delta[0] / 2.0,
                                              locations[i][1] + delta[1] / 2.0, 0.0};
    cell_locations[2] = std::array<double, 3>{locations[i][0] - delta[0] / 2.0,
                                              locations[i][1] - delta[1] / 2.0, 0.0};
    cell_locations[3] = std::array<double, 3>{locations[i][0] + delta[0] / 2.0,
                                              locations[i][1] - delta[1] / 2.0, 0.0};

    const std::array<double, 3> v_center{locations[i][0] - sphere_center[0],
                                         locations[i][1] - sphere_center[1], 0.0};
    const double r_center = sqrt(v_center[0] * v_center[0] + v_center[1] * v_center[1]);
    const double original_area = delta[0] * delta[1];
    // transform cell locations to spherical coordinates and find max
    std::array<double, 3> min{0.0, 0.0, 0.0};
    std::array<double, 3> max{0.0, 0.0, 0.0};
    size_t count = 0;
    for (auto &location : cell_locations) {
      const std::array<double, 3> v{location[0] - sphere_center[0], location[1] - sphere_center[1],
                                    0.0};
      const double r = sqrt(v[0] * v[0] + v[1] * v[1]);
      const double mag = sqrt(v[0] * v[0] + v[1] * v[1]);
      double cos_theta = mag > 0.0 ? std::max(std::min(v[1] / mag, 1.0), -1.0) : 0.0;
      location = std::array<double, 3>{
          r,
          location[0] < sphere_center[0] ? 2.0 * rtt_units::PI - acos(cos_theta) : acos(cos_theta),
          0.0};
      if (count == 0 || location[0] < min[0])
        min[0] = location[0];
      if (count == 0 || location[1] < min[1])
        min[1] = location[1];
      if (count == 0 || location[0] > max[0])
        max[0] = location[0];
      if (count == 0 || location[1] > max[1])
        max[1] = location[1];
      count++;
    }
    // update the delta based on min/max spherical bounding points
    delta[0] = max[0] - min[0];
    delta[1] = max[1] - min[1];
    // apply equal scaling (alpha) in each dimension for area preservation
    // new_A = (dtheta/2)*(rmax**2-rmin**2) ->
    // original_area = ((dtheta*alpha)/2)*((r0+(dr*alpha)/2)**2-(r0-(dr*alpha)/2)**2) ->
    // original_area = dtheta*r0*dr*alpha**2 ->
    // alpha = sqrt(original_area/(dtheta*r0*dr))
    const double alpha = std::sqrt(original_area / (delta[0] * r_center * delta[1]));
    delta[0] *= alpha;
    delta[1] *= alpha;
    i++;
  }

  return r_theta_deltas;
}
//================================================================================================//
//! Provide a hash like index of spatial distributed data along with simple mapping functions.
//================================================================================================//

class quick_index {
public:
  //! cartsian constructor
  quick_index(const size_t dim, const std::vector<std::array<double, 3>> &locations,
              const double max_window_size, const size_t bins_per_dimension,
              const bool domain_decomposed, const bool spherical = false,
              const std::array<double, 3> &sphere_center = {0.0, 0.0, 0.0});

  //! Structured constructor
  quick_index(const size_t dim, const std::vector<std::array<double, 3>> &locations,
              const std::vector<std::array<double, 3>> &deltas, const double max_window_size,
              const size_t bins_per_dimension, const bool domain_decomposed,
              const bool spherical = false,
              const std::array<double, 3> &sphere_center = {0.0, 0.0, 0.0});

  //! Collect Ghost Data
  void collect_ghost_data(const std::vector<double> &local_data,
                          std::vector<double> &local_ghost_data) const;

  //! Override function for integer vectors
  void collect_ghost_data(const std::vector<int> &local_data,
                          std::vector<int> &local_ghost_data) const;

  //! Override function of 3D array ghost data.
  void collect_ghost_data(const std::vector<std::array<double, 3>> &local_data,
                          std::vector<std::array<double, 3>> &local_ghost_data) const;

  //! Override function for vector<vector<double> array ghost data.
  void collect_ghost_data(const std::vector<std::vector<double>> &local_data,
                          std::vector<std::vector<double>> &local_ghost_data) const;

  //! Override function for vector<vector<double> array ghost data.
  void collect_ghost_data(const std::vector<std::vector<double>> &local_data,
                          std::vector<std::vector<double>> &local_ghost_data,
                          const size_t data_size) const;

  //! Fetch list of coarse index values bound by the window
  std::vector<size_t> window_coarse_index_list(const std::array<double, 3> &window_min,
                                               const std::array<double, 3> &window_max) const;

  //! Map local+ghost data to grid window
  void map_data_to_grid_window(const std::vector<double> &local_data,
                               const std::vector<double> &ghost_data,
                               std::vector<double> &grid_data,
                               const std::array<double, 3> &window_min,
                               const std::array<double, 3> &window_max,
                               const std::array<size_t, 3> &grid_bins, const std::string &map_type,
                               const bool normalize, const bool bias) const;

  //! Map local+ghost data to grid window for multi-dimensional data
  void map_data_to_grid_window(const std::vector<std::vector<double>> &local_data,
                               const std::vector<std::vector<double>> &ghost_data,
                               std::vector<std::vector<double>> &grid_data,
                               const std::array<double, 3> &window_min,
                               const std::array<double, 3> &window_max,
                               const std::array<size_t, 3> &grid_bins, const std::string &map_type,
                               const bool normalize, const bool bias) const;

  //! Calculate the orthogonal distance between to locations
  std::array<double, 3> calc_orthogonal_distance(const std::array<double, 3> &r0,
                                                 const std::array<double, 3> &r) const;

  // PUBLIC DATA
  // Quick index initialization data
  const size_t dim;
  const bool domain_decomposed;
  const bool spherical;
  const std::array<double, 3> sphere_center;
  const size_t coarse_bin_resolution;
  const std::vector<std::array<double, 3>> locations;
  const std::vector<std::array<double, 3>> deltas;
  const size_t n_locations;

  // Global bounds
  std::array<double, 3> bounding_box_min{0.0, 0.0, 0.0};
  std::array<double, 3> bounding_box_max{0.0, 0.0, 0.0};
  // Local Data map
  std::map<size_t, std::vector<size_t>> coarse_index_map{};
  std::map<size_t, std::array<double, 3>> coarse_index_center{};
  std::map<size_t, std::array<double, 3>> coarse_index_size{};

  // DOMAIN DECOMPOSED DATA
  // Local bounds
  std::array<double, 3> local_bounding_box_min{0.0, 0.0, 0.0};
  std::array<double, 3> local_bounding_box_max{0.0, 0.0, 0.0};
  // Ordered list of local bins (indexes values are based on the global bin structure)
  std::vector<size_t> local_bins{};
  // Size of ghost data buffer
  size_t local_ghost_buffer_size{0};
  // Map used to index into a local ghost buffer
  std::map<size_t, std::vector<size_t>> local_ghost_index_map{};
  // Local ghost locations (build at construction time)
  std::vector<std::array<double, 3>> local_ghost_locations{};
  // Local ghost deltas (build at construction time)
  std::vector<std::array<double, 3>> local_ghost_deltas{};

private:
  //! Data initialization function
  void initialize_data();

  // PRIVATE DATA
  // Map used to write local data to other processor ghost cells
  // put_window_map[global_id] = [put_rank, ghost_proc_buffer_size, ghost_proc_put_offset]
  // array is integers to accommodate mpi data types
  std::map<size_t, std::vector<std::array<int, 2>>> put_window_map{};
  // max put buffer size;
  size_t max_put_buffer_size;
  double max_window_size;
};

} // end namespace  rtt_kde

#endif // rtt_kde_quick_index_hh

//------------------------------------------------------------------------------------------------//
// end of kde/quick_index.hh
//------------------------------------------------------------------------------------------------//
