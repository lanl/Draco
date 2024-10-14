//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/CAMR_Mesh.hh
 * \author H. Park
 * \date   May 10 2023
 * \brief  CAMR_Mesh class definition.
 * \note   Copyright (C) 2023-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_mesh_CAMR_Mesh_hh
#define rtt_mesh_CAMR_Mesh_hh

#include "Comm_Patterns.hh"
#include "Parallel_Data_Operator.hh"
#include "Topology.hh"
#include "Topology_Builder.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "c4/C4_Functions.hh"
#include "mesh/CAMR_Layout.hh"
#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <sstream>

// cSpell:ignore topo, nfine, LCOV

namespace rtt_mesh {

//================================================================================================//
/*!
 * \brief Represents a cell-based AMR mesh
 *
 * Class is templated on the dimension of the mesh. Only values of 2 or 3 are allowed.
 *
 * The vertex data is assumed to be stored in x,y[,z] order, then by cell.
 *
 * This class imposes an extra constraint on the data ordering of neighbors in the
 * rtt_mesh::CAMR_Layout object. The neighbors of a cell should be listed in order of *increasing*
 * coordinate values, with the increasing coordinates in the order x, y [,z]. For example, face one
 * of a cell is in the yz plane in a 3D mesh.  The four neighbors in a refined cell should be listed
 * in the order shown below:
 *
 *\code
 *                     ^ z
 *                     |
 *                     +-----+-----+
 *                     |     |     |
 *                     | 3   | 4   |
 *                     +-----+-----+
 *                     |     |     |
 *                     | 1   | 2   |
 *                     +-----+-----+---> y
 *\endcode
 *
 * We note that the neighboring cells do *not* need to *numbered* in this order. They only need to
 * be *listed as neighbors* in this order.
 *
 * For neighbor cells which are also part of the mesh, we can verify this constraint with the
 * neighbor cell's geometric information. We cannot do this for off-processor neighbors, (without
 * parallel communication) and so the convention is implicitly assumed.
 *
 * \tparam D Dimensionality, one of {1, 2, 3}.
 *
 * \sa CAMR_Mesh.cc for detailed descriptions.
 */
//================================================================================================//

template <int D> class CAMR_Mesh {
public:
  // REQUIRED TYPEDEFS
  using SP_Comm_Patterns = std::shared_ptr<Comm_Patterns>;
  using SP_Topology = std::shared_ptr<Topology>;

  using sf_int = std::vector<int>;
  using vf_int = std::vector<std::vector<int>>;
  using sf_double = std::vector<double>;
  using vf_double = std::vector<std::vector<double>>;

  //replacing rtt_tracking::U3.
  using U3 = std::array<double, 3>;

  //! Default constructors.
  template <typename I_IT, typename D_IT>
  CAMR_Mesh(size_t const size_, double const width_, D_IT vertex_data, I_IT level_data,
            I_IT neighbor_count, I_IT neighbor_values, bool const submesh_ = false);

  // this construct stores ghost info.
  template <typename I_IT, typename D_IT>
  CAMR_Mesh(int64_t const num_global_cells_, size_t const active_size_, size_t const full_size_,
            size_t const ghost_size_, double const width_, D_IT vertex_data, I_IT level_data,
            I_IT neighbor_count, I_IT neighbor_values, I_IT ghost_cells, I_IT ghost_procs,
            I_IT full2active, I_IT active2full, bool const submesh_ = false);

  //  template <typename I_IT, typename D_IT>
  //  CAMR_Mesh(double const width_, D_IT vertex_data, I_IT level_data,
  //            rtt_mesh::CAMR_Layout const &layout_, bool const submesh_ = false);

  //! Disabled copy/move constructor
  CAMR_Mesh(CAMR_Mesh const &rhs) = delete;
  CAMR_Mesh(CAMR_Mesh const &&rhs) = delete;

  //! Disable assignment
  CAMR_Mesh &operator=(CAMR_Mesh const &rhs) = delete;
  CAMR_Mesh &operator=(CAMR_Mesh const &&rhs) = delete;

  ~CAMR_Mesh() = default;

  // ACCESSORS
  NODISCARD sf_double const &get_widths() const { return widths; }
  NODISCARD size_t get_level(size_t cell) const;
  NODISCARD size_t get_maximum_level() const { return maximum_level; }
  NODISCARD size_t get_size() const { return active_size; }
  NODISCARD size_t num_cells() const { return active_size; }
  NODISCARD size_t get_num_global_cells() const { return num_global_cells; }
  NODISCARD size_t get_num_full_cells() const { return full_size; }
  NODISCARD size_t get_num_ghost_cells() const { return ghost_size; }

  NODISCARD double get_width(size_t level = 0) const;
  NODISCARD double get_cell_width(size_t cell) const;

  NODISCARD sf_double get_vertex(size_t cell) const;

  NODISCARD sf_int get_bnd_cells() const { return bnd_cells; }
  NODISCARD sf_int get_bnd_procs() const { return bnd_procs; }
  NODISCARD sf_int get_active2full() const { return active2full; }
  NODISCARD sf_int get_full2active() const { return full2active; }

  NODISCARD inline bool is_refined(int const cell, uint32_t const face) const;
  NODISCARD inline int is_coarse(int const cell, uint32_t const face) const;

  NODISCARD inline sf_int next_cells(int const cell, int const coarse_face) const;
  NODISCARD inline sf_int next_cells(size_t const cell, size_t const coarse_face) const;

  NODISCARD int64_t global_cell(size_t const cell) const { return topology_->global_cell(cell); }
  NODISCARD int64_t get_global_offset(int rank) const { return topology_->get_global_offset(rank); }
  NODISCARD int remote_proc(size_t boundary_cell) const {
    return topology_->remote_proc(boundary_cell);
  }

  //! Comparison operator
  bool operator==(CAMR_Mesh<D> const &rhs) const;

  //! Negative comparison operator
  bool operator!=(CAMR_Mesh<D> const &rhs) const { return !(*this == rhs); }

  NODISCARD rtt_mesh::CAMR_Layout const &get_Layout() const { return *layout; }

  //! Return the number of finite elements per cell (placeholder).
  static constexpr size_t get_num_fe() {
    return static_cast<uint32_t>(2) << static_cast<uint32_t>(D - 1);
  }

  //--------------------//
  // APPLICATION INTERFACE

  // Mesh-wide information:

  NODISCARD int get_spatial_dimension() const { return D; }
  NODISCARD bool full_Mesh() const { return !submesh; } //!< Full or partitioned mesh?
  NODISCARD double get_total_volume() const { return local_volume; }

  // Cell-by-cell information. Geometric and topological.

  NODISCARD double min(size_t d, size_t cell) const;
  NODISCARD double max(size_t d, size_t cell) const;

  NODISCARD inline double cell_width(size_t const coordinate, size_t const cell) const;
  NODISCARD bool in_cell(size_t const cell, U3 const &position) const;
  NODISCARD bool hard_in_cell(size_t const cell, U3 const &position) const;
  NODISCARD bool in_cell(size_t const cell, const uint32_t face, U3 const &position) const;
  NODISCARD int next_cell(int const cell, uint32_t const face, U3 const &position) const;
  NODISCARD inline size_t next_face(int const cell, const uint32_t coarse_face) const;
  NODISCARD inline size_t next_face(size_t const cell, const size_t coarse_face) const;
  NODISCARD int get_cell(U3 const &position) const;
  NODISCARD inline uint32_t get_face(U3 const &position, const size_t cell) const;
  NODISCARD inline U3 get_face_center(size_t const cell, const uint32_t face) const;
  NODISCARD inline U3 get_cell_centroid(size_t const cell) const;
  NODISCARD U3 get_shared_face_center(size_t const cell, uint32_t const face,
                                      int const next_cell) const;
  NODISCARD double get_db(U3 const &position, U3 const &direction, int const cell,
                          uint32_t &face) const;
  NODISCARD double get_orthogonal_dist_to_bnd(U3 const &position, int const cell) const;
  NODISCARD U3 get_normal(int const cell, uint32_t const face) const;
  NODISCARD U3 get_normal_in(size_t const cell, uint32_t const face) const;
  NODISCARD double face_area(size_t const cell, const size_t face) const;
  NODISCARD double mean_chord_length(int const cell) const;
  NODISCARD double get_min_width(size_t cell) const;
  // this returns dx*dy or dr*dz, not true volume for rz.
  NODISCARD double volume(size_t const cell) const;
  NODISCARD vf_double get_vertices(int const cell) const;
  NODISCARD vf_double get_vertices(int const cell, const uint32_t face) const;
  NODISCARD inline sf_double get_short_vertices(size_t cell) const;
  NODISCARD inline U3 get_intersection(const U3 &r, const U3 &omega, const double &distance,
                                       const int &cell, const uint32_t face) const;
  NODISCARD sf_int get_neighbors(uint32_t const cell) const;
  NODISCARD size_t get_number_of_faces() const { return layout->num_coarse_faces(); }
  NODISCARD inline sf_int get_valid_cell_faces() const;

  // Mesh boundary information:

  NODISCARD size_t get_bndface(std::string const &boundary_name, int const /*cell*/) const;
  NODISCARD sf_int get_surcells(std::string const &boundary_name) const;
  NODISCARD bool check_defined_surcells(std::string const &boundary_name,
                                        sf_int const &cell_list) const;

  // END APPLICATION INTERFACE
  //------------------------//

  //-----------------------------------//
  // Methods required for graphics dumps.
  NODISCARD sf_int get_cell_types() const;
  NODISCARD vf_double get_point_coord() const;
  NODISCARD vf_int get_cell_pair() const;

  NODISCARD const sf_double &get_vertices() const { return vertices; }
  NODISCARD const sf_int &get_levels() const { return levels; }
  //-----------------------------------//

private:
  // IMPLEMENTATION

  void build_widths();
  void compute_total_volume();

  template <typename I_IT, typename D_IT> void constructor_body(I_IT level_data, D_IT vertex_data);
  template <typename I_IT, typename D_IT>
  void constructor_body(I_IT level_data, D_IT vertex_data, I_IT neighbor_count);

  NODISCARD inline sf_double::const_iterator vertex_data(size_t cell) const;

  // DATA

  std::shared_ptr<rtt_mesh::CAMR_Layout> layout;

  //
  sf_double vertices;
  sf_int levels;
  sf_int neighbor_count;

  int64_t num_global_cells;
  int decomposition = rtt_mesh::Topology::DOMAIN_REPLICATION;

  size_t maximum_level;
  size_t number_of_levels;
  size_t active_size;
  size_t full_size;
  size_t ghost_size;
  size_t vertices_size;

  double width;        //<! Constant width of zero level cells.
  double local_volume; //<! Volume of all cells in this mesh.
  sf_double widths;    //<! Widths of all cell levels.
  bool submesh;        //<! True if this is a DD submesh.

  sf_int bnd_cells;
  sf_int bnd_procs;
  sf_int active2full;
  sf_int full2active;

  SP_Topology topology_;
};

//------------------------------------------------------------------------------------------------//
// INLINE FUNCTIONS
//------------------------------------------------------------------------------------------------//
//replacing rtt_tracking::U3.
using U3 = std::array<double, 3>;

//------------------------------------------------------------------------------------------------//
//! Private function which handles the operations common to all constructors.
template <int D>
template <typename I_IT, typename D_IT>
void CAMR_Mesh<D>::constructor_body(I_IT level_data, D_IT vertex_data) {

  // Copy the argument data into member data:
  levels.assign(level_data, level_data + active_size);
  vertices.assign(vertex_data, vertex_data + vertices_size);

  build_widths();

  compute_total_volume();

  // this create MPI communication pattern.
  Topology_Builder topo_builder(decomposition, num_global_cells, bnd_cells.size(), bnd_cells.data(),
                                bnd_procs.data(), full_size, active_size, full2active.data());

  SP_Topology const topo_ = topo_builder.build_Topology();
  topology_ = topo_;

  SP_Comm_Patterns comm_pat_;
  // Make the communications patterns
  comm_pat_ = std::make_shared<rtt_mesh::Comm_Patterns>();
  comm_pat_->calc_patterns(topo_);

  // to identify ghost level and AMR configuration, we need to exchange level and vertex.
  std::shared_ptr<rtt_mesh::Parallel_Data_Operator> par_op_;
  par_op_ = std::make_shared<rtt_mesh::Parallel_Data_Operator>(topo_);
  par_op_->fill_ghost_cell_data(comm_pat_, levels);

  par_op_->fill_ghost_cell_data(comm_pat_, vertices, D);
}

//------------------------------------------------------------------------------------------------//
//! Private function which handles the operations common to all constructors.
template <int D>
template <typename I_IT, typename D_IT>
void CAMR_Mesh<D>::constructor_body(I_IT level_data, D_IT vertex_data, I_IT neighbor_count_data) {

  // Copy the argument data into member data:
  levels.assign(level_data, level_data + active_size);
  vertices.assign(vertex_data, vertex_data + vertices_size);
  neighbor_count.assign(neighbor_count_data, neighbor_count_data + active_size * 2 * D);

  build_widths();

  compute_total_volume();

  // this create MPI communication pattern.
  Topology_Builder topo_builder(decomposition, num_global_cells,
                                static_cast<uint32_t>(bnd_cells.size()), bnd_cells.data(),
                                bnd_procs.data(), static_cast<uint32_t>(full_size),
                                static_cast<uint32_t>(active_size), full2active.data());

  SP_Topology const topo_ = topo_builder.build_Topology();
  topology_ = topo_;

  SP_Comm_Patterns comm_pat_;
  // Make the communications patterns
  comm_pat_ = std::make_shared<rtt_mesh::Comm_Patterns>();
  comm_pat_->calc_patterns(topo_);

  // to identify ghost level and AMR configuration, we need to exchange level and vertex.
  std::shared_ptr<rtt_mesh::Parallel_Data_Operator> par_op_;
  par_op_ = std::make_shared<rtt_mesh::Parallel_Data_Operator>(topo_);
  par_op_->fill_ghost_cell_data(comm_pat_, levels);
  par_op_->fill_ghost_cell_data(comm_pat_, vertices, D);

  par_op_->fill_ghost_cell_data(comm_pat_, neighbor_count, 2 * D);
}
//------------------------------------------------------------------------------------------------//
/**
 * \brief Construct a CAMR_Mesh object from a full set of required data.
 *
 * This constructs a ready-to use mesh object. No other initialization is required.
 *
 * \tparam I_IT - integer type
 * \tparam D_IT - double type
 *
 * \param[in] num_global_cells_ - The number of global cells in the mesh e.g., sum_p(active_size_)
 * \param[in] active_size_ - The number of active cells in the mesh
 * \param[in] full_size_   - The number of full cells in the mesh
 * \param[in] ghost_size_  - The number of ghost cells in the mesh (counting active cell only)
 * \param[in] width_  - The width of all of the cells on level 0 of the mesh.
 * \param vertex_data - An iterator into the data containing the vertex coordinates. See below for
 *                      ordering info.
 * \param level_data  - An iterator in to a integer array containing the refinement level of each
 *                      cell in the mesh
 * \param neighbor_count_in - An iterator into the data containing the number of neighbors for each
 *                      cell x face.
 * \param neighbor_values - An iterator into the data containing the cell indices for the neighbors
 *                      for each cell x face.
 * \param ghost_cells - An iterator into the data containing ghost cell IDs.
 * \param ghost_procs - An iterator into the data containing ghost processor IDs.
 * \param full2active_ - An iterator into the data containing map from full to active mesh.
 *                      full mesh includes parent elements, while active mesh only includes the top
 *                      level elements.
 * \param active2full_ - An iterator into the data containing map from active to full mesh.
 * \param[in] submesh_ - A boolean, true if the mesh is a domain decomposed part of a global mesh.
 *
 * For argument vertex_data, the data is ordered by coordinate dimension, then by cell. e.g. {x0,
 * y0, z0, x1, y1, z1, ... }
 */

//! \bug this can be fixed when 4127 is fixed via C++17 'constexpr if'
#if defined(MSVC)
#pragma warning(push)
#pragma warning(disable : 4127)
#endif

template <int D>
template <typename I_IT, typename D_IT>
CAMR_Mesh<D>::CAMR_Mesh(int64_t const num_global_cells_, size_t const active_size_,
                        size_t const full_size_, size_t const ghost_size_, double const width_,
                        D_IT vertex_data, I_IT level_data, I_IT neighbor_count_in,
                        I_IT neighbor_values, I_IT ghost_cells, I_IT ghost_procs, I_IT full2active_,
                        I_IT active2full_, bool const submesh_)
    : vertices(), levels(), maximum_level(0), number_of_levels(0), active_size(active_size_),
      full_size(full_size_), ghost_size(ghost_size_), vertices_size(active_size_ * D),
      width(width_), local_volume(0.0), widths(), submesh(submesh_) {

  num_global_cells = num_global_cells_;
  if (rtt_c4::nodes() > 1)
    decomposition = rtt_mesh::Topology::DOMAIN_DECOMPOSITION;

  // copy full2active, active2full array.
  full2active.assign(full2active_, full2active_ + full_size);
  active2full.assign(active2full_, active2full_ + active_size);

  // first find duplicates in ghost cells,
  std::map<std::pair<int, int>, int> ghost_cell_center;

  int ghost_size_new(0);
  std::map<int, int> ghost_map;

  std::vector<size_t> opposite_face{2, 1, 4, 3, 6, 5};

  // check and create list without duplicates.
  for (int i = 0; i < static_cast<int>(ghost_size); ++i) {

    if (ghost_cell_center.find(std::make_pair(ghost_cells[i], ghost_procs[i])) !=
        ghost_cell_center.end()) {
      // New entry
      // place local cell-id.

      int const val = ghost_cell_center[std::make_pair(ghost_cells[i], ghost_procs[i])];
      ghost_map[-i - 1] = -val - 1;
    } else {
      //  entry exists
      ghost_cell_center[std::make_pair(ghost_cells[i], ghost_procs[i])] = ghost_size_new;

      bnd_cells.push_back(ghost_cells[i]);
      bnd_procs.push_back(ghost_procs[i]);

      // this means -(i+1) in neighbor_values needs to be replaced with -(ghost_size_new+1),
      ghost_map[-i - 1] = -ghost_size_new - 1;
      ghost_size_new++;
    }
  } //i

  // * replace neighbor values.
  // * set the layout given the neighbor count and neighbor values
  int center(0);
  auto neighbor_begin = neighbor_count_in;
  uint32_t const num_coarse_faces = 2 * D;
  auto num_amr_faces = std::pow(2, (D - 1));

  for (int32_t cell = 1; static_cast<uint32_t>(cell) <= active_size; ++cell) {
    // Loop over the faces
    for (uint32_t face = 1; face <= num_coarse_faces; ++face) {
      // Get next value (number of neighbor faces..)
      uint32_t const number_faces = *neighbor_count_in++;

      if (number_faces <= num_amr_faces) { // face is refined.
        for (uint32_t amr_face = 1; amr_face <= number_faces; ++amr_face) {
          if (neighbor_values[center] < 0) {
            int const id = neighbor_values[center];
            neighbor_values[center] = ghost_map[id];
          }
          center++;
        }
      } else {
        std::ostringstream out;
        out << "Invalid number of neighbors(0): " << number_faces << " specified in CAMR_Layout";
        Insist(0, out.str());
      }
    }
  }
  neighbor_count_in = neighbor_begin;

  // do some communication, this fct fill levsls, vertex and neighbor_count for parallel.
  constructor_body(level_data, vertex_data, neighbor_count_in);

  layout =
      std::make_shared<rtt_mesh::CAMR_Layout>(active_size_, ghost_size_, D, neighbor_count.data());

  // set the layout given the neighbor count and neighbor values
  for (int32_t cell = 1; static_cast<uint32_t>(cell) <= active_size; ++cell) {
    size_t const my_level = get_level(cell);
    U3 my_centroid = get_cell_centroid(cell);

    // Loop over the faces
    for (uint32_t face = 1; face <= layout->num_coarse_faces(); ++face) {
      // Get next value (number of neighbor faces..)
      uint32_t const number_faces =
          neighbor_count[(cell - 1) * layout->num_coarse_faces() + (face - 1)];

      if (number_faces == 1) { //no amr.
        /////////////////////////
        // neighbor is ghost cell.
        if (*neighbor_values < 0) {

          // neighbor id. 1-base (appended to active cells).
          int const cid = static_cast<int>(active_size) - *neighbor_values;

          // we need to identify, this is cell is  bottom/top (lower/upper)
          size_t const neighbor_level = get_level(cid);

          // check the levels between my cell and neighbor.
          if (neighbor_level == my_level) // the same level.
          {
            for (uint32_t f = 1; f <= layout->num_amr_faces(); ++f) {
              ((*layout))(cell, face, f) = *neighbor_values;
              (*layout)(cid, opposite_face[face - 1], f) = cell;
            }

          } else if (neighbor_level + 1 == my_level) // neighbor is coarser
          {
            U3 neighbor_centroid = get_cell_centroid(cid);

            if (D == 2) {
              if (face <= 2) { // edge along y-axis
                if (my_centroid[1] < neighbor_centroid[1]) {
                  (*layout)(cid, opposite_face[face - 1], 1) = cell;

                } else {
                  (*layout)(cid, opposite_face[face - 1], 2) = cell;
                }
              } else { //edge along x-axis
                if (my_centroid[0] < neighbor_centroid[0]) {
                  (*layout)(cid, opposite_face[face - 1], 1) = cell;
                } else {
                  (*layout)(cid, opposite_face[face - 1], 2) = cell;
                }
              }

            } else //3D
            {
              // identify axis.
              int axis0 = 1;
              int axis1 = 2;

              if (face <= 2) {        //y-z
              } else if (face <= 4) { //x-z
                axis0 = 0;
              } else { //x-y
                axis0 = 0;
                axis1 = 1;
              }

              if (my_centroid[axis0] < neighbor_centroid[axis0]) {
                if (my_centroid[axis1] < neighbor_centroid[axis1])
                  (*layout)(cid, opposite_face[face - 1], 1) = cell;
                else
                  (*layout)(cid, opposite_face[face - 1], 3) = cell;
              } else {
                if (my_centroid[axis1] < neighbor_centroid[axis1])
                  (*layout)(cid, opposite_face[face - 1], 2) = cell;
                else
                  (*layout)(cid, opposite_face[face - 1], 4) = cell;
              }
            }
          } else {
            std::cout << "should not be here.. " << neighbor_level << " " << my_level << '\n';
            throw;
          }
        }
        // * this is terrible, but layout has a size of nfine_face*cells.
        // * and if the face is not refined, we need to fill the all of them with the same neighbor
        //   values..
        for (uint32_t f = 1; f <= layout->num_amr_faces(); ++f)
          (*layout)(cell, face, f) = *neighbor_values;
        neighbor_values++;
      } else if (number_faces ==
                 layout->num_amr_faces()) { // face is refined. (neighbor is refined)
        for (uint32_t amr_face = 1; amr_face <= layout->num_amr_faces(); ++amr_face) {
          //        for (uint32_t amr_face = 1; amr_face <= number_faces; ++amr_face){
          if (*neighbor_values < 0) {
            // neighbor is refined, so cell is the only neighbor of cid.
            int const cid = static_cast<int>(active_size) - *neighbor_values;
            //            layout(cid,opposite_face[face-1],amr_face) = cell;
            (*layout)(cid, opposite_face[face - 1], 1) = cell;
          }

          (*layout)(cell, face, amr_face) = *neighbor_values++;
        }

      } else {
        std::ostringstream out;
        out << "Invalid number of neighbors(1): " << number_faces << " specified in CAMR_Layout";
        Insist(0, out.str());
      }
    }
  }

  Ensure(layout->validate() == 0);
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Construct a CAMR_Mesh object from a full set of required data.
 *
 * This constructs a ready-to use mesh object. No other initialization is required.
 *
 * \tparam I_IT - integer type
 * \tparam D_IT - double type
 *
 * \param[in] size_   - The number of cells in the mesh
 * \param[in] width_  - The width of all of the cells on level 0 of the mesh.
 * \param vertex_data - An iterator into the data containing the vertex coordinates. See below for
 *                      ordering info.
 * \param level_data  - An iterator in to a integer array containing the refinement level of each
 *                      cell in the mesh
 * \param neighbor_count_in - An iterator into the data containing the number of neighbors for each
 *                      cell x face.
 * \param neighbor_values - An iterator into the data containing the cell indices for the neighbors
 *                      for each cell x face.
 * \param[in] submesh_ - A boolean, true if the mesh is a domain decomposed part of a global mesh.
 *
 * For argument vertex_data, the data is ordered by coordinate dimension, then by cell. e.g. {x0,
 * y0, z0, x1, y1, z1, ... }
 */
template <int D>
template <typename I_IT, typename D_IT>
CAMR_Mesh<D>::CAMR_Mesh(size_t const size_, double const width_, D_IT vertex_data, I_IT level_data,
                        I_IT neighbor_count_in, I_IT neighbor_values, bool const submesh_)
    //    : layout(size_, 2 * D, static_cast<uint32_t>(std::pow(2, D - 1))),
    : vertices(), levels(), maximum_level(0), number_of_levels(0), active_size(size_),
      full_size(size_), ghost_size(0), vertices_size(size_ * D), width(width_), local_volume(0.0),
      widths(), submesh(submesh_) {
  num_global_cells = size_;

  layout = std::make_shared<rtt_mesh::CAMR_Layout>(active_size, ghost_size, D, neighbor_count_in);

  // set the layout given the neighbor count and neighbor values
  for (int32_t cell = 1; static_cast<uint32_t>(cell) <= active_size; ++cell) {
    // Loop over the faces
    for (uint32_t face = 1; face <= layout->num_coarse_faces(); ++face) {
      // Get next value (number of neighbor faces..)
      uint32_t const number_faces = *neighbor_count_in++;

      if (number_faces == 1) { //no amr.
        for (uint32_t amr_face = 1; amr_face <= layout->num_amr_faces(); ++amr_face)
          (*layout)(cell, face, amr_face) = *neighbor_values;
        neighbor_values++;
      } else if (number_faces == layout->num_amr_faces()) { // face is refined.
        for (uint32_t amr_face = 1; amr_face <= layout->num_amr_faces(); ++amr_face)
          (*layout)(cell, face, amr_face) = *neighbor_values++;
      } else {
        std::ostringstream out;
        out << "Invalid number of neighbors: " << number_faces << " specified in CAMR_Layout";
        Insist(0, out.str());
      }
    }
  }
  Ensure(layout->validate() == 0);

  constructor_body(level_data, vertex_data, neighbor_count_in);
}
#if 0
//------------------------------------------------------------------------------------------------//
/**
 * \brief Construct a CAMR_Mesh object from geometric data and a layout object.
 *
 * This constructs a ready-to use mesh object. No other initialization is required.
 *
 * \param[in] width_      - The width of all of the cells on level 0 of the mesh.
 * \param     vertex_data - An iterator into the data containing the vertex coordinates. See below
 *                          for ordering info.
 * \param     level_data  - An iterator in to a integer array containing the refinement level of
 *                          each cell in the mesh
 * \param     layout_     - An object of type rtt_tracking::AMR_Layout.
 * \param[in] submesh_    - A boolean, true if the mesh is a domain decomposed part of a global
 *                          mesh
 *
 * For argument vertex_data, the data is ordered by coordinate dimension, then by cell. e.g. {x0,
 * y0, z0, x1, y1, z1, ... }
 */
template <int D>
template <typename I_IT, typename D_IT>
CAMR_Mesh<D>::CAMR_Mesh(double const width_, D_IT vertex_data, I_IT level_data,
                        rtt_mesh::CAMR_Layout const &layout_, bool const submesh_)
    : layout(layout_), vertices(), levels(), maximum_level(0), number_of_levels(0),
      active_size(layout_.get_num_owned_cells()), // number of cells in the mesh
      full_size(layout_.get_num_owned_cells()), ghost_size(0),
      vertices_size(layout_.get_num_owned_cells() * D),
      width(width_), local_volume(0.0), widths(), submesh(submesh_) {
  Ensure(layout.validate() == 0);
  constructor_body(level_data, vertex_data);
}
#endif
//------------------------------------------------------------------------------------------------//
/**
 * \brief Return true if the given cell face is refined (i.e. has more than one cell on the other
 *        side)
 * \param[in] cell The cell
 * \param[in] coarse_face The coarse face of the cell.
 */
template <int D> bool CAMR_Mesh<D>::is_refined(int const cell, uint32_t const coarse_face) const {
  return (layout->num_cells_across(cell, coarse_face) > 1);
}

//------------------------------------------------------------------------------------------------//
/*!
  * \brief Determines if the specified cell is coarse relative to its neighbor across a given face.
  *
  * This function checks whether the neighboring cell across the specified face is refined or
  * coarse. It handles cases for boundary conditions, off-processor neighbors, and different
  * dimensions. The function returns:
  * - 0 if the cell across the face is refined or if it is a boundary cell.
  * - 1 if the neighbor cell is coarser and this is a 2D mesh with the same cell.
  * - 2 if the neighbor cell is coarser and this is a 2D mesh with a different cell.
  * - A positive integer indicating the index of the cell in the next cell list if the neighbor is
  *   coarser.
  *
  * \tparam D The dimension of the mesh (e.g., 2 for 2D meshes, 3 for 3D meshes).
  *
  * \param[in] cell The identifier of the current cell being examined.
  * \param[in] face The identifier of the face across which the neighbor cell is being examined.

  * \return An integer indicating the relationship between the current cell and its neighbor:
  *         - 0 if refined or boundary,
  *         - 1 or 2 if coarser in 2D,
  *         - or a positive index if coarser in higher dimensions.
  */
template <int D> int CAMR_Mesh<D>::is_coarse(int const cell, uint32_t const face) const {

  std::vector<int> opposite_face{2, 1, 4, 3, 6, 5};
  auto neighbor_cells = next_cells(cell, face);

  // Determine if the neighbor cell is valid and refined
  if (neighbor_cells.size() > 1) {
    return 0;
  }

  int cid = neighbor_cells[0];
  // Handle boundary or reflection case
  if (cid == 0 || cid == cell) {
    return 0;
  }

  // Handle off-proc neighbor case
  if (cid < 0) {
    cid = static_cast<int>(active_size) - cid;
  }

  // List of cells across face (from cid, neighbor cell)
  auto next_cell = next_cells(cid, opposite_face[face - 1]);

  // Return early if next_cell size is 1
  if (next_cell.size() == 1) {
    return 0;
  }

  // Determine if neighbor is coarser
  if (D == 2) {
    return next_cell[0] == cell ? 1 : 2;
  }

  // Search for the cell in next_cell
  for (size_t j = 0; j < next_cell.size(); ++j) {
    if (next_cell[j] == cell) {
      return static_cast<int>(j + 1);
    }
  }

  return 0; // Default case
}

#if defined(MSVC)
#pragma warning(pop)
#endif

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return a vector of the cells across a given cell face
 * \param[in] cell The cell
 * \param[in] coarse_face The coarse face
 */
template <int D>
typename CAMR_Mesh<D>::sf_int CAMR_Mesh<D>::next_cells(const int cell,
                                                       const int coarse_face) const {
  return (*layout)(cell, coarse_face);
}
template <int D>
typename CAMR_Mesh<D>::sf_int CAMR_Mesh<D>::next_cells(const size_t cell,
                                                       const size_t coarse_face) const {
  Check(cell < INT32_MAX && coarse_face < INT32_MAX);
  return (*layout)(static_cast<int>(cell), static_cast<int>(coarse_face));
}

//------------------------------------------------------------------------------------------------//
//! Return a pointer to the coordinate data for a given cell.
template <int D>
typename CAMR_Mesh<D>::sf_double::const_iterator CAMR_Mesh<D>::vertex_data(size_t cell) const {
  Require(cell >= 1 && cell <= layout->get_num_total_cells());
  return vertices.begin() + (cell - 1) * D; //XXX-1
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return the requested dimension of the cell
 * \arg coordinate The dimension. One-based.
 * \arg cell The cell
 *
 * This function dispatches to get_width(cell). We do not need the coordinate information since the
 * cells are all uniformly sized in all D dimensions.
 */
template <int D>
double CAMR_Mesh<D>::cell_width(size_t const Remember(coord), size_t const cell) const {
  Check(coord >= 1);
  Check(coord <= D);
  return get_cell_width(cell);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Calculate the face index corresponding to the cell across the face of a given cell. This
 *        function cannot be called for vacuum or reflective boundaries
 *
 * \param[in] cell        - cell index
 * \param[in] coarse_face - coarse_face index
 * \return index of corresponding coarse_face
 */
template <int D>
size_t CAMR_Mesh<D>::next_face(int const Remember(cell), uint32_t const coarse_face) const {
  using std::count;

  Remember(sf_int neighbors((*layout)(cell, coarse_face)));
  Require(count(neighbors.begin(), neighbors.end(), 0) == 0);
  Require(count(neighbors.begin(), neighbors.end(), cell) == 0);

  switch (coarse_face) {
  case 1:
    return 2; // low-x face  -> high-x face
  case 2:
    return 1; // high-x face -> low-x face
  case 3:
    return 4; // low-y face  -> high-y face
  case 4:
    return 3; // high-y face -> low-y face
  case 5:
    return 6; // low-z face  -> high-z face
  case 6:
    return 5; // high-z face -> low-z face
  default:
    Insist(0, "Invalid Face Index!"); // LCOV_EXCL_LINE
    return 0;                         // LCOV_EXCL_LINE
  }
}
template <int D> size_t CAMR_Mesh<D>::next_face(size_t const cell, size_t const coarse_face) const {
  Check(cell < INT32_MAX && coarse_face < UINT32_MAX);
  return CAMR_Mesh<D>::next_face(static_cast<int>(cell), static_cast<uint32_t>(coarse_face));
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Find face index given position and cell index
 *
 * \param[in] position - position in the cell
 * \param[in] cell     - cell index
 * \return face index, else 0 if r is not on a face
 */
template <int D> uint32_t CAMR_Mesh<D>::get_face(U3 const &position, size_t const cell) const {
  using rtt_dsxx::soft_equiv;

  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_owned_cells());
  Check(cell < INT32_MAX);
  Require(in_cell(static_cast<int>(cell), position));

  // get cell vertex and size
  sf_double const cell_vertex = get_vertex(cell);
  double const delta = get_cell_width(cell);

  // high and low sides of cell (for a given dimension)
  auto high = static_cast<double>(FLT_MIN);
  auto low = static_cast<double>(FLT_MIN);

  // loop over dimensions and find face
  for (int d = 0; d < D; ++d) {
    // assign high and low sides for this dimension
    high = cell_vertex[d] + delta;
    low = cell_vertex[d];

    // check if position is on low face
    if (soft_equiv(position[d], low))
      return 2 * d + 1;
    // check if position is on high face
    else if (soft_equiv(position[d], high))
      return 2 * d + 2;
  }
  // face was not found
  return 0;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return a vector<vector> containing the vertices of the cell.
 *
 * The data is in vertex, dimension order, i.e. inner vector is over the vertices, the outer is over
 * the dimensions.
 *
 * \param[in] cell The cell
 * \return A vector<vector> containing the vertices of the cell.
 */
template <int D>
typename CAMR_Mesh<D>::sf_double CAMR_Mesh<D>::get_short_vertices(size_t cell) const {

  // for distance to boundary we only need the low dimension values and the delta
  Require(cell > 0 && static_cast<uint32_t>(cell) <= layout->get_num_total_cells());

  sf_double myVertices = get_vertex(cell);
  myVertices.push_back(get_cell_width(cell));
  return myVertices;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief calculate the center of a cell face
 *
 * This function calculates the geometric center of a cell face.
 *
 * \param[in] cell - cell index
 * \param[in] face - face index
 * \return face center
 */
template <int D> U3 CAMR_Mesh<D>::get_face_center(const size_t cell, const uint32_t face) const {
  Require(cell >= 1 && static_cast<uint32_t>(cell) <= layout->get_num_owned_cells());
  Require(face >= 1 && face <= layout->num_coarse_faces());

  // get the lowest cell vertex
  sf_double const cell_vertex = get_vertex(cell);

  // get the cell width
  double const delta = get_cell_width(cell);

  // get the face data
  std::pair<int, bool> const face_data = layout->face_axis_sign(face);

  // get the dimension (x, y, or z) of the face
  int const &face_dimension = face_data.first;

  // get whether the face is high (true) or low (false)
  bool const &face_positive = face_data.second;

  // calculate the face center
  U3 center{{0.0, 0.0, 0.0}};
  for (int d = 1; d <= D; ++d) {
    if (d == face_dimension) {
      // this is the same dimension as the face, set the coordinate to face coordinate
      center[d - 1] = cell_vertex[d - 1] + (face_positive ? delta : 0);
    } else {
      // for other dimensions, the coordinates correspond to the center of the cell
      center[d - 1] = cell_vertex[d - 1] + 0.5 * delta;
    }
  }

  Ensure(get_face(center, cell) == face);
  return center;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief calculate the center of a cell
 *
 * This function calculates the geometric center of a cell.
 *
 * \param[in] cell index to cell to be queried
 * \return cell center
 */
template <int D> U3 CAMR_Mesh<D>::get_cell_centroid(const size_t cell) const {
  Require(cell >= 1 && static_cast<uint32_t>(cell) <= layout->get_num_total_cells());

  // get the lowest cell vertex
  sf_double const cell_vertex = get_vertex(cell);

  // get the cell width
  double const delta = get_cell_width(cell);

  U3 center{{0.0, 0.0, 0.0}};
  for (int d = 1; d <= D; ++d) {
    center[d - 1] = cell_vertex[d - 1] + 0.5 * delta;
  }
  return center;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Returns a vector of valid face indices
 *
 * \return vector of valid face indices
 */
template <int D> typename CAMR_Mesh<D>::sf_int CAMR_Mesh<D>::get_valid_cell_faces() const {
  constexpr std::array<int, 6> f = {1, 2, 3, 4, 5, 6};
  return {f.begin(), f.begin() + D * 2};
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Calculate an intersection location on an CAMR_Mesh face
 *
 * \param[in] r position
 * \param[in] omega direction
 * \param[in] distance minimum distance to boundary
 * \param[in] cell index to the cell containing position r
 * \param[in] face the pre-calculated face containing the intersection
 *
 * \return intersection point on face
 */
template <int D>
U3 CAMR_Mesh<D>::get_intersection(const U3 &r, const U3 &omega, const double &distance,
                                  const int &cell, const uint32_t face) const {
  sf_double const cell_vertex = get_vertex(cell);
  U3 intersection{0, 0, 0};
  switch (face) {
  case 1:
    intersection[0] = cell_vertex[0];
    intersection[1] = r[1] + distance * omega[1];
    intersection[2] = r[2] + distance * omega[2];
    break;
  case 2:
    intersection[0] = cell_vertex[0] + get_cell_width(cell);
    intersection[1] = r[1] + distance * omega[1];
    intersection[2] = r[2] + distance * omega[2];
    break;
  case 3:
    intersection[0] = r[0] + distance * omega[0];
    intersection[1] = cell_vertex[1];
    intersection[2] = r[2] + distance * omega[2];
    break;
  case 4:
    intersection[0] = r[0] + distance * omega[0];
    intersection[1] = cell_vertex[1] + get_cell_width(cell);
    intersection[2] = r[2] + distance * omega[2];
    break;
  case 5:
    intersection[0] = r[0] + distance * omega[0];
    intersection[1] = r[1] + distance * omega[1];
    intersection[2] = cell_vertex[2];
    break;
  case 6:
    intersection[0] = r[0] + distance * omega[0];
    intersection[1] = r[1] + distance * omega[1];
    intersection[2] = cell_vertex[2] + get_cell_width(cell);
    break;
  default:
    Insist(0, "Invalid Face Index!"); // LCOV_EXCL_LINE
    break;
  }

  return intersection;
}

} // namespace rtt_mesh

#endif // mesh_CAMR_Mesh_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/CAMR_Mesh.hh
//------------------------------------------------------------------------------------------------//
