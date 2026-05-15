//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/CAMR_Layout.hh
 * \author H. Park
 * \brief  CAMR_Layout class definition.
 * \note   Copyright (C) 2023-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_mesh_CAMR_Layout_hh
#define rtt_mesh_CAMR_Layout_hh

#include "dsxx/Assert.hh"
#include <climits> // UINT32_MAX
#include <cstdint> // uint32_t
#include <iosfwd>
#include <map>
#include <memory>
#include <utility>
#include <vector>

// cSpell:ignore nghosts, ndim

namespace rtt_mesh {

//================================================================================================//
/*!
 * \class CAMR_Layout

 * \brief Provides cell connectivity information for 1-N meshes.

 * Class that defines the cell-face-cell connectivity for 1-to-N cell meshes.  1-to-N describes
 * meshes whose cells may have up to N cell neighbors across a face.  CAMR_Layout, like Layout, is
 * independent of the geometry of the mesh.

 * This class is designed for use in particular types of meshes; particularly those that have 1-to-N
 * cells across faces.  Therefore, we have not designed this into a inheritance hierarchy with
 * Layout.  Each mesh will use the appropriate services that are required by a particular Layout.
 * The actual interface to the outside world is through the mesh type, not the layout class.  The
 * only classes that see the layout are the mesh and its builder; thus, there is no need to restrain
 * the interface to a common set of functions declared by a base class.

 * The ordering of cells across a face is row-wise starting from the lowest coordinate: thus across
 * the hi x face in XY, 1=lo y and 2=hi y.  Similarly, across the hi y face in XYZ 1=lo x lo z, 2=hi
 * x lo z, 3=lo x hi z, and 4=hi x hi z.  There is no formal requirement to follow this methodology;
 * however, for consistency one should attempt to use this convention.

 * The CAMR_Layout uses the convention of coarse faces.  For example, in a 2-D RZ mesh each cell has
 * 4 coarse faces.  However, if one of the sides of a cell borders a refined cell then there are 2
 * fine faces (2 cell neighbors) on that coarse face.  This presents a "hierarchical" notion of a
 * refined mesh.
 */
//================================================================================================//

class CAMR_Layout {

public:
  // Important typedefs.
  using sf_int = std::vector<int>;
  using SP_Layout = std::shared_ptr<CAMR_Layout>;

private:
  // min number of faces per cell.
  uint32_t n_coarse_faces;
  // number of refined faces per coarse face.
  uint32_t n_fine_faces;
  // max number of faces per cell.
  uint32_t n_faces_per_cell;

  size_t num_owned_cells;
  size_t num_ghost_cells;

  //! Cell-face data, indexing is 0:N*n_faces_per_cell-1.
  sf_int cell_face;

  // how many cells across face?
  sf_int num_cells_across_face;

public:
  /*! \brief Constructor
   * \param[in] num_active_cells number of active cells
   * \param[in] nghosts number of ghost cells, zero as default.
   * \param[in] ndim dimensionality?
   * \param[in] neighbor_count pointer to neighbor cell mapping data? */
  CAMR_Layout(size_t num_active_cells, size_t nghosts, uint32_t ndim, int *neighbor_count);

  // >>> ACCESSORS

  //! Get the number of cells in the layout.
  //  size_t get_num_cells() const { return (cell_face.size() / n_faces_per_cell); }
  [[nodiscard]] size_t get_num_owned_cells() const { return num_owned_cells; }
  [[nodiscard]] size_t get_num_ghost_cells() const { return num_ghost_cells; }
  [[nodiscard]] size_t get_num_total_cells() const { return num_owned_cells + num_ghost_cells; }

  //! Get the number of coarse faces for this layout.
  [[nodiscard]] inline uint32_t num_coarse_faces() const { return n_coarse_faces; }

  //! Get the number of AMR faces on a coarse faces for this layout
  [[nodiscard]] inline uint32_t num_amr_faces() const { return n_fine_faces; }

  [[nodiscard]] inline uint32_t num_max_faces_per_cell() const { return n_faces_per_cell; }

  //! Number of cells across a coarse face (Fortran Indexing)
  [[nodiscard]] inline uint32_t num_cells_across(uint32_t const cell,
                                                 uint32_t const coarse_face) const;
  [[nodiscard]] inline uint32_t num_cells_across(size_t const cell, size_t const coarse_face) const;

  //! Get the indices of cells on the boundary at a given direction
  [[nodiscard]] inline sf_int get_edge_cells(const uint32_t face) const {
    Require(face >= 1 && face <= n_coarse_faces);
    sf_int boundary_cells;
    for (uint32_t cell = 1; cell <= get_num_owned_cells(); ++cell) {
      int const neighbor =
          cell_face[(static_cast<size_t>(cell) - 1) * n_faces_per_cell + n_fine_faces * (face - 1)];
      if (neighbor == 0 || std::cmp_equal(neighbor, cell)) {
        boundary_cells.push_back(cell);
      }
    }
    return boundary_cells;
  }

  //! Return the axis and direction of a face
  [[nodiscard]] inline std::pair<uint32_t, bool> face_axis_sign(const uint32_t face) const {
    Require(face >= 1 && face <= n_coarse_faces);
    return std::make_pair((face + 1) / 2, static_cast<bool>((face + 1) % 2));
  }

  //! Make sure that a meshes neighbor relationship is valid (graph is cyclic)
  [[nodiscard]] uint32_t validate() const;

  [[nodiscard]] const std::vector<int> &get_internal_data() const { return cell_face; }

  // >>> OVERLOADED OPERATORS

  // Overloaded subscripting operators for assignment and retrieval.
  inline int operator()(size_t cell, const size_t coarse_face, const size_t fine_face) const;
  inline int &operator()(size_t cell, const size_t coarse_face, const size_t fine_face);
  inline sf_int operator()(uint32_t cell, const uint32_t coarse_face) const;

  // Overloaded equality operators.
  bool operator==(const CAMR_Layout &rhs) const;
  bool operator!=(const CAMR_Layout &rhs) const { return !(*this == rhs); }

  // Diagnostic functions.
  void print(std::ostream &out, uint32_t cell) const;
};

//------------------------------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Overload operator for stream output.
 *
 * \param[in,out] out    std::ostream object that receives text streams.
 * \param[in]     object Report data from this AMR_Layout object.
 */
std::ostream &operator<<(std::ostream &out, const CAMR_Layout &object);

//------------------------------------------------------------------------------------------------//
// INLINE FUNCTIONS FOR AMR_LAYOUT
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
//! Return the number of cells across a coarse face.
uint32_t CAMR_Layout::num_cells_across(uint32_t const cell, uint32_t const coarse_face) const {
  Require(cell > 0 && cell <= get_num_total_cells());
  Require(coarse_face > 0 && coarse_face <= n_coarse_faces);
  Check((cell - 1) * n_faces_per_cell < INT32_MAX && (coarse_face - 1) * n_fine_faces < INT32_MAX);

  uint32_t const start = (cell - 1) * n_coarse_faces + (coarse_face - 1);

  Check(start < num_cells_across_face.size());

  return num_cells_across_face[start];
}
uint32_t CAMR_Layout::num_cells_across(size_t const cell, size_t const coarse_face) const {
  Check(cell < UINT32_MAX && coarse_face < UINT32_MAX);
  return CAMR_Layout::num_cells_across(static_cast<uint32_t>(cell),
                                       static_cast<uint32_t>(coarse_face));
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Access the cell across a coarse/fine face combination.
 * \param cell cell index
 * \param coarse_face coarse face index
 * \param fine_face fine face index
 * \return cell index on other side of cell, coarse_face, fine_face triplet
 */
int CAMR_Layout::operator()(size_t const cell, const size_t coarse_face,
                            const size_t fine_face) const {
  Require(cell <= get_num_total_cells());
  Require(coarse_face <= n_coarse_faces);
  Require(fine_face <= n_fine_faces);
  return cell_face[(cell - 1) * n_faces_per_cell + (coarse_face - 1) * n_fine_faces +
                   (fine_face - 1)];
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Mutable access for the cell across a coarse/fine face combination.

 * This overloaded operator allows the user to change the value of the cell across a face.  This
 * operator should be used to build the CAMR_Layout.

 * \param cell cell index
 * \param coarse_face coarse face index
 * \param fine_face fine face index

 * \return mutable cell index on other side of cell, coarse_face, fine_face triplet

 */
int &CAMR_Layout::operator()(const size_t cell, const size_t coarse_face, const size_t fine_face) {
  Require(cell <= get_num_total_cells());
  Require(coarse_face <= n_coarse_faces);
  Require(fine_face <= n_fine_faces);
  return cell_face[(cell - 1) * n_faces_per_cell + (coarse_face - 1) * n_fine_faces +
                   (fine_face - 1)];
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return a vector of cells across a coarse face.
 *
 * This is a non-mutable look-up option for determining the cells across a coarse face.
 *
 * \param cell cell index
 * \param coarse_face coarse face index
 * \return vector of cells across the coarse face
 */
CAMR_Layout::sf_int CAMR_Layout::operator()(uint32_t cell, const uint32_t coarse_face) const {

  Require(cell <= get_num_total_cells());
  Require(coarse_face <= n_coarse_faces);

  int const start = (cell - 1) * n_faces_per_cell + (coarse_face - 1) * n_fine_faces;
  int const end = start + num_cells_across(cell, coarse_face);

  return {cell_face.begin() + start, cell_face.begin() + end};
}

} // namespace rtt_mesh

#endif // rtt_mesh_CAMR_Layout_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/CAMR_Layout.hh
//------------------------------------------------------------------------------------------------//
