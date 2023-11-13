//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/CAMR_Layout.cc
 * \author H. Park
 * \date   May 10 2023
 * \brief  AMR_Layout member definitions.
 * \note   Copyright (C) 2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "CAMR_Layout.hh"
#include "c4/C4_Functions.hh"
#include "ds++/Endian.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
namespace rtt_mesh {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//------------------------------------------------------------------------------------------------//

/*!
 * \brief Constructor when the problem dimension is known up front.
 *
 * This constructor allows the user to set up the data inside of CAMR_Layout based on the number of
 * cells in the layout. The number of faces per cell is \c (ndim*2)*2^(ndim-1) to work with AMR
 * meshes. Current AMR hydro codes at LANL do not allow more than one level of refinement across a
 * cell face. If we need more flexibility we should make an Unstructured_Layout class. The one thing
 * standing in the way of this is that RZWedge is 3D and only has 2 AMR cells across each face. When
 * this class goes away we can simplify the constructor to use only num_cells and n_dim.
 *
 * \param[in] num_cells number of cells in the layout
 * \param[in] faces coarse_faces (generally 2*ndim)
 * \param[in] amr_faces generally 2^(ndim-1)
 * \param[in] nghosts number of ghost cells, zero as default.
 */
CAMR_Layout::CAMR_Layout(size_t num_cells, uint32_t faces, uint32_t amr_faces, size_t nghosts)
    : n_coarse_faces(faces), n_fine_faces(amr_faces),
      n_faces_per_cell(n_coarse_faces * n_fine_faces), num_owned_cells(num_cells),
      num_ghost_cells(nghosts), cell_face(sf_int((num_cells + nghosts) * n_faces_per_cell)) {
  num_cells_across_face.resize(n_coarse_faces * (num_cells + nghosts));
  Ensure(faces > 0);
  Ensure(amr_faces > 0);
  Ensure(cell_face.size() == (num_cells + nghosts) * n_faces_per_cell);
}

CAMR_Layout::CAMR_Layout(size_t num_active_cells, size_t nghosts, uint32_t ndim,
                         int *neighbor_counts)
    : n_coarse_faces(2 * ndim), n_fine_faces(static_cast<uint32_t>(std::pow(2, ndim - 1))),
      n_faces_per_cell(n_coarse_faces * n_fine_faces), num_owned_cells(num_active_cells),
      num_ghost_cells(nghosts), cell_face(sf_int((num_active_cells + nghosts) * n_faces_per_cell)) {

  size_t tot_ncaf = n_coarse_faces * (num_owned_cells + nghosts);

  num_cells_across_face.assign(neighbor_counts, neighbor_counts + tot_ncaf);

  Ensure(ndim > 0);
  Ensure(cell_face.size() == (num_owned_cells + nghosts) * n_faces_per_cell);
  Ensure(num_cells_across_face.size() == tot_ncaf);
}

//------------------------------------------------------------------------------------------------//
// PUBLIC MEMBER FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Print out the AMR_Layout.
 *
 * \param[in,out] out \c std::ostream that will capture the report
 * \param cell print out layout for given cell index
 *
 * Return a updated \c std::ostream, appended with the AMR_Layout report.
 */
void CAMR_Layout::print(std::ostream &out, uint32_t cell) const {
  using std::endl;
  using std::setw;

  Require(cell > 0 && cell <= get_num_total_cells());

  // Preserve the "classic" non-AMR printing
  {
    // print out the cell info
    out << "====================" << endl;
    out << "Cell         : " << cell << endl;
    out << "Coarse Faces : " << n_coarse_faces << endl;
    out << "--------------------" << endl;
    for (uint32_t i = 0; i < n_coarse_faces; ++i) {
      out << setw(4) << i + 1 << endl;
      sf_int coarse_nbrs = operator()(cell, i + 1);
      int j = 0;
      for (auto nbr : coarse_nbrs) {
        out << setw(8) << j + 1 << setw(5) << nbr << endl;
        j++;
      }
    }
    out << "====================" << endl;
  }
}

//------------------------------------------------------------------------------------------------//
//! Test for equality.
bool CAMR_Layout::operator==(CAMR_Layout const &rhs) const {
  // if the cell data is equal then return true
  return cell_face == rhs.cell_face;
}

//------------------------------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//------------------------------------------------------------------------------------------------//

//! Print out all cells using operator<<.
std::ostream &operator<<(std::ostream &out, const CAMR_Layout &object) {
  auto const num_cells(object.get_num_total_cells());
  for (uint32_t i = 1; i <= num_cells; ++i)
    object.print(out, i);
  return out;
}

//------------------------------------------------------------------------------------------------//
//! Make sure that a meshes neighbor relationship is valid (graph is cyclic)
uint32_t CAMR_Layout::validate() const {

  std::array<uint32_t, 7> constexpr opposites = {0, 2, 1, 4, 3, 6, 5};

  for (uint32_t cell = 1; cell <= get_num_owned_cells(); ++cell) {
    for (uint32_t face = 1; face <= n_coarse_faces; ++face) {
      uint32_t const opposite_face = opposites[face];
      sf_int neighbors = (*this)(cell, face); // negative value indicates off-processor index.
      for (auto n : neighbors) {
        // don't check for edge cells or boundary cells
        if (n != static_cast<int>(cell) && n >= 1) {
          Check(static_cast<uint32_t>(n) < UINT32_MAX);
          sf_int opp_neighbors = (*this)(static_cast<uint32_t>(n), opposite_face);
          // count the number of times that the current cell shows up in the cell across a face's
          // neighbor list. It should be more than 0.
          if (std::count(opp_neighbors.begin(), opp_neighbors.end(), cell) == 0) {
            return cell;
          }
        }
      }
    }
  }
  return 0;
}

} // end namespace rtt_mesh

//------------------------------------------------------------------------------------------------//
// end of CAMR_Layout.cc
//------------------------------------------------------------------------------------------------//
