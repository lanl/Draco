//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/test/tstCAMR_Mesh.cc
 * \author H. Park
 * \date   Mon Jul 10 2023
 * \brief  Unit test for CAMR_Mesh
 * \note   Copyright (C) 2023-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "c4/ParallelUnitTest.hh"
#include "mesh/CAMR_Layout.hh"
#include "mesh/CAMR_Mesh.hh"
#include "viz/Ensight_Translator.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

using namespace std;
using namespace rtt_mesh;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
struct CAMR_Mesh_Data {
  //! RZ data and AMR connectivity
  int *num_across_face = nullptr;  //!< Num. cells across a coarse face. (RZ, XYZ)
  int *cell_across_face = nullptr; //!< Cells across a face. (All)
  double *rx_extent = nullptr;     //!< low/hi X/R dimensions per cell, cm (R, RZ).
  double *z_extent = nullptr;      //!< low/hi Z dimensions per cell, cm (RZ only).

  //! XY/XYZ only data:
  double *xyz_extents = nullptr; //!< Lowest vertex of each cell (XY/XYZ only).
  int *cell_levels = nullptr;    //!< Refinement levels of active cells (XY/XYZ).
  double coarse_cell_size = 0.0; //!< Width of cells in coarse mesh,cm (XY/XYZ).

  //! Mapping from full mesh to active mesh. (index is assumed to be local to each proc)
  int num_full_cells = 0;             //!< Number of local cells on full mesh. (All)
  int num_active_cells = 0;           //!< Number of local active cells (not ghost cells)
  int *full_to_active_cell = nullptr; //!< Mapping from full-active mesh (GHOST_CELLS)
  int *active_to_full_cell = nullptr; //!< Mapping from active-full mesh (no GHOST_CELLS)

  //! Domain decomposed topology information
  int num_dims = 0;

  //! Number of total active cells (DOMAIN_DECOMPOSITION)
  int64_t number_global_cells = 0;
  //! Number of boundary cell on processor (DOMAIN_DECOMPOSITION)
  size_t num_bnd_cells = 0;
  //! Local, "full cell" indices (DOMAIN_DECOMPOSITION)
  int *processor_boundary_cells = nullptr;
  //! Processors where boundary cells live (DOMAIN_DECOMPOSITION)
  int *processor_boundary_procs = nullptr;

  //! Default constructor sets all pointers to nullptr.
  CAMR_Mesh_Data() = default;
  //! Checks for validity of pointers given the problem type
};

void test_2D_mesh(rtt_c4::ParallelUnitTest &ut) {
  /*
   * Create a simple 2D CAMR_Mesh
   *
   * ~~~
   *
   * // --------------------
   * // |       |       || -4
   * // |       |       ||---
   * // | 2     | 3     || -3
   * // ---------------------
   * // |       | 5 | 7 || -2
   * // |       |-------||---
   * // | 1     | 4 | 6 || -1
   * // --------------------
   * ////////////////////
   *
   * ~~~
   *
   * Top level dimension is 1.0cm. Lower left corner is at (0,0).
   */

  std::array<int, 28> neighbor_count = {1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  std::array<int, 31> neighbor_values = {1, 4, 5, 1, 2, 2, 3, 1, 0,  2, -3, -4, 5,  7, 0, 1,
                                         6, 4, 5, 1, 7, 4, 3, 4, -1, 6, 7,  5,  -2, 6, 3};

  std::array<int, 7> level_values = {0, 0, 0, 1, 1, 1, 1};

  double constexpr dimension = 1.0;

  std::array<double, 14> vertex_data = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
                                        0.0, 1.0, 0.5, 1.5, 0.0, 1.5, 0.5};

  CAMR_Mesh<2> const mesh(7, dimension, vertex_data.begin(), level_values.data(),
                          neighbor_count.data(), neighbor_values.data());

  // Test the equality operator
  FAIL_IF_NOT((mesh == mesh));
  FAIL_IF(mesh != mesh);
  FAIL_IF_NOT(soft_equiv(mesh.get_total_volume(), 4.0));
  FAIL_IF_NOT(mesh.full_Mesh());

  // ----------------------------------------------------------------
  // check that the mesh returns proper coordinate system information
  // ----------------------------------------------------------------
  FAIL_IF_NOT(mesh.get_spatial_dimension() == 2);

  // ----------------------------------------------
  // Check the vertices and dimensions of the cells
  // ----------------------------------------------
  {
    std::array<double, 2> dim = {1.0, 0.5};

    FAIL_IF_NOT(mesh.get_spatial_dimension() == 2);
    FAIL_IF_NOT(mesh.get_size() == 7);
    FAIL_IF_NOT(mesh.get_maximum_level() == 1);
    FAIL_IF_NOT(soft_equiv(mesh.get_width(), 1.0));
    FAIL_IF_NOT(soft_equiv(mesh.get_width(0), 1.0));
    FAIL_IF_NOT(soft_equiv(mesh.get_width(1), 0.5));
    FAIL_IF_NOT(
        soft_equiv(mesh.get_widths().begin(), mesh.get_widths().end(), dim.begin(), dim.end()));

    FAIL_IF_NOT(mesh.get_level(1) == 0);
    FAIL_IF_NOT(mesh.get_level(7) == 1);

    FAIL_IF_NOT(soft_equiv(mesh.get_cell_width(1), 1.0));
    FAIL_IF_NOT(soft_equiv(mesh.cell_width(1, 7), 0.5));
    FAIL_IF_NOT(soft_equiv(mesh.cell_width(2, 7), 0.5));

    // Check the vertices of _all_ of the cells, because I can:
    for (int cell = 1; cell <= 7; ++cell) {
      std::vector<double> const vertex = mesh.get_vertex(cell);
      FAIL_IF_NOT(soft_equiv(vertex.begin(), vertex.end(), vertex_data.begin() + 2 * (cell - 1),
                             vertex_data.begin() + 2 * cell));
    }
  }

  // ---------------
  // Test is_in_cell
  // ---------------
  {
    int cell;
    std::array<double, 3> position{0., 0., 0.};

    cell = 1;
    position[0] = 0.5;
    position[1] = 0.5;
    FAIL_IF(mesh.in_cell(cell, position) != true);

    position[0] = 1.0;
    FAIL_IF_NOT(mesh.in_cell(cell, position));
    FAIL_IF_NOT(mesh.hard_in_cell(cell, position));

    position[0] = 1.0 + 1.0e-14;
    FAIL_IF_NOT(mesh.in_cell(cell, position));
    FAIL_IF(mesh.hard_in_cell(cell, position));

    cell = 2;
    FAIL_IF(mesh.in_cell(cell, position) != false);

    FAIL_IF_NOT(soft_equiv(mesh.min(1, cell), 0.));
    FAIL_IF_NOT(soft_equiv(mesh.max(1, cell), 1.));

    std::array<double, 3> pos{0.1, 0., 0.};
    FAIL_IF_NOT(mesh.in_cell(1, 1, pos));
    pos[0] = 1.1;
    FAIL_IF_NOT(!mesh.in_cell(1, pos));

    FAIL_IF_NOT(soft_equiv(mesh.get_min_width(cell), 1.));
  }

  // ---------
  // Next cell:
  // ---------
  {
    std::array<double, 3> position{0., 0., 0.};

    position[0] = 1.0;
    position[1] = 0.75;
    FAIL_IF(mesh.next_cell(1, 2, position) != 5);

    position[0] = 1.0;
    position[1] = 1.5;
    FAIL_IF(mesh.next_cell(2, 2, position) != 3);

    position[0] = 1.25;
    position[1] = 1.0;
    FAIL_IF(mesh.next_cell(3, 3, position) != 5);

    position[0] = 2.0;
    position[1] = 1.75;
    FAIL_IF(mesh.next_cell(3, 2, position) != -4);

    position[0] = 1.5;
    position[1] = 2.0;
    FAIL_IF(mesh.next_cell(3, 4, position) != 0);

    position[0] = 1.25;
    position[1] = 0.0;
    FAIL_IF(mesh.next_cell(4, 3, position) != 4);
  }

  // ---------
  // Next face:
  // ---------
  {
    FAIL_IF(mesh.next_face(1, 2) != 1);
    FAIL_IF(mesh.next_face(1, 4) != 3);
    FAIL_IF(mesh.next_face(2, 2) != 1);
    FAIL_IF(mesh.next_face(2, 3) != 4);
    FAIL_IF(mesh.next_face(3, 1) != 2);
    FAIL_IF(mesh.next_face(3, 3) != 4);
    FAIL_IF(mesh.next_face(4, 1) != 2);
    FAIL_IF(mesh.next_face(4, 2) != 1);
    FAIL_IF(mesh.next_face(4, 4) != 3);
    FAIL_IF(mesh.next_face(5, 1) != 2);
    FAIL_IF(mesh.next_face(5, 2) != 1);
    FAIL_IF(mesh.next_face(5, 3) != 4);
    FAIL_IF(mesh.next_face(5, 4) != 3);
    FAIL_IF(mesh.next_face(6, 1) != 2);
    FAIL_IF(mesh.next_face(6, 2) != 1);
    FAIL_IF(mesh.next_face(6, 4) != 3);
    FAIL_IF(mesh.next_face(7, 1) != 2);
    FAIL_IF(mesh.next_face(7, 2) != 1);
    FAIL_IF(mesh.next_face(7, 3) != 4);
    FAIL_IF(mesh.next_face(7, 4) != 3);
    // For completeness, be sure you can't access a cell across an invalid face:
    bool caught(false);
    try {
      size_t const nf = mesh.next_face(7, 7);
      FAIL_IF(nf == 0); // should never get here.
    } catch (...) {
      caught = true;
    }
    FAIL_IF_NOT(caught);
  }

  // --------
  // Get cell:
  // --------
  {
    std::array<double, 3> position{0., 0., 0.};

    position[0] = 0.5;
    position[1] = 0.5;
    FAIL_IF(mesh.get_cell(position) != 1);

    position[0] = 1.25;
    position[1] = 0.75;
    FAIL_IF(mesh.get_cell(position) != 5);

    position[0] = -100;
    position[1] = 0.75;
    FAIL_IF(mesh.get_cell(position) != -1);
  }

  // --------
  // Get face:
  // --------
  {
    std::array<double, 3> position{0., 0., 0.};
    position[0] = 1.5;
    position[1] = 0.75;
    FAIL_IF(mesh.get_face(position, 7) != 1);

    position[0] = 2.0;
    FAIL_IF(mesh.get_face(position, 7) != 2);

    position[0] = 1.75;
    position[1] = 0.5;
    FAIL_IF(mesh.get_face(position, 7) != 3);

    position[1] = 1.0;
    FAIL_IF(mesh.get_face(position, 7) != 4);

    position[1] = 0.75;
    FAIL_IF(mesh.get_face(position, 7) != 0);
  }

  //-----------------
  // Get cell center:
  //-----------------
  {
    std::array<double, 3> center{0., 0., 0.};

    center = mesh.get_cell_centroid(1);
    FAIL_IF_NOT(soft_equiv(center[0], 0.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_cell_centroid(2);
    FAIL_IF_NOT(soft_equiv(center[0], 0.5));
    FAIL_IF_NOT(soft_equiv(center[1], 1.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_cell_centroid(3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 1.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_cell_centroid(4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_cell_centroid(5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_cell_centroid(6);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_cell_centroid(7);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));
  }

  //-----------------
  // Get face center:
  //-----------------
  {
    std::array<double, 3> center{0., 0., 0.};

    center = mesh.get_face_center(1, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 0.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(1, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(1, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 0.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(1, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 0.5));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(2, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 0.0));
    FAIL_IF_NOT(soft_equiv(center[1], 1.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(2, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 1.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(2, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 0.5));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(2, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 0.5));
    FAIL_IF_NOT(soft_equiv(center[1], 2.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(3, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 1.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(3, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 2.0));
    FAIL_IF_NOT(soft_equiv(center[1], 1.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(3, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(3, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 2.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(4, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(4, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(4, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(4, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(5, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(5, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(5, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(5, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(6, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(6, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 2.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(6, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(6, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(7, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(7, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 2.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(7, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(7, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));
  }

  if (ut.numFails == 0)
    PASSMSG("2D CAMR_Mesh test passed");
  return;
}

//--------------------------------------------------------------------------------------------//
void test_2D_mesh2(rtt_c4::ParallelUnitTest &ut) {
  /*
   * Create a simple 2D CAMR_Mesh
   *
   * ~~~
   *
   * // --------------------
   * // |       |       || -4
   * // |       |       ||---
   * // | 2     | 3     || -3
   * // ---------------------
   * // |       | 5 | 7 || -2
   * // |       |-------||---
   * // | 1     | 4 | 6 || -1
   * // --------------------
   * ////////////////////
   *
   * ~~~
   *
   * Top level dimension is 1.0cm. Lower left corner is at (0,0).
   */

  std::array<int, 28> neighbor_count = {1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  std::array<int, 31> neighbor_values = {1, 4, 5, 1, 2, 2, 3, 1, 0,  2, -3, -4, 5,  7, 0, 1,
                                         6, 4, 5, 1, 7, 4, 3, 4, -1, 6, 7,  5,  -2, 6, 3};

  std::array<int, 7> level_values = {0, 0, 0, 1, 1, 1, 1};

  double constexpr dimension = 1.0;

  std::array<double, 14> vertex_data = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
                                        0.0, 1.0, 0.5, 1.5, 0.0, 1.5, 0.5};

  CAMR_Mesh<2> const mesh(7, dimension, vertex_data.begin(), level_values.data(),
                          neighbor_count.data(), neighbor_values.data());

  // Test the equality operator
  FAIL_IF_NOT((mesh == mesh));
  FAIL_IF(mesh != mesh);
  FAIL_IF_NOT(soft_equiv(mesh.get_total_volume(), 4.0));
  FAIL_IF_NOT(mesh.full_Mesh());

  // ----------------------------------------------------------------
  // check that the mesh returns proper coordinate system information
  // ----------------------------------------------------------------
  FAIL_IF_NOT(mesh.get_spatial_dimension() == 2);

  //-----------------
  // Get shared face center:
  //-----------------
  {
    std::array<double, 3> center{0., 0., 0.};

    center = mesh.get_shared_face_center(1, 2, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(1, 2, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(1, 4, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 0.5));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(2, 2, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 1.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(2, 3, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 0.5));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(3, 1, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 1.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(3, 2, -3);
    FAIL_IF_NOT(soft_equiv(center[0], 2.0));
    FAIL_IF_NOT(soft_equiv(center[1], 1.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(3, 2, -4);
    FAIL_IF_NOT(soft_equiv(center[0], 2.0));
    FAIL_IF_NOT(soft_equiv(center[1], 1.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(3, 3, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(3, 3, 7);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(4, 1, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(4, 2, 6);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(4, 4, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(5, 1, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(5, 2, 7);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(5, 3, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(5, 4, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(6, 1, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(6, 2, -1);
    FAIL_IF_NOT(soft_equiv(center[0], 2.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(6, 4, 7);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(7, 1, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(7, 2, -2);
    FAIL_IF_NOT(soft_equiv(center[0], 2.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(7, 3, 6);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_shared_face_center(7, 4, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));
  }

  // ------------------------
  // Get information about refinement level
  // ------------------------
  {
    FAIL_IF(mesh.is_refined(1, 1));
    FAIL_IF_NOT(mesh.is_refined(1, 2));
    FAIL_IF(mesh.is_refined(1, 3));
    FAIL_IF(mesh.is_refined(1, 4));
    FAIL_IF(mesh.is_refined(3, 1));
    FAIL_IF_NOT(mesh.is_refined(3, 2));
    FAIL_IF_NOT(mesh.is_refined(3, 3));
    FAIL_IF(mesh.is_refined(3, 4));
  }

  // ------------------------
  // Get Distance to Boundary:
  // ------------------------
  {
    std::array<double, 3> position{0, 0.0}, direction{0, 0, 0};
    int cell;
    uint32_t face;

    position[0] = 0.0;
    position[1] = 0.5;
    direction[0] = 1.0;
    direction[1] = 0.0;
    direction[2] = 0.0;
    cell = 1;

    FAIL_IF_NOT(soft_equiv(mesh.get_db(position, direction, cell, face), 1.0));
    FAIL_IF_NOT(face == 2);

    position[0] = 1.5;
    position[1] = 0.5;
    direction[0] = sqrt(2.0) / 2.0;
    direction[1] = sqrt(2.0) / 2;
    direction[2] = 0.0;
    cell = 7;

    FAIL_IF_NOT(soft_equiv(mesh.get_db(position, direction, cell, face), sqrt(2.0) / 2.0));
    FAIL_IF(face != 2 && face != 4); // Allow both for a tie

    // Since the position is a corner. Try another cell.
    cell = 4;
    FAIL_IF_NOT(soft_equiv(mesh.get_db(position, direction, cell, face), 0.0));
    FAIL_IF(face != 2 && face != 4); // Allow both for a tie

    // Try an out-of-cell computation, where the point is not inside of the given cell. This is one
    // of the capabilities which provides robustness in the presence of roundoff.
    position[0] = 1.0;
    position[1] = 0.0;
    direction[0] = sqrt(2.0) / 2.0;
    direction[1] = sqrt(2.0) / 2.0;
    direction[2] = 0.0;
    cell = 7;

    FAIL_IF_NOT(soft_equiv(mesh.get_db(position, direction, cell, face), sqrt(2.0)));
    FAIL_IF(face != 2 && face != 4); // Allow both for a tie

    // Test the other direction.
    position[0] = 1.5;
    position[1] = 0.5;
    direction[0] = -1.0;
    direction[1] = 0.0;
    direction[2] = 0.0;
    cell = 4;

    FAIL_IF_NOT(soft_equiv(mesh.get_db(position, direction, cell, face), 0.5));
    FAIL_IF_NOT(face == 1);

    position[0] = 1.5;
    position[1] = 0.5;
    direction[0] = 0.0;
    direction[1] = -1.0;
    direction[2] = 0.0;
    cell = 4;

    FAIL_IF_NOT(soft_equiv(mesh.get_db(position, direction, cell, face), 0.5));
    FAIL_IF(face != 3);
  }

  // ------------------
  // Orthogonal distance to the boundary
  // ------------------
  {
    std::array<double, 3> position{0., 0., 0.};
    int cell;

    position[0] = 0.5;
    position[1] = 0.5;
    cell = 1;
    FAIL_IF_NOT(soft_equiv(mesh.get_orthogonal_dist_to_bnd(position, cell), 0.5));

    position[0] = 1.0;
    position[1] = 0.5;
    cell = 1;
    FAIL_IF_NOT(soft_equiv(mesh.get_orthogonal_dist_to_bnd(position, cell), 0.0));

    position[0] = 1.1;
    position[1] = 0.75;
    cell = 5;
    FAIL_IF_NOT(soft_equiv(mesh.get_orthogonal_dist_to_bnd(position, cell), 0.1));
  }

  // -----------
  // Get normals
  // -----------
  {
    std::array<double, 3> result{0, 0, 0};
    int face;

    result[0] = -1.0;
    result[1] = 0.0;
    face = 1;
    FAIL_IF(mesh.get_normal(1, face) != result);

    result[0] = 1.0;
    result[1] = 0.0;
    face = 2;
    FAIL_IF(mesh.get_normal(2, face) != result);

    result[0] = 0.0;
    result[1] = -1.0;
    face = 3;
    FAIL_IF(mesh.get_normal(3, face) != result);

    result[0] = 0.0;
    result[1] = 1.0;
    face = 4;
    FAIL_IF(mesh.get_normal(4, face) != result);

    result[0] = 1.0;
    result[1] = 0.0;
    face = 1;
    FAIL_IF(mesh.get_normal_in(1, face) != result);

    result[0] = -1.0;
    result[1] = 0.0;
    face = 2;
    FAIL_IF(mesh.get_normal_in(2, face) != result);

    result[0] = 0.0;
    result[1] = 1.0;
    face = 3;
    FAIL_IF(mesh.get_normal_in(3, face) != result);

    result[0] = 0.0;
    result[1] = -1.0;
    face = 4;
    FAIL_IF(mesh.get_normal_in(4, face) != result);
  }

  // ---------------------------
  // Cell volumes and face areas, chord length.
  // ---------------------------
  {
    int cell;

    cell = 1;
    FAIL_IF_NOT(soft_equiv(mesh.volume(cell), 1.0));
    for (int face = 1; face <= 4; ++face)
      FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 1.0));
    FAIL_IF_NOT(soft_equiv(mesh.mean_chord_length(1), 1.0));

    cell = 4;
    FAIL_IF_NOT(soft_equiv(mesh.volume(cell), 0.25));
    for (int face = 1; face <= 4; ++face)
      FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(soft_equiv(mesh.mean_chord_length(4), 0.5));
  }
  // -------------------------
  // Complete vertex listings:
  // -------------------------
  {
    int cell;
    vector<vector<double>> vertices;

    cell = 1;
    vertices = mesh.get_vertices(cell);
    FAIL_IF(vertices.size() != 2);
    FAIL_IF(vertices[0].size() != 4);
    FAIL_IF(vertices[1].size() != 4);

    FAIL_IF_NOT(soft_equiv(vertices[0][0], 0.0));
    FAIL_IF_NOT(soft_equiv(vertices[1][3], 1.0));

    cell = 7;
    vertices = mesh.get_vertices(cell);
    FAIL_IF_NOT(soft_equiv(vertices[0][1], 2.0));
    FAIL_IF_NOT(soft_equiv(vertices[1][1], 0.5));
    FAIL_IF_NOT(soft_equiv(vertices[1][3], 1.0));
  }

  // ---------------------
  // Face vertex listings:
  // ---------------------
  {
    int face, cell;
    vector<vector<double>> vertices;

    cell = 1;
    face = 2;
    vertices = mesh.get_vertices(cell, face);
    FAIL_IF(vertices.size() != 2);
    FAIL_IF(vertices[0].size() != 2);
    FAIL_IF(vertices[1].size() != 2);

    FAIL_IF_NOT(soft_equiv(vertices[0][0], 1.0));
    FAIL_IF_NOT(soft_equiv(vertices[1][0], 0.0));

    cell = 5;
    face = 3;
    vertices = mesh.get_vertices(cell, face);
    FAIL_IF(vertices.size() != 2);

    FAIL_IF_NOT(soft_equiv(vertices[0][0], 1.0));
    FAIL_IF_NOT(soft_equiv(vertices[1][0], 0.5));
    FAIL_IF_NOT(soft_equiv(vertices[0][1], 1.5));
    FAIL_IF_NOT(soft_equiv(vertices[1][1], 0.5));
  }

  // ---------------
  // Boundary Labels
  // ---------------
  {
    int constexpr cell = 0;

    FAIL_IF(mesh.get_bndface("lox", cell) != 1);
    FAIL_IF(mesh.get_bndface("hix", cell) != 2);
    FAIL_IF(mesh.get_bndface("loy", cell) != 3);
    FAIL_IF(mesh.get_bndface("hiy", cell) != 4);
  }

  // ---------------
  // Boundary Cells
  // ---------------
  {
    std::vector<int> surface_cells;
    std::vector<int> surface_cells_ref;

    surface_cells_ref.resize(2);
    surface_cells_ref[0] = 1;
    surface_cells_ref[1] = 2;
    surface_cells = mesh.get_surcells("lox");
    std::sort(surface_cells.begin(), surface_cells.end());
    FAIL_IF(surface_cells != surface_cells_ref);

    surface_cells_ref.resize(0);
    surface_cells = mesh.get_surcells("hix");
    FAIL_IF(surface_cells != surface_cells_ref);

    surface_cells_ref.resize(3);
    surface_cells_ref[0] = 1;
    surface_cells_ref[1] = 4;
    surface_cells_ref[2] = 6;
    surface_cells = mesh.get_surcells("loy");
    std::sort(surface_cells.begin(), surface_cells.end());
    FAIL_IF(surface_cells != surface_cells_ref);

    surface_cells_ref.resize(2);
    surface_cells_ref[0] = 2;
    surface_cells_ref[1] = 3;
    surface_cells = mesh.get_surcells("hiy");
    std::sort(surface_cells.begin(), surface_cells.end());
    FAIL_IF(surface_cells != surface_cells_ref);

    // Checking surface cells:
    // -----------------------

    surface_cells.resize(2);
    surface_cells[0] = 2;
    surface_cells[1] = 3;
    FAIL_IF(mesh.check_defined_surcells("hiy", surface_cells) != true);

    surface_cells[1] = 1;
    FAIL_IF(mesh.check_defined_surcells("hiy", surface_cells) != false);
    FAIL_IF(mesh.check_defined_surcells("lox", surface_cells) != false);
  }

  // ---------------------------
  // Check the list of neighbors
  // ---------------------------
  {
    int cell;
    std::vector<int> neighbors, neighbors_ref;

    cell = 1;
    neighbors = mesh.get_neighbors(cell);
    std::sort(neighbors.begin(), neighbors.end());
    std::array<int, 5> neighbors_of_one = {1, 1, 2, 4, 5};
    neighbors_ref.assign(neighbors_of_one.begin(), neighbors_of_one.end());
    FAIL_IF_NOT(neighbors == neighbors_ref);

    cell = 7;
    neighbors = mesh.get_neighbors(cell);
    std::sort(neighbors.begin(), neighbors.end());
    std::array<int, 4> neighbors_of_seven = {-2, 3, 5, 6};
    neighbors_ref.assign(neighbors_of_seven.begin(), neighbors_of_seven.end());
    FAIL_IF_NOT(neighbors == neighbors_ref);
  }

  // -----------------------
  // Test the dump functions
  // -----------------------
  {
    vector<int> cell_types(mesh.get_cell_types());

    FAIL_IF_NOT(cell_types.size() == mesh.get_size());

    FAIL_IF(std::find_if(cell_types.begin(), cell_types.end(), [](const int i) {
              return i != rtt_viz::four_node_quadrangle;
            }) != cell_types.end());

    vector<vector<double>> vertices = mesh.get_point_coord();
    vector<vector<int>> indices = mesh.get_cell_pair();

    FAIL_IF(vertices.size() != mesh.get_size() * 4);

    // Cell three, first corner:
    int index = indices[3 - 1][0] - 1;
    FAIL_IF(index != 8);

    FAIL_IF_NOT(soft_equiv(vertices[index][0], 1.0));
    FAIL_IF_NOT(soft_equiv(vertices[index][1], 1.0));

    // Cell 4, upper right corner.
    index = indices[4 - 1][2] - 1;
    FAIL_IF(index != 14);

    FAIL_IF_NOT(soft_equiv(vertices[index][0], 1.5));
    FAIL_IF_NOT(soft_equiv(vertices[index][1], 0.5));

    // Cell 7, second corner:
    index = indices[7 - 1][1] - 1;
    FAIL_IF(index != 25);

    FAIL_IF_NOT(soft_equiv(vertices[index][0], 2.0));
    FAIL_IF_NOT(soft_equiv(vertices[index][1], 0.5));
  }

  // Check the valid faces
  vector<int> valid_faces(mesh.get_valid_cell_faces());
  FAIL_IF(valid_faces.size() != 4);
  FAIL_IF(valid_faces[0] != 1);
  FAIL_IF(valid_faces[1] != 2);
  FAIL_IF(valid_faces[2] != 3);
  FAIL_IF(valid_faces[3] != 4);

  // Check the number of faces
  for (size_t cell = 1; cell < mesh.num_cells(); ++cell)
    FAIL_IF(mesh.get_number_of_faces() != 4);

  if (ut.numFails == 0)
    PASSMSG("2D CAMR_Mesh2 test passed");
  return;
}

//------------------------------------------------------------------------------------------------//
void test_3D_mesh(rtt_c4::ParallelUnitTest &ut) {
  /* Create a simple 3D mesh
   *
   *          ---------------------
   *         /|       | 7 | 9 || -4
   *        / |       |-------||---
   *       /  | 1     | 6 | 8 || -3
   *      /   ---------------------
   *     /   /   /   /   //  /   /
   *    /---/---/---/---//---   /
   *    |  /    | 3/| 5 || -2  /
   *    | /     |-/-----||--- /
   *    |/  1   |/2 | 4 || -1/
   *    ---------------------
   *
   *   x  z    Reflecting conditions in negative directions.
   *   | /     Vacuum in positive y and z. Processor boundary in
   *   |/      positive x.
   *   +--y
   */

  std::array<int, 54> neighbor_count = {1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  std::array<int, 57> neighbor_values = {
      1, 2, 3, 6, 7, 1, 0, 1, 0, // Cell 1
      // No more refined faces...
      1, 4, 2, 3, 2, 6,  // Cell 2
      1, 5, 2, 0, 3, 7,  // Cell 3
      2, -1, 4, 5, 4, 8, // Cell 4
      3, -2, 4, 0, 5, 9, // Cell 5
      1, 8, 6, 7, 2, 0,  // Cell 6
      1, 9, 6, 0, 3, 0,  // Cell 7
      6, -3, 8, 9, 4, 0, // Cell 8
      7, -4, 8, 0, 5, 0  // Cell 9
  };

  std::array<int, 9> level_values = {0, 1, 1, 1, 1, 1, 1, 1, 1};

  double constexpr width = 1.0;

  std::array<double, 27> vertex_data = {
      0.0, 0.0, 0.0, // Cell 1
      // Refined cells at z=0...
      1.0, 0.0, 0.0, // Cell 2
      1.0, 0.5, 0.0, // Cell 3
      1.5, 0.0, 0.0, // Cell 4
      1.5, 0.5, 0.0, // Cell 5
      // Refined cells at z=0.5...
      1.0, 0.0, 0.5, // Cell 6
      1.0, 0.5, 0.5, // Cell 7
      1.5, 0.0, 0.5, // Cell 8
      1.5, 0.5, 0.5, // Cell 9
  };

  CAMR_Mesh<3> const mesh(9, width, vertex_data.begin(), level_values.data(), neighbor_count.data(),
                          neighbor_values.data());

  FAIL_IF_NOT(soft_equiv(mesh.get_total_volume(), 2.0));

  // ----------------------------------------------------------------
  // check that the mesh returns proper coordinate system information
  // ----------------------------------------------------------------
  FAIL_IF(mesh.get_spatial_dimension() != 3);

  // ----------------------------------------------
  // Check the vertices and widths of the cells
  // ----------------------------------------------
  {
    std::array<double, 2> dim = {1.0, 0.5};

    FAIL_IF(mesh.get_spatial_dimension() != 3);
    FAIL_IF(mesh.get_size() != 9);
    FAIL_IF(mesh.get_maximum_level() != 1);
    FAIL_IF_NOT(soft_equiv(mesh.get_width(), 1.0));
    FAIL_IF_NOT(soft_equiv(mesh.get_width(0), 1.0));
    FAIL_IF_NOT(soft_equiv(mesh.get_width(1), 0.5));
    FAIL_IF_NOT(
        soft_equiv(mesh.get_widths().begin(), mesh.get_widths().end(), dim.begin(), dim.end()));

    for (int cell = 1; cell <= 9; ++cell)
      FAIL_IF_NOT(static_cast<int>(mesh.get_level(cell)) == level_values[cell - 1]);

    // Check the vertices of _all_ of the cells, because I can:
    for (int cell = 1; cell <= 9; ++cell) {
      std::vector<double> vertex = mesh.get_vertex(cell);

      FAIL_IF_NOT(soft_equiv(vertex.begin(), vertex.end(), vertex_data.begin() + 3 * (cell - 1),
                             vertex_data.begin() + 3 * cell));
    }
  }

  // ---------------
  // Test in_cell
  // ---------------
  {
    int cell;
    std::array<double, 3> position{0., 0., 0.};

    cell = 1;
    position[0] = 0.5;
    position[1] = 0.5;
    position[2] = 0.5;
    FAIL_IF(mesh.in_cell(cell, position) != true);

    position[1] = 0.0 + 1.0e-14;
    FAIL_IF_NOT(mesh.in_cell(cell, position));
    FAIL_IF_NOT(mesh.hard_in_cell(cell, position));

    position[1] = 0.0;
    FAIL_IF_NOT(mesh.in_cell(cell, position));
    FAIL_IF_NOT(mesh.hard_in_cell(cell, position));

    position[1] = 0.0 - 1.0e-14;
    FAIL_IF_NOT(mesh.in_cell(cell, position));
    FAIL_IF(mesh.hard_in_cell(cell, position));

    position[1] = 0.0 - 1.0e-8;
    FAIL_IF(mesh.in_cell(cell, position));
    FAIL_IF(mesh.hard_in_cell(cell, position));

    cell = 2;
    FAIL_IF(mesh.in_cell(cell, position) != false);
  }

  // ---------
  // Next cell:
  // ---------
  {
    int cell, face;
    std::array<double, 3> position{0., 0., 0.};

    cell = 1;
    face = 2;
    position[0] = 1.0;
    position[1] = 0.25;
    position[2] = 0.25;
    FAIL_IF(mesh.next_cell(cell, face, position) != 2);

    position[0] = 1.0;
    position[1] = 0.75;
    position[2] = 0.25;
    FAIL_IF(mesh.next_cell(cell, face, position) != 3);

    position[0] = 1.0;
    position[1] = 0.25;
    position[2] = 0.75;
    FAIL_IF(mesh.next_cell(cell, face, position) != 6);

    position[0] = 1.0;
    position[1] = 0.75;
    position[2] = 0.75;
    FAIL_IF(mesh.next_cell(cell, face, position) != 7);

    vector<int> next_cells = mesh.next_cells(cell, face);

    FAIL_IF(next_cells.size() != 4);
    FAIL_IF(next_cells[0] != 2);
    FAIL_IF(next_cells[1] != 3);
    FAIL_IF(next_cells[2] != 6);
    FAIL_IF(next_cells[3] != 7);

    vector<int> next_cells2 = mesh.next_cells(static_cast<size_t>(cell), static_cast<size_t>(face));
    FAIL_IF(next_cells2.size() != 4);
    FAIL_IF(next_cells2[0] != 2);
    FAIL_IF(next_cells2[1] != 3);
    FAIL_IF(next_cells2[2] != 6);
    FAIL_IF(next_cells2[3] != 7);
  }

  // ---------
  // Next face:
  // ---------
  {
    FAIL_IF(mesh.next_face(static_cast<size_t>(1), static_cast<size_t>(2)) != 1);
    FAIL_IF(mesh.next_face(2, 1) != 2);
    FAIL_IF(mesh.next_face(2, 2) != 1);
    FAIL_IF(mesh.next_face(2, 4) != 3);
    FAIL_IF(mesh.next_face(2, 6) != 5);
    FAIL_IF(mesh.next_face(3, 1) != 2);
    FAIL_IF(mesh.next_face(3, 2) != 1);
    FAIL_IF(mesh.next_face(3, 3) != 4);
    FAIL_IF(mesh.next_face(3, 6) != 5);
    FAIL_IF(mesh.next_face(4, 1) != 2);
    FAIL_IF(mesh.next_face(4, 2) != 1);
    FAIL_IF(mesh.next_face(4, 4) != 3);
    FAIL_IF(mesh.next_face(4, 6) != 5);
    FAIL_IF(mesh.next_face(5, 1) != 2);
    FAIL_IF(mesh.next_face(5, 2) != 1);
    FAIL_IF(mesh.next_face(5, 3) != 4);
    FAIL_IF(mesh.next_face(5, 6) != 5);
    FAIL_IF(mesh.next_face(6, 1) != 2);
    FAIL_IF(mesh.next_face(6, 2) != 1);
    FAIL_IF(mesh.next_face(6, 4) != 3);
    FAIL_IF(mesh.next_face(6, 5) != 6);
    FAIL_IF(mesh.next_face(7, 1) != 2);
    FAIL_IF(mesh.next_face(7, 2) != 1);
    FAIL_IF(mesh.next_face(7, 3) != 4);
    FAIL_IF(mesh.next_face(7, 5) != 6);
    FAIL_IF(mesh.next_face(8, 1) != 2);
    FAIL_IF(mesh.next_face(8, 2) != 1);
    FAIL_IF(mesh.next_face(8, 4) != 3);
    FAIL_IF(mesh.next_face(8, 5) != 6);
    FAIL_IF(mesh.next_face(9, 1) != 2);
    FAIL_IF(mesh.next_face(9, 2) != 1);
    FAIL_IF(mesh.next_face(9, 3) != 4);
    FAIL_IF(mesh.next_face(9, 5) != 6);
  }

  // --------
  // Get cell:
  // --------
  {
    std::array<double, 3> position{0., 0., 0.};

    position[0] = 0.5;
    position[1] = 0.5;
    position[2] = 0.5;
    FAIL_IF(mesh.get_cell(position) != 1);

    position[0] = 1.25;
    position[1] = 0.75;
    position[2] = 0.75;
    FAIL_IF(mesh.get_cell(position) != 7);

    position[0] = -100;
    position[1] = 0.75;
    position[2] = 0.0;
    FAIL_IF(mesh.get_cell(position) != -1);
  }

  // --------
  // Get face:
  // --------
  {
    std::array<double, 3> position{0., 0., 0.};
    position[0] = 1.5;
    position[1] = 0.75;
    position[2] = 0.75;
    FAIL_IF(mesh.get_face(position, 9) != 1);

    position[0] = 2.0;
    FAIL_IF(mesh.get_face(position, 9) != 2);

    position[0] = 1.75;
    position[1] = 0.5;
    FAIL_IF(mesh.get_face(position, 9) != 3);

    position[1] = 1.0;
    FAIL_IF(mesh.get_face(position, 9) != 4);

    position[1] = 0.75;
    position[2] = 0.5;
    FAIL_IF(mesh.get_face(position, 9) != 5);

    position[2] = 1.0;
    FAIL_IF(mesh.get_face(position, 9) != 6);

    position[2] = 0.75;
    FAIL_IF(mesh.get_face(position, 9) != 0);
  }

  if (ut.numFails == 0)
    PASSMSG("3D CAMR_Mesh test passed");
  return;
}

//------------------------------------------------------------------------------------------------//
void test_3D_mesh1(rtt_c4::ParallelUnitTest &ut) {
  /* Create a simple 3D mesh
   *
   *          ---------------------
   *         /|       | 7 | 9 || -4
   *        / |       |-------||---
   *       /  | 1     | 6 | 8 || -3
   *      /   ---------------------
   *     /   /   /   /   //  /   /
   *    /---/---/---/---//---   /
   *    |  /    | 3/| 5 || -2  /
   *    | /     |-/-----||--- /
   *    |/  1   |/2 | 4 || -1/
   *    ---------------------
   *
   *   x  z    Reflecting conditions in negative directions.
   *   | /     Vacuum in positive y and z. Processor boundary in
   *   |/      positive x.
   *   +--y
   */

  std::array<int, 54> neighbor_count = {1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  std::array<int, 57> neighbor_values = {
      1, 2, 3, 6, 7, 1, 0, 1, 0, // Cell 1
      // No more refined faces...
      1, 4, 2, 3, 2, 6,  // Cell 2
      1, 5, 2, 0, 3, 7,  // Cell 3
      2, -1, 4, 5, 4, 8, // Cell 4
      3, -2, 4, 0, 5, 9, // Cell 5
      1, 8, 6, 7, 2, 0,  // Cell 6
      1, 9, 6, 0, 3, 0,  // Cell 7
      6, -3, 8, 9, 4, 0, // Cell 8
      7, -4, 8, 0, 5, 0  // Cell 9
  };

  std::array<int, 9> level_values = {0, 1, 1, 1, 1, 1, 1, 1, 1};

  double constexpr width = 1.0;

  std::array<double, 27> vertex_data = {
      0.0, 0.0, 0.0, // Cell 1
      // Refined cells at z=0...
      1.0, 0.0, 0.0, // Cell 2
      1.0, 0.5, 0.0, // Cell 3
      1.5, 0.0, 0.0, // Cell 4
      1.5, 0.5, 0.0, // Cell 5
      // Refined cells at z=0.5...
      1.0, 0.0, 0.5, // Cell 6
      1.0, 0.5, 0.5, // Cell 7
      1.5, 0.0, 0.5, // Cell 8
      1.5, 0.5, 0.5, // Cell 9
  };

  CAMR_Mesh<3> const mesh(9, width, vertex_data.begin(), level_values.data(), neighbor_count.data(),
                          neighbor_values.data());

  FAIL_IF_NOT(soft_equiv(mesh.get_total_volume(), 2.0));

  // ----------------------------------------------------------------
  // check that the mesh returns proper coordinate system information
  // ----------------------------------------------------------------
  FAIL_IF_NOT(mesh.get_spatial_dimension() == 3);

  //-----------------
  // Get face center:
  //-----------------
  {
    std::array<double, 3> center{0., 0., 0.};

    center = mesh.get_face_center(1, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 0.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_face_center(1, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_face_center(1, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 0.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_face_center(1, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 0.5));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_face_center(1, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 0.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(1, 6);
    FAIL_IF_NOT(soft_equiv(center[0], 0.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 1.0));

    center = mesh.get_face_center(2, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(2, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(2, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(2, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(2, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(2, 6);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_face_center(3, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(3, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(3, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(3, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(3, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(3, 6);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_face_center(4, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(4, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 2.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(4, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(4, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(4, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(4, 6);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_face_center(5, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(5, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 2.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(5, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(5, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_face_center(5, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.0));

    center = mesh.get_face_center(5, 6);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_face_center(6, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(6, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(6, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(6, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(6, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_face_center(6, 6);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 1.0));

    center = mesh.get_face_center(7, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(7, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(7, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(7, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(7, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_face_center(7, 6);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 1.0));

    center = mesh.get_face_center(8, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(8, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 2.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(8, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(8, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(8, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_face_center(8, 6);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 1.0));

    center = mesh.get_face_center(9, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(9, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 2.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(9, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(9, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 1.0));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_face_center(9, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_face_center(9, 6);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 1.0));
  }
  if (ut.numFails == 0)
    PASSMSG("3D CAMR_Mesh1 test passed");
  return;
}

//------------------------------------------------------------------------------------------------//
void test_3D_mesh2(rtt_c4::ParallelUnitTest &ut) {
  /* Create a simple 3D mesh
   *
   * ~~~
   *
   *          ---------------------
   *         /|       | 7 | 9 || -4
   *        / |       |-------||---
   *       /  | 1     | 6 | 8 || -3
   *      /   ---------------------
   *     /   /   /   /   //  /   /
   *    /---/---/---/---//---   /
   *    |  /    | 3/| 5 || -2  /
   *    | /     |-/-----||--- /
   *    |/  1   |/2 | 4 || -1/
   *    ---------------------
   *
   *   x  z    Reflecting conditions in negative directions.
   *   | /     Vacuum in positive y and z. Processor boundary in
   *   |/      positive x.
   *   +--y
   *
   * ~~~
   */

  std::array<int, 54> neighbor_count = {1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  std::array<int, 57> neighbor_values = {
      1, 2, 3, 6, 7, 1, 0, 1, 0, // Cell 1
      // No more refined faces...
      1, 4, 2, 3, 2, 6,  // Cell 2
      1, 5, 2, 0, 3, 7,  // Cell 3
      2, -1, 4, 5, 4, 8, // Cell 4
      3, -2, 4, 0, 5, 9, // Cell 5
      1, 8, 6, 7, 2, 0,  // Cell 6
      1, 9, 6, 0, 3, 0,  // Cell 7
      6, -3, 8, 9, 4, 0, // Cell 8
      7, -4, 8, 0, 5, 0  // Cell 9
  };

  std::array<int, 9> level_values = {0, 1, 1, 1, 1, 1, 1, 1, 1};

  double constexpr width = 1.0;

  std::array<double, 27> vertex_data = {
      0.0, 0.0, 0.0, // Cell 1
      // Refined cells at z=0...
      1.0, 0.0, 0.0, // Cell 2
      1.0, 0.5, 0.0, // Cell 3
      1.5, 0.0, 0.0, // Cell 4
      1.5, 0.5, 0.0, // Cell 5
      // Refined cells at z=0.5...
      1.0, 0.0, 0.5, // Cell 6
      1.0, 0.5, 0.5, // Cell 7
      1.5, 0.0, 0.5, // Cell 8
      1.5, 0.5, 0.5, // Cell 9
  };

  CAMR_Mesh<3> const mesh(9, width, vertex_data.begin(), level_values.data(), neighbor_count.data(),
                          neighbor_values.data());

  FAIL_IF_NOT(soft_equiv(mesh.get_total_volume(), 2.0));

  // ----------------------------------------------------------------
  // check that the mesh returns proper coordinate system information
  // ----------------------------------------------------------------
  FAIL_IF_NOT(mesh.get_spatial_dimension() == 3);

  //------------------------
  // Get shared face center:
  //------------------------
  {
    std::array<double, 3> center{0., 0., 0.};

    center = mesh.get_shared_face_center(1, 2, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_shared_face_center(1, 2, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_shared_face_center(1, 2, 6);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_shared_face_center(1, 2, 7);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_shared_face_center(2, 1, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_shared_face_center(2, 2, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_shared_face_center(2, 4, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_shared_face_center(2, 6, 6);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_shared_face_center(3, 1, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_shared_face_center(3, 2, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_shared_face_center(3, 3, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_shared_face_center(3, 6, 7);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_shared_face_center(4, 1, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_shared_face_center(4, 2, -1);
    FAIL_IF_NOT(soft_equiv(center[0], 2.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_shared_face_center(4, 4, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_shared_face_center(4, 6, 8);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_shared_face_center(5, 1, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_shared_face_center(5, 2, -2);
    FAIL_IF_NOT(soft_equiv(center[0], 2.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_shared_face_center(5, 3, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.25));

    center = mesh.get_shared_face_center(5, 6, 9);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_shared_face_center(6, 1, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_shared_face_center(6, 2, 8);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_shared_face_center(6, 4, 7);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_shared_face_center(6, 5, 2);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_shared_face_center(7, 1, 1);
    FAIL_IF_NOT(soft_equiv(center[0], 1.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_shared_face_center(7, 2, 9);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_shared_face_center(7, 3, 6);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_shared_face_center(7, 5, 3);
    FAIL_IF_NOT(soft_equiv(center[0], 1.25));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_shared_face_center(8, 1, 6);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_shared_face_center(8, 2, -3);
    FAIL_IF_NOT(soft_equiv(center[0], 2.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_shared_face_center(8, 4, 9);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_shared_face_center(8, 5, 4);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.25));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));

    center = mesh.get_shared_face_center(9, 1, 7);
    FAIL_IF_NOT(soft_equiv(center[0], 1.5));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_shared_face_center(9, 2, -4);
    FAIL_IF_NOT(soft_equiv(center[0], 2.0));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_shared_face_center(9, 3, 8);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.5));
    FAIL_IF_NOT(soft_equiv(center[2], 0.75));

    center = mesh.get_shared_face_center(9, 5, 5);
    FAIL_IF_NOT(soft_equiv(center[0], 1.75));
    FAIL_IF_NOT(soft_equiv(center[1], 0.75));
    FAIL_IF_NOT(soft_equiv(center[2], 0.5));
  }

  if (ut.numFails == 0)
    PASSMSG("3D CAMR_Mesh2 test passed");
  return;
}
//------------------------------------------------------------------------------------------------//
void test_3D_mesh3(rtt_c4::ParallelUnitTest &ut) {
  /* Create a simple 3D mesh
   *
   * ~~~
   *          ---------------------
   *         /|       | 7 | 9 || -4
   *        / |       |-------||---
   *       /  | 1     | 6 | 8 || -3
   *      /   ---------------------
   *     /   /   /   /   //  /   /
   *    /---/---/---/---//---   /
   *    |  /    | 3/| 5 || -2  /
   *    | /     |-/-----||--- /
   *    |/  1   |/2 | 4 || -1/
   *    ---------------------
   *
   *   x  z    Reflecting conditions in negative directions.
   *   | /     Vacuum in positive y and z. Processor boundary in
   *   |/      positive x.
   *   +--y
   * ~~~
   */

  std::array<int, 54> neighbor_count = {1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  std::array<int, 57> neighbor_values = {
      1, 2, 3, 6, 7, 1, 0, 1, 0, // Cell 1
      // No more refined faces...
      1, 4, 2, 3, 2, 6,  // Cell 2
      1, 5, 2, 0, 3, 7,  // Cell 3
      2, -1, 4, 5, 4, 8, // Cell 4
      3, -2, 4, 0, 5, 9, // Cell 5
      1, 8, 6, 7, 2, 0,  // Cell 6
      1, 9, 6, 0, 3, 0,  // Cell 7
      6, -3, 8, 9, 4, 0, // Cell 8
      7, -4, 8, 0, 5, 0  // Cell 9
  };

  std::array<int, 9> level_values = {0, 1, 1, 1, 1, 1, 1, 1, 1};

  double constexpr width = 1.0;

  std::array<double, 27> vertex_data = {
      0.0, 0.0, 0.0, // Cell 1
      // Refined cells at z=0...
      1.0, 0.0, 0.0, // Cell 2
      1.0, 0.5, 0.0, // Cell 3
      1.5, 0.0, 0.0, // Cell 4
      1.5, 0.5, 0.0, // Cell 5
      // Refined cells at z=0.5...
      1.0, 0.0, 0.5, // Cell 6
      1.0, 0.5, 0.5, // Cell 7
      1.5, 0.0, 0.5, // Cell 8
      1.5, 0.5, 0.5, // Cell 9
  };

  CAMR_Mesh<3> const mesh(9, width, vertex_data.begin(), level_values.data(), neighbor_count.data(),
                          neighbor_values.data());

  FAIL_IF_NOT(soft_equiv(mesh.get_total_volume(), 2.0));

  // ----------------------------------------------------------------
  // check that the mesh returns proper coordinate system information
  // ----------------------------------------------------------------
  FAIL_IF_NOT(mesh.get_spatial_dimension() == 3);

  // ------------------------
  // Get Distance to Boundary:
  // ------------------------
  {
    std::array<double, 3> position{0, 0, 0}, direction{0, 0, 0};
    uint32_t face;
    int cell;

    position[0] = 0.0;
    position[1] = 0.5;
    position[2] = 0.5;
    direction[0] = 1.0;
    direction[1] = 0.0;
    direction[2] = 0.0;
    cell = 1;

    FAIL_IF_NOT(soft_equiv(mesh.get_db(position, direction, cell, face), 1.0));
    FAIL_IF_NOT(face == 2);

    position[0] = 1.5;
    position[1] = 0.5;
    position[2] = 0.0;
    direction[0] = sqrt(2.0) / 2.0;
    direction[1] = sqrt(2.0) / 2.0;
    direction[2] = 0.0;
    cell = 5;

    FAIL_IF_NOT(soft_equiv(mesh.get_db(position, direction, cell, face), sqrt(2.0) / 2.0));
    FAIL_IF(face != 2 && face != 4); // Allow both for a tie.

    // Since the position is a corner. Try another cell.
    cell = 2;
    FAIL_IF_NOT(soft_equiv(mesh.get_db(position, direction, cell, face), 0.0));
    FAIL_IF(face != 2 && face != 4); // Allow both for a tie.

    // Try an out-of-cell computation, where the point is not inside of the given cell. This is one
    // of the capabilities which provides robustness in the presence of roundoff.
    position[0] = 1.0;
    position[1] = 0.0;
    position[2] = 0.0;
    direction[0] = sqrt(3.0) / 3.0;
    direction[1] = sqrt(3.0) / 3.0;
    direction[2] = sqrt(3.0) / 3.0;
    cell = 9;

    FAIL_IF_NOT(soft_equiv(mesh.get_db(position, direction, cell, face), sqrt(3.0)));
    // Allow all three in a tie.
    FAIL_IF(face != 2 && face != 4 && face != 6);

    // Test the other direction.
    position[0] = 1.5;
    position[1] = 0.5;
    position[2] = 0.5;
    direction[0] = -1.0;
    direction[1] = 0.0;
    direction[2] = 0.0;
    cell = 2;

    FAIL_IF_NOT(soft_equiv(mesh.get_db(position, direction, cell, face), 0.5));
    FAIL_IF(face != 1);

    position[0] = 1.5;
    position[1] = 0.5;
    position[2] = 0.5;
    direction[0] = 0.0;
    direction[1] = -1.0;
    direction[2] = 0.0;
    cell = 2;

    FAIL_IF_NOT(soft_equiv(mesh.get_db(position, direction, cell, face), 0.5));
    FAIL_IF(face != 3);
  }

  // ------------------
  // Orthogonal distance to boundary:
  // ------------------
  {
    std::array<double, 3> position{0., 0., 0.};
    int cell;

    position[0] = 0.5;
    position[1] = 0.5;
    position[2] = 0.5;
    cell = 1;
    FAIL_IF_NOT(soft_equiv(mesh.get_orthogonal_dist_to_bnd(position, cell), 0.5));

    position[0] = 1.0;
    position[1] = 0.5;
    position[2] = 0.75;
    cell = 1;
    FAIL_IF_NOT(soft_equiv(mesh.get_orthogonal_dist_to_bnd(position, cell), 0.0));

    position[0] = 1.1;
    position[1] = 0.75;
    position[2] = 0.25;
    cell = 3;
    FAIL_IF_NOT(soft_equiv(mesh.get_orthogonal_dist_to_bnd(position, cell), 0.1));
  }

  // -----------
  // Get normals
  // -----------
  {
    std::array<double, 3> result{0., 0., 0.};
    int constexpr cell = 1;
    int face;

    result[0] = -1.0;
    result[1] = 0.0;
    result[2] = 0.0;
    face = 1;
    FAIL_IF(mesh.get_normal(cell, face) != result);

    result[0] = 1.0;
    result[1] = 0.0;
    result[2] = 0.0;
    face = 2;
    FAIL_IF(mesh.get_normal(cell, face) != result);

    result[0] = 0.0;
    result[1] = -1.0;
    result[2] = 0.0;
    face = 3;
    FAIL_IF(mesh.get_normal(cell, face) != result);

    result[0] = 0.0;
    result[1] = 1.0;
    result[2] = 0.0;
    face = 4;
    FAIL_IF(mesh.get_normal(cell, face) != result);

    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = -1.0;
    face = 5;
    FAIL_IF(mesh.get_normal(cell, face) != result);

    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = 1.0;
    face = 6;
    FAIL_IF(mesh.get_normal(cell, face) != result);

    result[0] = 1.0;
    result[1] = 0.0;
    result[2] = 0.0;
    face = 1;
    FAIL_IF(mesh.get_normal_in(cell, face) != result);

    result[0] = -1.0;
    result[1] = 0.0;
    result[2] = 0.0;
    face = 2;
    FAIL_IF(mesh.get_normal_in(cell, face) != result);

    result[0] = 0.0;
    result[1] = 1.0;
    result[2] = 0.0;
    face = 3;
    FAIL_IF(mesh.get_normal_in(cell, face) != result);

    result[0] = 0.0;
    result[1] = -1.0;
    result[2] = 0.0;
    face = 4;
    FAIL_IF(mesh.get_normal_in(cell, face) != result);

    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = 1.0;
    face = 5;
    FAIL_IF(mesh.get_normal_in(cell, face) != result);

    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = -1.0;
    face = 6;
    FAIL_IF(mesh.get_normal_in(cell, face) != result);
  }

  // ---------------------------
  // Cell volumes and face areas.
  // ---------------------------
  {
    int cell;

    cell = 1;
    FAIL_IF_NOT(soft_equiv(mesh.volume(cell), 1.0));
    for (int face = 1; face <= 6; ++face)
      FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 1.0));

    FAIL_IF_NOT(soft_equiv(mesh.mean_chord_length(cell), (2.0 / 3.0)));

    cell = 4;
    FAIL_IF_NOT(soft_equiv(mesh.volume(cell), 0.125));
    for (int face = 1; face <= 6; ++face)
      FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));

    FAIL_IF_NOT(soft_equiv(mesh.mean_chord_length(cell), (1.0 / 3.0)));
  }

  // -------------------------
  // Complete vertex listings:
  // -------------------------
  {
    int cell;
    vector<vector<double>> vertices;

    cell = 1;
    vertices = mesh.get_vertices(cell);
    FAIL_IF(vertices.size() != 3);
    FAIL_IF(vertices[0].size() != 8);
    FAIL_IF(vertices[1].size() != 8);
    FAIL_IF(vertices[2].size() != 8);

    FAIL_IF_NOT(soft_equiv(vertices[0][0], 0.0));
    FAIL_IF_NOT(soft_equiv(vertices[1][3], 1.0));

    cell = 9;
    vertices = mesh.get_vertices(cell);
    FAIL_IF_NOT(soft_equiv(vertices[0][1], 2.0));
    FAIL_IF_NOT(soft_equiv(vertices[1][1], 0.5));
    FAIL_IF_NOT(soft_equiv(vertices[2][3], 0.5));
  }

  // ---------------------
  // Face vertex listings:
  // ---------------------
  {
    int face, cell;
    vector<vector<double>> vertices;

    cell = 1;
    face = 2;
    vertices = mesh.get_vertices(cell, face);
    FAIL_IF(vertices.size() != 3);
    FAIL_IF(vertices[0].size() != 4);
    FAIL_IF(vertices[1].size() != 4);

    FAIL_IF_NOT(soft_equiv(vertices[0][0], 1.0));
    FAIL_IF_NOT(soft_equiv(vertices[1][0], 0.0));
    FAIL_IF_NOT(soft_equiv(vertices[2][0], 0.0));

    cell = 9;
    face = 3;
    vertices = mesh.get_vertices(cell, face);
    FAIL_IF(vertices.size() != 3);

    FAIL_IF_NOT(soft_equiv(vertices[0][0], 1.5));
    FAIL_IF_NOT(soft_equiv(vertices[1][0], 0.5));
    FAIL_IF_NOT(soft_equiv(vertices[2][0], 0.5));

    FAIL_IF_NOT(soft_equiv(vertices[0][3], 2.0));
    FAIL_IF_NOT(soft_equiv(vertices[1][3], 0.5));
    FAIL_IF_NOT(soft_equiv(vertices[2][3], 1.0));
  }

  // ---------------
  // Boundary Labels
  // ---------------
  {
    int constexpr cell = 0;

    FAIL_IF(mesh.get_bndface("lox", cell) != 1);
    FAIL_IF(mesh.get_bndface("hix", cell) != 2);
    FAIL_IF(mesh.get_bndface("loy", cell) != 3);
    FAIL_IF(mesh.get_bndface("hiy", cell) != 4);
    FAIL_IF(mesh.get_bndface("loz", cell) != 5);
    FAIL_IF(mesh.get_bndface("hiz", cell) != 6);
  }

  // ---------------
  // Boundary Cells
  // ---------------
  {
    std::vector<int> surface_cells;
    std::vector<int> surface_cells_ref;

    surface_cells_ref.resize(1);
    surface_cells_ref[0] = 1;
    surface_cells = mesh.get_surcells("lox");
    std::sort(surface_cells.begin(), surface_cells.end());
    FAIL_IF(surface_cells != surface_cells_ref);

    surface_cells_ref.resize(0);
    surface_cells = mesh.get_surcells("hix");
    FAIL_IF(surface_cells != surface_cells_ref);

    surface_cells_ref.resize(5);
    surface_cells_ref[0] = 1;
    surface_cells_ref[1] = 2;
    surface_cells_ref[2] = 4;
    surface_cells_ref[3] = 6;
    surface_cells_ref[4] = 8;
    surface_cells = mesh.get_surcells("loy");
    std::sort(surface_cells.begin(), surface_cells.end());
    FAIL_IF(surface_cells != surface_cells_ref);

    surface_cells_ref.resize(5);
    surface_cells_ref[0] = 1;
    surface_cells_ref[1] = 3;
    surface_cells_ref[2] = 5;
    surface_cells_ref[3] = 7;
    surface_cells_ref[4] = 9;
    surface_cells = mesh.get_surcells("hiy");
    std::sort(surface_cells.begin(), surface_cells.end());
    FAIL_IF(surface_cells != surface_cells_ref);

    surface_cells_ref.resize(5);
    surface_cells_ref[0] = 1;
    surface_cells_ref[1] = 2;
    surface_cells_ref[2] = 3;
    surface_cells_ref[3] = 4;
    surface_cells_ref[4] = 5;
    surface_cells = mesh.get_surcells("loz");
    std::sort(surface_cells.begin(), surface_cells.end());
    FAIL_IF(surface_cells != surface_cells_ref);

    surface_cells_ref.resize(5);
    surface_cells_ref[0] = 1;
    surface_cells_ref[1] = 6;
    surface_cells_ref[2] = 7;
    surface_cells_ref[3] = 8;
    surface_cells_ref[4] = 9;
    surface_cells = mesh.get_surcells("hiz");
    std::sort(surface_cells.begin(), surface_cells.end());
    FAIL_IF(surface_cells != surface_cells_ref);

    // Checking surface cells:
    // -----------------------

    surface_cells.resize(3);
    surface_cells[0] = 1;
    surface_cells[1] = 3;
    surface_cells[2] = 9;
    FAIL_IF(mesh.check_defined_surcells("hiy", surface_cells) != true);

    surface_cells[2] = 4;
    FAIL_IF(mesh.check_defined_surcells("loz", surface_cells) != false);
    FAIL_IF(mesh.check_defined_surcells("hiy", surface_cells) != false);
    FAIL_IF(mesh.check_defined_surcells("lox", surface_cells) != false);
  }

  // ---------------------------
  // Check the list of neighbors
  // ---------------------------
  {
    int cell;
    std::vector<int> neighbors, neighbors_ref;

    cell = 1;
    neighbors = mesh.get_neighbors(cell);
    std::sort(neighbors.begin(), neighbors.end());
    std::array<int, 9> neighbors_of_one = {0, 0, 1, 1, 1, 2, 3, 6, 7};
    neighbors_ref.assign(neighbors_of_one.begin(), neighbors_of_one.end());
    FAIL_IF_NOT(neighbors == neighbors_ref);

    cell = 9;
    neighbors = mesh.get_neighbors(cell);
    std::sort(neighbors.begin(), neighbors.end());
    std::array<int, 6> neighbors_of_nine = {-4, 0, 0, 5, 7, 8};
    neighbors_ref.assign(neighbors_of_nine.begin(), neighbors_of_nine.end());
    FAIL_IF_NOT(neighbors == neighbors_ref);
  }

  // -----------------------
  // Test the dump functions
  // -----------------------
  {
    vector<int> cell_types(mesh.get_cell_types());

    FAIL_IF_NOT(cell_types.size() == mesh.get_size());

    FAIL_IF(std::find_if(cell_types.begin(), cell_types.end(), [](const int i) {
              return i != rtt_viz::eight_node_hexahedron;
            }) != cell_types.end());

    vector<vector<double>> vertices = mesh.get_point_coord();
    vector<vector<int>> indices = mesh.get_cell_pair();

    FAIL_IF(vertices.size() != mesh.get_size() * 8);

    // Cell 1, first vertex:
    int index = indices[1 - 1][0] - 1;
    FAIL_IF(index != 0);

    FAIL_IF_NOT(soft_equiv(vertices[index][0], 0.0));
    FAIL_IF_NOT(soft_equiv(vertices[index][1], 0.0));
    FAIL_IF_NOT(soft_equiv(vertices[index][2], 0.0));

    // Cell 6, fourth vertex:
    index = indices[6 - 1][3] - 1;
    FAIL_IF(index != 43);

    FAIL_IF_NOT(soft_equiv(vertices[index][0], 1.0));
    FAIL_IF_NOT(soft_equiv(vertices[index][1], 0.5));
    FAIL_IF_NOT(soft_equiv(vertices[index][2], 0.5));

    // Cell 9, 7th vertex:
    index = indices[9 - 1][6] - 1;
    FAIL_IF(index != 70);

    FAIL_IF_NOT(soft_equiv(vertices[index][0], 2.0));
    FAIL_IF_NOT(soft_equiv(vertices[index][1], 1.0));
    FAIL_IF_NOT(soft_equiv(vertices[index][2], 1.0));
  }

  // Check the valid faces
  vector<int> valid_faces(mesh.get_valid_cell_faces());
  FAIL_IF(valid_faces.size() != 6);
  FAIL_IF(valid_faces[0] != 1);
  FAIL_IF(valid_faces[1] != 2);
  FAIL_IF(valid_faces[2] != 3);
  FAIL_IF(valid_faces[3] != 4);
  FAIL_IF(valid_faces[4] != 5);
  FAIL_IF(valid_faces[5] != 6);

  // Check the number of faces
  for (size_t cell = 1; cell < mesh.num_cells(); ++cell)
    FAIL_IF_NOT(mesh.get_number_of_faces() == 6);

  if (ut.numFails == 0)
    PASSMSG("3D CAMR_Mesh3 test passed");
  return;
}

//--------------------------------------------------//
void test_2D_host(rtt_c4::ParallelUnitTest &ut) {
  /*
   * Create a simple 2D CAMR_Mesh
   *
   * ~~~
   * // --------------------
   * // |       |       || -4
   * // |       |       ||---
   * // | 2     | 3     || -3
   * // ---------------------
   * // |       | 5 | 7 || -2
   * // |       |-------||---
   * // | 1     | 4 | 6 || -1
   * // --------------------
   * ////////////////////
   * ~~~
   * Top level dimension is 1.0cm. Lower left corner is at (0,0).
   */

  std::array<int, 28> neighbor_count = {1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  std::array<int, 31> neighbor_values = {1, 4, 5, 1, 2, 2, 3, 1, 0,  2, -3, -4, 5,  7, 0, 1,
                                         6, 4, 5, 1, 7, 4, 3, 4, -1, 6, 7,  5,  -2, 6, 3};

  std::array<int, 7> level_values = {0, 0, 0, 1, 1, 1, 1};

  double const dimension = 1.0;

  std::array<double, 14> vertex_data = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
                                        0.0, 1.0, 0.5, 1.5, 0.0, 1.5, 0.5};

  CAMR_Mesh_Data md;

  md.num_across_face = neighbor_count.data();
  md.xyz_extents = vertex_data.data();
  md.cell_levels = level_values.data();
  md.cell_across_face = neighbor_values.data();

  CAMR_Mesh<2> const mesh(7, dimension, md.xyz_extents, md.cell_levels, md.num_across_face,
                          md.cell_across_face);

  // mesh.get_Layout().print(std::cout, 1);

  CAMR_Layout const &clo = mesh.get_Layout();
  std::cout << clo;
  {
    CAMR_Layout const &clo2 = mesh.get_Layout();
    FAIL_IF_NOT(clo == clo2);
  }

  // Test the equality operator
  FAIL_IF_NOT((mesh == mesh));
  FAIL_IF(mesh != mesh);
  FAIL_IF_NOT(soft_equiv(mesh.get_total_volume(), 4.0));
  FAIL_IF_NOT(mesh.full_Mesh());

  size_t const nc = mesh.get_num_full_cells();
  size_t const ng = mesh.get_num_ghost_cells();
  FAIL_IF_NOT(nc == 7);
  FAIL_IF_NOT(ng == 0);

  if (ut.numFails == 0)
    PASSMSG("2D CAMR_host test passed");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_2D_mesh(ut);
    test_2D_mesh2(ut);
    test_3D_mesh(ut);
    test_3D_mesh1(ut);
    test_3D_mesh2(ut);
    test_3D_mesh3(ut);

    test_2D_host(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstRage_Mesh.cc
//------------------------------------------------------------------------------------------------//
