//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/test/tstCAMR_Mesh.cc
 * \author H. Park
 * \date   Mon Jul 10 2023
 * \brief  Unit test for CAMR_Mesh
 * \note   Copyright (C) 2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "mesh/CAMR_Mesh.hh"
#include "viz/Ensight_Translator.hh"

using namespace std;
using namespace rtt_mesh;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
struct CAMR_Mesh_Data {
  //! RZ data and AMR connectivities
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

//--------------------------------------------------//
void test_2D_host(rtt_c4::ParallelUnitTest &ut) {
  /*
--
-- Create a simple 2D Rage_Mesh (test for domian decomposed,
--  non-orthgonal dd boundary
--   * 
--   * // ------------------
--   * // |       |        |
--   * // |       | -1     |
--   * // | 3     |        |
--   * // ------------------
--   * // |       |        |
--   * // |       |        |
--   * // | 1     | 2      |
--   * // ------------------
--   * ////////////////////
--   *
--   * // ------------------
--   * // |       |        |
--   * // |       | 1      |
--   * // | -2    |        |
--   * // ------------------
--   * // |       |        |
--   * // |       |        |
--   * // |       |-1      |
--   * // ------------------
--   * ////////////////////
--   *
-- mesh sizes
--NOTE about cell_across_face indexing
--We assume (in 2D ordering of  mesh.next_cells() to be (lower, upper)
--e.g., see the pic below
--    !  RZ:
--    !
--    !      | 1 | 2 |
--    !      ---------
--    !    2 |       | 2
--    !  R --|       |--
--    !    1 |       | 1
--    !      ---------
--    !      | 1 | 2 |
--    !
--    !          Z
   * Top level dimension is 1.0cm. Lower left corner is at (0,0).
   */

  std::cout << "in dd host" << std::endl;

  CAMR_Mesh_Data md;

  //    size_t full_size=4;
  size_t ghost_size = 0;

  double dx = 0.5;
  double x0 = 0;
  double y0 = 0;

  std::vector<int> full2active;
  std::vector<int> active2full;

  std::vector<double> vertex_data;
  std::vector<int> neighbor_count;
  std::vector<int> level_values;
  std::vector<int> neighbor_values;
  std::vector<int> bnd_cells;
  std::vector<int> bnd_procs;

  if (rtt_c4::node() == 0) {
    full2active.resize(3);
    full2active[0] = 1;
    full2active[1] = 2;
    full2active[2] = 3;
    active2full.resize(3);
    active2full[0] = 1;
    active2full[1] = 2;
    active2full[2] = 3;

    std::vector<double> vertex_data0 = {x0, y0, x0 + dx, y0, x0, y0 + dx};

    std::vector<int> neighbor_count0 = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    std::vector<int> level_values0 = {0, 0, 0};

    std::vector<int> neighbor_values0 = {0, 2, 0, 3, 1, 0, 0, -1, 0, -1, 1, 0};
    std::vector<int> bnd_cells0 = {1};

    std::vector<int> bnd_procs0 = {1};

    vertex_data = std::move(vertex_data0);
    neighbor_count = std::move(neighbor_count0);
    level_values = std::move(level_values0);
    neighbor_values = std::move(neighbor_values0);
    bnd_cells = std::move(bnd_cells0);
    bnd_procs = std::move(bnd_procs0);

  } else if (rtt_c4::node() == 1) {
    full2active.resize(1);
    full2active[0] = 1;
    active2full.resize(1);
    active2full[0] = 1;

    std::vector<double> vertex_data0 = {x0 + dx, y0 + dx};
    std::vector<int> neighbor_count0 = {1, 1, 1, 1};
    std::vector<int> level_values0 = {0};
    std::vector<int> neighbor_values0 = {-2, 0, -1, 0};
    std::vector<int> bnd_cells0 = {2, 3};
    std::vector<int> bnd_procs0 = {0, 0};

    vertex_data = std::move(vertex_data0);
    neighbor_count = std::move(neighbor_count0);
    level_values = std::move(level_values0);
    neighbor_values = std::move(neighbor_values0);
    bnd_cells = std::move(bnd_cells0);
    bnd_procs = std::move(bnd_procs0);
  }

  md.num_across_face = neighbor_count.data();
  md.xyz_extents = vertex_data.data();
  md.cell_levels = level_values.data();
  md.cell_across_face = neighbor_values.data();
  md.processor_boundary_cells = bnd_cells.data();
  md.processor_boundary_procs = bnd_procs.data();

  md.num_active_cells = static_cast<int>(level_values.size());
  md.num_full_cells = md.num_active_cells;
  ghost_size = bnd_cells.size();

  CAMR_Mesh<2> mesh(4, md.num_active_cells, md.num_full_cells, ghost_size, dx, md.xyz_extents,
                    md.cell_levels, md.num_across_face, md.cell_across_face,
                    md.processor_boundary_cells, md.processor_boundary_procs, full2active.data(),
                    active2full.data());

  // mesh.get_Layout().print(std::cout, 1);

  //std::cout << mesh.get_Layout();

  // Test the equality operator
  FAIL_IF(!(mesh == mesh));
  FAIL_IF(mesh != mesh);
  FAIL_IF_NOT(mesh.full_Mesh());

  std::vector<double> sv = mesh.get_short_vertices(1);

  if (rtt_c4::node() == 0) {
    FAIL_IF_NOT(sv.size() == 3);
    FAIL_IF_NOT(mesh.get_bnd_cells()[0] == 1);
    FAIL_IF_NOT(mesh.get_bnd_procs()[0] == 1);
    FAIL_IF_NOT(mesh.get_active2full().size() == 3);
    FAIL_IF_NOT(mesh.get_full2active().size() == 3);
    FAIL_IF_NOT(mesh.get_vertices().size() == 8);
    FAIL_IF_NOT(mesh.get_levels()[2] == 0);
    FAIL_IF_NOT(mesh.get_num_fe() == 4);
  } else {
    FAIL_IF_NOT(mesh.get_bnd_cells()[1] == 3);
    FAIL_IF_NOT(mesh.get_bnd_procs()[1] == 0);
    FAIL_IF_NOT(mesh.get_active2full().size() == 1);
    FAIL_IF_NOT(mesh.get_full2active().size() == 1);
    FAIL_IF_NOT(mesh.get_vertices().size() == 6);
    FAIL_IF_NOT(mesh.get_levels()[0] == 0);
  }

  if (ut.numFails == 0) {
    if (rtt_c4::node() == 0)
      PASSMSG("rank=0 2D CAMR_host (1-level) test passed");
    else
      PASSMSG("rank=1 2D CAMR_host (1-level) test passed");
  }

  return;
}
//--------------------------------------------------//
void test_2D_amr_host(rtt_c4::ParallelUnitTest &ut) {
  /*
--
--
-- Create a simple 2D Rage_Mesh (test for domian decomposed,
--  non-orthgonal dd boundary
--   * P1
--   * // ------------------
--   * // |       |-2  |   |
--   * // |       |--------|
--   * // | 3     | 4  |-1 |
--   * // ------------------
--   * // |       |        |
--   * // |       |        |
--   * // | 1     | 2      |
--   * // ------------------
--   * ////////////////////
--   *
--   * P0
--   * // ------------------
--   * // |       | 3  | 1 |
--   * // |       |--------|
--   * // | -3    | -1 | 2 |
--   * // ------------------
--   * // |       |        |
--   * // |       |        |
--   * // |       | -2     |
--   * // ------------------
--   * ////////////////////
--   *
-- mesh sizes
--NOTE about cell_across_face indexing
--We assume (in 2D ordering of  mesh.next_cells() to be (lower, upper)
--e.g., see the pic below
--    !  RZ:
--    !
--    !      | 1 | 2 |
--    !      ---------
--    !    2 |       | 2
--    !  R --|       |--
--    !    1 |       | 1
--    !      ---------
--    !      | 1 | 2 |
--    !
--    !          Z
   * Top level dimension is 1.0cm. Lower left corner is at (0,0).
   */

  std::cout << "in dd host" << std::endl;

  CAMR_Mesh_Data md;

  //    size_t full_size=4;
  size_t ghost_size = 0;

  double dx = 0.5;
  double x0 = 0;
  double y0 = 0;
  double x1 = dx;
  double y1 = dx;
  double dx1 = dx / 2.;
  double dy1 = dx1;

  std::vector<int> full2active;
  std::vector<int> active2full;

  std::vector<double> vertex_data;
  std::vector<int> neighbor_count;
  std::vector<int> level_values;
  std::vector<int> neighbor_values;
  std::vector<int> bnd_cells;
  std::vector<int> bnd_procs;

  if (rtt_c4::node() == 0) {
    full2active.resize(3);
    full2active[0] = 1;
    full2active[1] = 2;
    full2active[2] = 3;
    active2full.resize(3);
    active2full[0] = 1;
    active2full[1] = 2;
    active2full[2] = 3;

    std::vector<double> vertex_data0 = {x1 + dx1, y1 + dy1, x1 + dx1, y1, x1, y1 + dy1};

    std::vector<int> neighbor_count0 = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    std::vector<int> level_values0 = {1, 1, 1};

    std::vector<int> neighbor_values0 = {3, 0, 2, 0, -1, 0, -2, 1, -3, 1, -1, 0};
    std::vector<int> bnd_cells0 = {4, 2, 3};

    std::vector<int> bnd_procs0 = {1, 1, 1};

    vertex_data = std::move(vertex_data0);
    neighbor_count = std::move(neighbor_count0);
    level_values = std::move(level_values0);
    neighbor_values = std::move(neighbor_values0);
    bnd_cells = std::move(bnd_cells0);
    bnd_procs = std::move(bnd_procs0);

  } else if (rtt_c4::node() == 1) {
    full2active.resize(4);
    full2active[0] = 1;
    full2active[1] = 2;
    full2active[2] = 3;
    full2active[3] = 4;
    active2full.resize(4);
    active2full[0] = 1;
    active2full[1] = 2;
    active2full[2] = 3;
    active2full[3] = 4;

    std::vector<double> vertex_data0 = {x0, y0, x0 + dx, y0, x0, y0 + dx, x1, y1};

    std::vector<int> neighbor_count0 = {1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1};
    std::vector<int> level_values0 = {0, 0, 0, 1};
    std::vector<int> neighbor_values0 = {0, 2, 0, 3, 1, 0, 0, 4, -1, 0, 4, -2, 1, 0, 3, -1, 2, -2};
    std::vector<int> bnd_cells0 = {2, 3};
    std::vector<int> bnd_procs0 = {0, 0};

    vertex_data = std::move(vertex_data0);
    neighbor_count = std::move(neighbor_count0);
    level_values = std::move(level_values0);
    neighbor_values = std::move(neighbor_values0);
    bnd_cells = std::move(bnd_cells0);
    bnd_procs = std::move(bnd_procs0);
  }

  md.num_across_face = neighbor_count.data();
  md.xyz_extents = vertex_data.data();
  md.cell_levels = level_values.data();
  md.cell_across_face = neighbor_values.data();
  md.processor_boundary_cells = bnd_cells.data();
  md.processor_boundary_procs = bnd_procs.data();

  md.num_active_cells = static_cast<int>(level_values.size());
  md.num_full_cells = md.num_active_cells;
  ghost_size = bnd_cells.size();

  CAMR_Mesh<2> mesh(7, md.num_active_cells, md.num_full_cells, ghost_size, dx, md.xyz_extents,
                    md.cell_levels, md.num_across_face, md.cell_across_face,
                    md.processor_boundary_cells, md.processor_boundary_procs, full2active.data(),
                    active2full.data());

  // mesh.get_Layout().print(std::cout, 1);
  if (rtt_c4::node() == 0)
    std::cout << mesh.get_Layout();

  // Test the equality operator
  FAIL_IF(!(mesh == mesh));
  FAIL_IF(mesh != mesh);
  FAIL_IF_NOT(mesh.full_Mesh());

  std::vector<double> sv = mesh.get_short_vertices(1);

  if (rtt_c4::node() == 0) {
    FAIL_IF_NOT(sv.size() == 3);
    FAIL_IF_NOT(mesh.get_bnd_cells()[0] == 4);
    FAIL_IF_NOT(mesh.get_bnd_procs()[0] == 1);
    FAIL_IF_NOT(mesh.get_active2full().size() == 3);
    FAIL_IF_NOT(mesh.get_full2active().size() == 3);
    FAIL_IF_NOT(mesh.get_vertices().size() == 12);
    FAIL_IF_NOT(mesh.get_levels()[2] == 1);
    FAIL_IF_NOT(mesh.get_num_fe() == 4);

    FAIL_IF_NOT(mesh.is_coarse(2, 3) == 2);
  } else {
    FAIL_IF_NOT(mesh.get_bnd_cells()[1] == 3);
    FAIL_IF_NOT(mesh.get_bnd_procs()[1] == 0);
    FAIL_IF_NOT(mesh.get_active2full().size() == 4);
    FAIL_IF_NOT(mesh.get_full2active().size() == 4);
    FAIL_IF_NOT(mesh.get_vertices().size() == 12);
    FAIL_IF_NOT(mesh.get_levels()[0] == 0);
  }

  if (ut.numFails == 0) {
    if (rtt_c4::node() == 0)
      PASSMSG("rank=0 2D CAMR_host (AMR) test passed");
    else
      PASSMSG("rank=1 2D CAMR_host (AMR) test passed");
  }

  return;
}
//------------------------------------------------------------------------------------------------//
//--------------------------------------------------//
void test_3D_amr_host(rtt_c4::ParallelUnitTest &ut) {
  /*
--
--
-- Create a simple 2D Rage_Mesh (test for domian decomposed,
--  
--  non-orthgonal dd boundary
--   * P0 y @z=0
--        ^
--   * // ------------------
--   * // |       |-2  |   |
--   * // |       |--------|
--   * // | 1     | -1 |   |
--   * // ------------------  x
--   * ////////////////////
--   *
--   * P1 y
--   * // ------------------
--   * // |       | 3 |  4 |
--   * // |       |--------|
--   * // | -1    | 1 | 2  |
--   * // ------------------  x
--   * ////////////////////
--   *
--   * P0 y @z= 0.5
--        ^
--   * // ------------------
--   * // |       |  2 | -4|
--   * // |       |--------|
--   * // | 1     | -3 |   |
--   * // ------------------  x
--   * ////////////////////
--   *
--   * P1 y
--   * // ------------------
--   * // |       |-2 |  7 |
--   * // |       |--------|
--   * // | -1    | 5 | 6  |
--   * // ------------------  x
--   * ////////////////////
--   *
-- mesh sizes
   *
   * Top level dimension is 1.0cm. Lower left corner is at (0,0).
   */

  CAMR_Mesh_Data md;

  //    size_t full_size=4;
  size_t ghost_size = 0;

  double dx = 1.;
  double x0 = 0;
  double y0 = 0;
  double z0 = 0;

  double x1 = dx;
  double dx1 = dx / 2.;
  double dy1 = dx1;
  double dz1 = dx1;

  std::vector<int> full2active;
  std::vector<int> active2full;

  std::vector<double> vertex_data;
  std::vector<int> neighbor_count;
  std::vector<int> level_values;
  std::vector<int> neighbor_values;
  std::vector<int> bnd_cells;
  std::vector<int> bnd_procs;

  if (rtt_c4::node() == 0) {
    full2active.resize(2);
    full2active[0] = 1;
    full2active[1] = 2;
    active2full.resize(2);
    active2full[0] = 1;
    active2full[1] = 2;

    std::vector<double> vertex_data0 = {x0, y0, z0, x1, y0 + dy1, z0 + dz1};

    std::vector<int> neighbor_count0 = {1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    std::vector<int> level_values0 = {0, 1};

    std::vector<int> neighbor_values0 = {0, -1, -2, -3, 2, 0, 0, 0, 0, 1, -4, -3, 0, -2, 0};
    std::vector<int> bnd_cells0 = {1, 3, 5, 7};

    std::vector<int> bnd_procs0 = {1, 1, 1, 1};

    vertex_data = std::move(vertex_data0);
    neighbor_count = std::move(neighbor_count0);
    level_values = std::move(level_values0);
    neighbor_values = std::move(neighbor_values0);
    bnd_cells = std::move(bnd_cells0);
    bnd_procs = std::move(bnd_procs0);

  } else if (rtt_c4::node() == 1) {
    full2active.resize(7);
    full2active[0] = 1;
    full2active[1] = 2;
    full2active[2] = 3;
    full2active[3] = 4;
    full2active[4] = 5;
    full2active[5] = 6;
    full2active[6] = 7;
    active2full.resize(7);
    active2full[0] = 1;
    active2full[1] = 2;
    active2full[2] = 3;
    active2full[3] = 4;
    active2full[4] = 5;
    active2full[5] = 6;
    active2full[6] = 7;

    std::vector<double> vertex_data0 = {
        x1, y0, z0, x1 + dx1, y0,       z0, x1,       y0 + dx1, z0,       x1 + dx1, y0 + dx1,
        z0, x1, y0, z0 + dx1, x1 + dx1, y0, z0 + dx1, x1 + dx1, y0 + dx1, z0 + dx1};

    std::vector<int> neighbor_count0 = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    std::vector<int> level_values0 = {1, 1, 1, 1, 1, 1, 1};

    std::vector<int> neighbor_values0 = {-1, 2, 0, 3,  0, 5, 1, 0, 0,  4, 0,  6, -1, 4,
                                         1,  0, 0, -2, 3, 0, 2, 0, 0,  7, -1, 6, 0,  -2,
                                         1,  0, 5, 0,  0, 7, 2, 0, -2, 0, 6,  0, 4,  0};
    std::vector<int> bnd_cells0 = {1, 2};
    std::vector<int> bnd_procs0 = {0, 0};

    vertex_data = std::move(vertex_data0);
    neighbor_count = std::move(neighbor_count0);
    level_values = std::move(level_values0);
    neighbor_values = std::move(neighbor_values0);
    bnd_cells = std::move(bnd_cells0);
    bnd_procs = std::move(bnd_procs0);
  }

  md.num_across_face = neighbor_count.data();
  md.xyz_extents = vertex_data.data();
  md.cell_levels = level_values.data();
  md.cell_across_face = neighbor_values.data();
  md.processor_boundary_cells = bnd_cells.data();
  md.processor_boundary_procs = bnd_procs.data();

  md.num_active_cells = static_cast<int>(level_values.size());
  md.num_full_cells = md.num_active_cells;
  ghost_size = bnd_cells.size();

  CAMR_Mesh<3> mesh(9, md.num_active_cells, md.num_full_cells, ghost_size, dx, md.xyz_extents,
                    md.cell_levels, md.num_across_face, md.cell_across_face,
                    md.processor_boundary_cells, md.processor_boundary_procs, full2active.data(),
                    active2full.data());

  // mesh.get_Layout().print(std::cout, 1);

  std::cout << mesh.get_Layout();

  // Test the equality operator
  FAIL_IF(!(mesh == mesh));
  FAIL_IF(mesh != mesh);
  FAIL_IF_NOT(mesh.full_Mesh());

  std::vector<double> sv = mesh.get_short_vertices(1);

  if (rtt_c4::node() == 0) {
    FAIL_IF_NOT(sv.size() == 4);
    FAIL_IF_NOT(mesh.get_bnd_cells()[0] == 1);
    FAIL_IF_NOT(mesh.get_bnd_procs()[0] == 1);
    FAIL_IF_NOT(mesh.get_active2full().size() == 2);
    FAIL_IF_NOT(mesh.get_full2active().size() == 2);
    FAIL_IF_NOT(mesh.get_vertices().size() == 18);
    FAIL_IF_NOT(mesh.get_levels()[2] == 1);
    FAIL_IF_NOT(mesh.get_levels()[3] == 1);
    FAIL_IF_NOT(mesh.get_num_fe() == 8);

    FAIL_IF_NOT(mesh.is_coarse(2, 1) == 4);

  } else {
    FAIL_IF_NOT(mesh.get_bnd_cells()[1] == 2);
    FAIL_IF_NOT(mesh.get_bnd_procs()[1] == 0);
    FAIL_IF_NOT(mesh.get_active2full().size() == 7);
    FAIL_IF_NOT(mesh.get_full2active().size() == 7);
    FAIL_IF_NOT(mesh.get_vertices().size() == 27);
    FAIL_IF_NOT(mesh.get_levels()[0] == 1);
    FAIL_IF_NOT(mesh.get_levels()[2] == 1);
    FAIL_IF_NOT(mesh.get_levels()[7] == 0);
    FAIL_IF_NOT(mesh.get_levels()[8] == 1);

    FAIL_IF_NOT(mesh.is_coarse(1, 1) == 1);
    FAIL_IF_NOT(mesh.is_coarse(3, 1) == 2);
    FAIL_IF_NOT(mesh.is_coarse(5, 1) == 3);
  }

  if (ut.numFails == 0)
    PASSMSG("3D CAMR_host test passed");
  return;
}

void test_3D_amr_host2(rtt_c4::ParallelUnitTest &ut) {
  /*
--
--
-- Create a simple 2D Rage_Mesh (test for domian decomposed,
--  
--  non-orthgonal dd boundary
--   * P0 y @z=0
--        ^
--   * // ------------------
--   * // |       | 2  |-2 |
--   * // |       |--------|
--   * // | 1     | -1 |   |
--   * // ------------------  x
--   * ////////////////////
--   *
--   * P1 y
--   * // ------------------
--   * // |       | -2 | 3 |
--   * // |       |--------|
--   * // | -1    | 1 | 2  |
--   * // ------------------  x
--   * ////////////////////
--   *
--   * P0 y @z= 0.5
--        ^
--   * // ------------------
--   * // |       |  -4|   |
--   * // |       |--------|
--   * // | 1     | -3 |   |
--   * // ------------------  x
--   * ////////////////////
--   *
--   * P1 y
--   * // ------------------
--   * // |       | 6 |  7 |
--   * // |       |--------|
--   * // | -1    | 4 | 5  |
--   * // ------------------  x
--   * ////////////////////
--   *
-- mesh sizes
   *
   * Top level dimension is 1.0cm. Lower left corner is at (0,0).
   */

  CAMR_Mesh_Data md;

  //    size_t full_size=4;
  size_t ghost_size = 0;

  double dx = 1.;
  double x0 = 0;
  double y0 = 0;
  double z0 = 0;

  double x1 = dx;
  double dx1 = dx / 2.;
  double dy1 = dx1;

  std::vector<int> full2active;
  std::vector<int> active2full;

  std::vector<double> vertex_data;
  std::vector<int> neighbor_count;
  std::vector<int> level_values;
  std::vector<int> neighbor_values;
  std::vector<int> bnd_cells;
  std::vector<int> bnd_procs;

  if (rtt_c4::node() == 0) {
    full2active.resize(2);
    full2active[0] = 1;
    full2active[1] = 2;
    active2full.resize(2);
    active2full[0] = 1;
    active2full[1] = 2;

    std::vector<double> vertex_data0 = {x0, y0, z0, x1, y0 + dy1, z0};

    std::vector<int> neighbor_count0 = {1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    std::vector<int> level_values0 = {0, 1};

    std::vector<int> neighbor_values0 = {0, -1, 2, -3, -4, 0, 0, 0, 0, 1, -2, -1, 0, 0, -4};
    std::vector<int> bnd_cells0 = {1, 3, 4, 6};

    std::vector<int> bnd_procs0 = {1, 1, 1, 1};

    vertex_data = std::move(vertex_data0);
    neighbor_count = std::move(neighbor_count0);
    level_values = std::move(level_values0);
    neighbor_values = std::move(neighbor_values0);
    bnd_cells = std::move(bnd_cells0);
    bnd_procs = std::move(bnd_procs0);

  } else if (rtt_c4::node() == 1) {
    full2active.resize(7);
    full2active[0] = 1;
    full2active[1] = 2;
    full2active[2] = 3;
    full2active[3] = 4;
    full2active[4] = 5;
    full2active[5] = 6;
    full2active[6] = 7;
    active2full.resize(7);
    active2full[0] = 1;
    active2full[1] = 2;
    active2full[2] = 3;
    active2full[3] = 4;
    active2full[4] = 5;
    active2full[5] = 6;
    active2full[6] = 7;

    std::vector<double> vertex_data0 = {
        x1,       y0,       z0, x1 + dx1, y0, z0,       x1 + dx1, y0 + dx1, z0,       x1,      y0,
        z0 + dx1, x1 + dx1, y0, z0 + dx1, x1, y0 + dx1, z0 + dx1, x1 + dx1, y0 + dx1, z0 + dx1};

    std::vector<int> neighbor_count0 = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    std::vector<int> level_values0 = {1, 1, 1, 1, 1, 1, 1};

    std::vector<int> neighbor_values0 = {-1, 2, 0,  -2, 0,  4, 1,  0, 0, 3, 0, 5, -2, 0,
                                         2,  0, 0,  7,  -1, 5, 0,  6, 1, 0, 4, 0, 0,  7,
                                         2,  0, -1, 7,  4,  0, -2, 0, 6, 0, 5, 0, 3,  0};
    std::vector<int> bnd_cells0 = {1, 2};
    std::vector<int> bnd_procs0 = {0, 0};

    vertex_data = std::move(vertex_data0);
    neighbor_count = std::move(neighbor_count0);
    level_values = std::move(level_values0);
    neighbor_values = std::move(neighbor_values0);
    bnd_cells = std::move(bnd_cells0);
    bnd_procs = std::move(bnd_procs0);
  }

  md.num_across_face = neighbor_count.data();
  md.xyz_extents = vertex_data.data();
  md.cell_levels = level_values.data();
  md.cell_across_face = neighbor_values.data();
  md.processor_boundary_cells = bnd_cells.data();
  md.processor_boundary_procs = bnd_procs.data();

  md.num_active_cells = static_cast<int>(level_values.size());
  md.num_full_cells = md.num_active_cells;
  ghost_size = bnd_cells.size();

  CAMR_Mesh<3> mesh(9, md.num_active_cells, md.num_full_cells, ghost_size, dx, md.xyz_extents,
                    md.cell_levels, md.num_across_face, md.cell_across_face,
                    md.processor_boundary_cells, md.processor_boundary_procs, full2active.data(),
                    active2full.data());

  // mesh.get_Layout().print(std::cout, 1);

  std::cout << mesh.get_Layout();

  // Test the equality operator
  FAIL_IF(!(mesh == mesh));
  FAIL_IF(mesh != mesh);
  FAIL_IF_NOT(mesh.full_Mesh());

  std::vector<double> sv = mesh.get_short_vertices(1);

  if (rtt_c4::node() == 0) {
    FAIL_IF_NOT(sv.size() == 4);
    FAIL_IF_NOT(mesh.get_bnd_cells()[1] == 3);
    FAIL_IF_NOT(mesh.get_bnd_procs()[2] == 1);
    FAIL_IF_NOT(mesh.get_active2full().size() == 2);
    FAIL_IF_NOT(mesh.get_full2active().size() == 2);
    FAIL_IF_NOT(mesh.get_vertices().size() == 18);
    FAIL_IF_NOT(mesh.get_levels()[2] == 1);
    FAIL_IF_NOT(mesh.get_levels()[3] == 1);
    FAIL_IF_NOT(mesh.get_num_fe() == 8);

    FAIL_IF_NOT(mesh.is_coarse(2, 1) == 2);

  } else {
    FAIL_IF_NOT(mesh.get_bnd_cells()[1] == 2);
    FAIL_IF_NOT(mesh.get_bnd_procs()[1] == 0);
    FAIL_IF_NOT(mesh.get_active2full().size() == 7);
    FAIL_IF_NOT(mesh.get_full2active().size() == 7);
    FAIL_IF_NOT(mesh.get_vertices().size() == 27);
    FAIL_IF_NOT(mesh.get_levels()[0] == 1);
    FAIL_IF_NOT(mesh.get_levels()[2] == 1);
    FAIL_IF_NOT(mesh.get_levels()[7] == 0);
    FAIL_IF_NOT(mesh.get_levels()[8] == 1);

    FAIL_IF_NOT(mesh.is_coarse(1, 1) == 1);
    FAIL_IF_NOT(mesh.is_coarse(4, 1) == 3);
    FAIL_IF_NOT(mesh.is_coarse(6, 1) == 4);
  }

  if (ut.numFails == 0) {
    if (rtt_c4::node() == 0)
      PASSMSG("rank=0 3D CAMR_host2 test passed");
    else if (rtt_c4::node() == 1)
      PASSMSG("rank=1 3D CAMR_host2 test passed");
  }

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    if (rtt_c4::nodes() == 2) {
      test_2D_host(ut);
      test_2D_amr_host(ut);
      test_3D_amr_host(ut);
      test_3D_amr_host2(ut);
    } else {
      Insist(false, "This test only uses 2  PE.");
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstRage_Mesh.cc
//------------------------------------------------------------------------------------------------//
