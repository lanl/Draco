//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/test/Test_Topology.cc
 * \author Mathew Cleveland
 * \date   Nov. 2023
 * \brief  Build DD Meshes for testing
 * \note   Copyright (C) 2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Test_Topology.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include "mesh/DD_Topology.hh"
#include <vector>

namespace rtt_mesh_test {
//------------------------------------------------------------------------------------------------//
// One processor interface, arbitrary 6 cell mesh
//------------------------------------------------------------------------------------------------//

int One_Proc::decomposition_scheme() const { return rtt_mesh::Topology::DOMAIN_REPLICATION; }

//------------------------------------------------------------------------------------------------//
// Two processor interface; see autodoc/2_proc_DD.fig
//------------------------------------------------------------------------------------------------//

Two_Proc::Two_Proc()
    : full_active(), active_full(), active_full_int(), bnd_cells(), bnd_procs(),
      num_cells_across_face(), cells_across_face(), r(), z() {

  if (rtt_c4::node() == 0) {
    full_active.resize(7);
    active_full.resize(6);
    active_full_int.resize(6);
    size_t const Nboundary(3);
    bnd_cells.resize(Nboundary);
    bnd_procs.resize(Nboundary);

    // Full to active map
    full_active[0] = -1; // mother cell, does not exist for active mesh
    full_active[1] = 5;
    full_active[2] = 6;
    full_active[3] = 1;
    full_active[4] = 2;
    full_active[5] = 3;
    full_active[6] = 4;

    for (int j = 0; j < 7; j++) {
      if (full_active[j] != -1)
        active_full_int[full_active[j] - 1] = j + 1;
    }

    bnd_cells[0] = 4;
    bnd_cells[1] = 6;
    bnd_cells[2] = 2;

    bnd_procs[0] = 1;
    bnd_procs[1] = 1;
    bnd_procs[2] = 1;

    for (size_t t = 0; t < active_full_int.size(); t++) {
      active_full[t] = static_cast<int64_t>(active_full_int[t]);
    }

    // build mesh (assume 1 face per side and fixup refined faces below)
    num_cells_across_face = std::vector<int>(4 * active_full.size(), 1);

    // These are for the active mesh.
    // No refinement for cell 1 2 3 4
    // refine cell 5 high y and low y face
    num_cells_across_face[4 * 4 + 2] = 2;
    num_cells_across_face[4 * 4 + 3] = 2;

    // refine cell 6 low x face
    num_cells_across_face[5 * 4] = 2;

    cells_across_face.resize(27);

    // cell 1
    cells_across_face[0] = 1; // reflection, use own index.
    cells_across_face[1] = 3;
    cells_across_face[2] = 1; // reflection, use own index.
    cells_across_face[3] = 2;

    // cell 2
    cells_across_face[4] = 2; // reflection, use own index.
    cells_across_face[5] = 4;
    cells_across_face[6] = 1;
    cells_across_face[7] = 5;

    // cell 3
    cells_across_face[8] = 1;
    cells_across_face[9] = 6;
    cells_across_face[10] = 3; // reflection, use own index.
    cells_across_face[11] = 4;

    // cell 4
    cells_across_face[12] = 2;
    cells_across_face[13] = 6;
    cells_across_face[14] = 3;
    cells_across_face[15] = 5;

    // cell 5
    cells_across_face[16] = 5;  // reflection, use own index.
    cells_across_face[17] = -3; // negative value == other processor
    cells_across_face[18] = 2;
    cells_across_face[19] = 4;
    cells_across_face[20] = -1;
    cells_across_face[21] = -2;

    // cell 6
    cells_across_face[22] = 3;
    cells_across_face[23] = 4;
    cells_across_face[14] = 6; // reflection, use own index.
    cells_across_face[25] = 6; // reflection, use own index.
    cells_across_face[26] = -3;

    r.resize(12);
    z.resize(12);

    // cell 1
    r[0] = 0.0;
    r[1] = 0.5;

    // cell 2
    r[2] = 0.0;
    r[3] = 0.5;

    // cell 3
    r[4] = 0.5;
    r[5] = 1.0;

    // cell 4
    r[6] = 0.5;
    r[7] = 1.0;

    // cell 5
    r[8] = 0.0;
    r[9] = 1.0;

    // cell 6
    r[10] = 1.0;
    r[11] = 2.0;

    // cell 1
    z[0] = 0.0;
    z[1] = 0.5;

    // cell 2
    z[2] = 0.5;
    z[3] = 1.0;

    // cell 3
    z[4] = 0.0;
    z[5] = 0.5;

    // cell 4
    z[6] = 0.5;
    z[7] = 1.0;

    // cell 5
    z[8] = 1.0;
    z[9] = 2.0;

    // cell 6
    z[10] = 0.0;
    z[11] = 1.0;
  }

  if (rtt_c4::node() == 1) {
    full_active.resize(11);
    active_full.resize(9);
    active_full_int.resize(9);
    size_t const Nboundary(2);
    bnd_cells.resize(Nboundary);
    bnd_procs.resize(Nboundary);

    // Full to active map
    full_active[0] = -1; // active mesh does not include mother cells.
    full_active[1] = 5;
    full_active[2] = -1; // active mesh does not include mother cells.
    full_active[3] = 1;
    full_active[4] = 2;
    full_active[5] = 3;
    full_active[6] = 4;
    full_active[7] = 6;
    full_active[8] = 7;
    full_active[9] = 8;
    full_active[10] = 9;

    for (int j = 0; j < 11; j++) {
      if (full_active[j] != -1)
        active_full_int[full_active[j] - 1] = j + 1;
    }

    bnd_cells[0] = 3;
    bnd_cells[1] = 2;

    bnd_procs[0] = 0;
    bnd_procs[1] = 0;

    for (size_t t = 0; t < active_full_int.size(); t++) {
      active_full[t] = static_cast<int64_t>(active_full_int[t]);
    }

    // build mesh
    num_cells_across_face = std::vector<int>(4 * active_full.size(), 1);

    // cell 1 2 3 4 are not refined
    // cell 5 uppwer y face is refined
    num_cells_across_face[4 * 4 + 3] = 2;
    // cell 6 7 8 9 are not refined

    cells_across_face.resize(37);

    // cell 1
    cells_across_face[0] = 1; // own index -> reflectiing
    cells_across_face[1] = 3;
    cells_across_face[2] = -2; // negative -> boundary processor.
    cells_across_face[3] = 2;

    // cell 2
    cells_across_face[4] = 2; // own index -> reflectiing
    cells_across_face[5] = 4;
    cells_across_face[6] = 1;
    cells_across_face[7] = 0; // zero -> vacuum

    // cell 3
    cells_across_face[8] = 1;
    cells_across_face[9] = 6;
    cells_across_face[10] = -2; // negative -> boundary processor.
    cells_across_face[11] = 4;

    // cell 4
    cells_across_face[12] = 2;
    cells_across_face[13] = 7;
    cells_across_face[14] = 3;
    cells_across_face[15] = 0; // zero -> vacuum

    // cell 5
    cells_across_face[16] = -2; // negative -> boundary processor.
    cells_across_face[17] = 5;  // own index -> reflectiing
    cells_across_face[18] = -1; // negative -> boundary processor.
    cells_across_face[19] = 6;
    cells_across_face[20] = 8;

    // cell 6
    cells_across_face[21] = 3;
    cells_across_face[22] = 8;
    cells_across_face[23] = 5;
    cells_across_face[24] = 7;

    // cell 7
    cells_across_face[25] = 4;
    cells_across_face[26] = 9;
    cells_across_face[27] = 6;
    cells_across_face[28] = 0; // zero -> vacuum

    // cell 8
    cells_across_face[29] = 6;
    cells_across_face[30] = 8; // own index -> reflectiing
    cells_across_face[31] = 5;
    cells_across_face[32] = 9;

    // cell 9
    cells_across_face[33] = 7;
    cells_across_face[34] = 9; // own index -> reflectiing
    cells_across_face[35] = 8;
    cells_across_face[36] = 0; // zero -> vacuum

    r.resize(18);
    z.resize(18);

    // cell 1
    r[0] = 0.0;
    r[1] = 0.5;

    // cell 2
    r[2] = 0.0;
    r[3] = 0.5;

    // cell 3
    r[4] = 0.5;
    r[5] = 1.0;

    // cell 4
    r[6] = 0.5;
    r[7] = 1.0;

    // cell 5
    r[8] = 1.0;
    r[9] = 2.0;

    // cell 6
    r[10] = 1.0;
    r[11] = 1.5;

    // cell 7
    r[12] = 1.0;
    r[13] = 1.5;

    // cell 8
    r[14] = 1.5;
    r[15] = 2.0;

    // cell 9
    r[16] = 1.5;
    r[17] = 2.0;

    // cell 1
    z[0] = 2.0;
    z[1] = 2.5;

    // cell 2
    z[2] = 2.5;
    z[3] = 3.0;

    // cell 3
    z[4] = 2.0;
    z[5] = 2.5;

    // cell 4
    z[6] = 2.5;
    z[7] = 3.0;

    // cell 5
    z[8] = 1.0;
    z[9] = 2.0;

    // cell 6
    z[10] = 2.0;
    z[11] = 2.5;

    // cell 7
    z[12] = 2.5;
    z[13] = 3.0;

    // cell 8
    z[14] = 2.0;
    z[15] = 2.5;

    // cell 9
    z[16] = 2.5;
    z[17] = 3.0;
  }
}

//------------------------------------------------------------------------------------------------//

int Two_Proc::decomposition_scheme() const { return rtt_mesh::Topology::DOMAIN_DECOMPOSITION; }

//------------------------------------------------------------------------------------------------//

int Two_Proc::get_num_bnd_cells() const {
  if (rtt_c4::node() == 0)
    return 3;
  if (rtt_c4::node() == 1)
    return 2;

  return 0;
}

//------------------------------------------------------------------------------------------------//

int Two_Proc::get_num_cells() const {
  if (rtt_c4::node() == 0)
    return 6;
  if (rtt_c4::node() == 1)
    return 9;

  return 0;
}

//------------------------------------------------------------------------------------------------//

int Two_Proc::get_num_full_cells() const {
  if (rtt_c4::node() == 0)
    return 7;
  if (rtt_c4::node() == 1)
    return 11;

  return 0;
}

//------------------------------------------------------------------------------------------------//

const int *Two_Proc::get_num_cells_across_face() const { return &num_cells_across_face[0]; }

//------------------------------------------------------------------------------------------------//

const int64_t *Two_Proc::get_cells_across_face() const { return &cells_across_face[0]; }

//------------------------------------------------------------------------------------------------//

const double *Two_Proc::get_r_dimensions() const { return &r[0]; }

//------------------------------------------------------------------------------------------------//

const double *Two_Proc::get_z_dimensions() const { return &z[0]; }

//================================================================================================//
// BUILD DD TOPOLOGIES
// ================================================================================================//

// build a DD Topology for 9 cell mesh described above; see autodoc/DD_Mesh.fig

std::shared_ptr<rtt_mesh::Topology> build_DD_Topology() {
  using rtt_mesh::DD_Topology;
  using rtt_mesh::Topology;
  using std::vector;

  Require(rtt_c4::nodes() == 4);

  std::shared_ptr<Topology> topology;

  // num_global = 9;
  vector<int> bc;
  vector<int> bp;

  // cells per pe
  vector<int> cpp(rtt_c4::nodes());
  cpp[0] = 2;
  cpp[1] = 2;
  cpp[2] = 2;
  cpp[3] = 3;

  if (rtt_c4::node() == 0) {
    bc.resize(3);
    bp.resize(3);

    bc[0] = 1;
    bc[1] = 2;
    bc[2] = 1;

    bp[0] = 1;
    bp[1] = 1;
    bp[2] = 2;
  } else if (rtt_c4::node() == 1) {
    bc.resize(5);
    bp.resize(5);

    bc[0] = 1;
    bc[1] = 2;
    bc[2] = 1;
    bc[3] = 2;
    bc[4] = 1;

    bp[0] = 0;
    bp[1] = 0;
    bp[2] = 2;
    bp[3] = 2;
    bp[4] = 3;
  } else if (rtt_c4::node() == 2) {
    bc.resize(5);
    bp.resize(5);

    bc[0] = 2;
    bc[1] = 1;
    bc[2] = 2;
    bc[3] = 2;
    bc[4] = 3;

    bp[0] = 0;
    bp[1] = 1;
    bp[2] = 1;
    bp[3] = 3;
    bp[4] = 3;
  } else if (rtt_c4::node() == 3) {
    bc.resize(3);
    bp.resize(3);

    bc[0] = 2;
    bc[1] = 1;
    bc[2] = 2;

    bp[0] = 1;
    bp[1] = 2;
    bp[2] = 2;
  }

  topology = std::make_shared<rtt_mesh::DD_Topology>(cpp, bc, bp);
  return topology;
}

//------------------------------------------------------------------------------------------------//
// Two processor interface no ghost cells; see autodoc/2_proc_DD.fig
//------------------------------------------------------------------------------------------------//

Two_Proc_NGC::Two_Proc_NGC() : bnd_cells(), bnd_procs() {
  if (rtt_c4::node() == 0) {
    bnd_cells.resize(3);
    bnd_procs.resize(3);

    bnd_cells[0] = 1;
    bnd_cells[1] = 3;
    bnd_cells[2] = 5;

    bnd_procs[0] = 1;
    bnd_procs[1] = 1;
    bnd_procs[2] = 1;
  }

  if (rtt_c4::node() == 1) {
    bnd_cells.resize(2);
    bnd_procs.resize(2);

    bnd_cells[0] = 6;
    bnd_cells[1] = 5;

    bnd_procs[0] = 0;
    bnd_procs[1] = 0;
  }
}

//------------------------------------------------------------------------------------------------//

int Two_Proc_NGC::decomposition_scheme() const { return rtt_mesh::Topology::DOMAIN_DECOMPOSITION; }

//------------------------------------------------------------------------------------------------//

int Two_Proc_NGC::get_num_bnd_cells() const {
  if (rtt_c4::node() == 0)
    return 3;
  if (rtt_c4::node() == 1)
    return 2;

  return 0;
}

//------------------------------------------------------------------------------------------------//

int Two_Proc_NGC::get_num_cells() const {
  if (rtt_c4::node() == 0)
    return 6;
  if (rtt_c4::node() == 1)
    return 9;

  return 0;
}

//------------------------------------------------------------------------------------------------//
// Four processor interface; see autodoc/4_proc_DD.fig
//------------------------------------------------------------------------------------------------//

Four_Proc::Four_Proc()
    : full_active(), active_full(), bnd_cells(), bnd_procs(), num_cells_across_face(),
      cells_across_face(), r(), z() {

  if (rtt_c4::node() == 0) {
    full_active.resize(2);
    active_full.resize(2);
    active_full_int.resize(2);
    bnd_cells.resize(8);
    bnd_procs.resize(8);

    full_active[0] = 1;
    full_active[1] = 2;

    active_full_int[0] = 1;
    active_full_int[1] = 2;

    bnd_cells[0] = 5;
    bnd_cells[1] = 7;
    bnd_cells[2] = 2;
    bnd_cells[3] = 1;
    bnd_cells[4] = 3;
    bnd_cells[5] = 5;
    bnd_cells[6] = 7;
    bnd_cells[7] = 9;

    bnd_procs[0] = 3;
    bnd_procs[1] = 3;
    bnd_procs[2] = 3;
    bnd_procs[3] = 1;
    bnd_procs[4] = 2;
    bnd_procs[5] = 2;
    bnd_procs[6] = 2;
    bnd_procs[7] = 2;

    for (size_t t = 0; t < active_full_int.size(); t++) {
      active_full[t] = static_cast<int64_t>(active_full_int[t]);
    }

    // make the mesh
    num_cells_across_face = std::vector<int>(8, 1);

    // cell 1 refine upper and lower y
    num_cells_across_face[2] = 2;
    num_cells_across_face[3] = 2;

    // cell 2 refine lower y
    num_cells_across_face[1 * 4 + 3] = 2;

    cells_across_face.resize(11);

    // cell 1
    cells_across_face[0] = 1;
    cells_across_face[1] = 2;
    cells_across_face[2] = -1;
    cells_across_face[3] = -2;
    cells_across_face[4] = -5;
    cells_across_face[5] = -6;

    // cell 2
    cells_across_face[6] = 1;
    cells_across_face[7] = -4;
    cells_across_face[8] = -3;
    cells_across_face[9] = -7;
    cells_across_face[10] = -8;

    r.resize(4);
    z.resize(4);

    // cell 1
    r[0] = 0.0;
    r[1] = 1.0;

    // cell 2
    r[2] = 1.0;
    r[3] = 2.0;

    // cell 1
    z[0] = 1.0;
    z[1] = 2.0;

    // cell 2
    z[2] = 1.0;
    z[3] = 2.0;
  }

  if (rtt_c4::node() == 1) {
    full_active.resize(2);
    active_full.resize(2);
    active_full_int.resize(2);
    bnd_cells.resize(4);
    bnd_procs.resize(4);

    full_active[0] = 1;
    full_active[1] = 2;

    active_full_int[0] = 1;
    active_full_int[1] = 2;

    bnd_cells[0] = 3;
    bnd_cells[1] = 2;
    bnd_cells[2] = 9;
    bnd_cells[3] = 10;

    bnd_procs[0] = 3;
    bnd_procs[1] = 0;
    bnd_procs[2] = 2;
    bnd_procs[3] = 2;

    for (size_t t = 0; t < active_full_int.size(); t++) {
      active_full[t] = static_cast<int64_t>(active_full_int[t]);
    }

    // make the mesh
    num_cells_across_face = std::vector<int>(8, 1);

    // cell 1 is not refined
    // cell 2 low x face is refined
    num_cells_across_face[1 * 4] = 2;

    cells_across_face.resize(9);

    // cell 1
    cells_across_face[0] = -2;
    cells_across_face[1] = 1;
    cells_across_face[2] = -1;
    cells_across_face[3] = 2;

    // cell 2
    cells_across_face[4] = -3;
    cells_across_face[5] = -4;
    cells_across_face[6] = 2;
    cells_across_face[7] = 1;
    cells_across_face[8] = 0;

    r.resize(4);
    z.resize(4);

    // cell 1
    r[0] = 2.0;
    r[1] = 3.0;

    // cell 2
    r[2] = 2.0;
    r[3] = 3.0;

    // cell 1
    z[0] = 1.0;
    z[1] = 2.0;

    // cell 2
    z[2] = 2.0;
    z[3] = 3.0;
  }

  if (rtt_c4::node() == 2) {
    full_active.resize(10);
    active_full.resize(8);
    active_full_int.resize(8);
    bnd_cells.resize(3);
    bnd_procs.resize(3);

    full_active[0] = -1;
    full_active[1] = -1;
    full_active[2] = 1;
    full_active[3] = 2;
    full_active[4] = 3;
    full_active[5] = 4;
    full_active[6] = 5;
    full_active[7] = 6;
    full_active[8] = 7;
    full_active[9] = 8;

    active_full_int[0] = 3;
    active_full_int[1] = 4;
    active_full_int[2] = 5;
    active_full_int[3] = 6;
    active_full_int[4] = 7;
    active_full_int[5] = 8;
    active_full_int[6] = 9;
    active_full_int[7] = 10;

    bnd_cells[0] = 1;
    bnd_cells[1] = 2;
    bnd_cells[2] = 2;

    bnd_procs[0] = 0;
    bnd_procs[1] = 0;
    bnd_procs[2] = 1;

    for (size_t t = 0; t < active_full_int.size(); t++) {
      active_full[t] = static_cast<int64_t>(active_full_int[t]);
    }

    // make the mesh
    num_cells_across_face = std::vector<int>(32, 1);
    cells_across_face.resize(32);

    // cell 1
    cells_across_face[0] = 1;
    cells_across_face[1] = 3;
    cells_across_face[2] = -1;
    cells_across_face[3] = 2;

    // cell 2
    cells_across_face[4] = 2;
    cells_across_face[5] = 4;
    cells_across_face[6] = 1;
    cells_across_face[7] = 0;

    // cell 3
    cells_across_face[8] = 1;
    cells_across_face[9] = 5;
    cells_across_face[10] = -1;
    cells_across_face[11] = 4;

    // cell 4
    cells_across_face[12] = 2;
    cells_across_face[13] = 6;
    cells_across_face[14] = 3;
    cells_across_face[15] = 0;

    // cell 5
    cells_across_face[16] = 3;
    cells_across_face[17] = 7;
    cells_across_face[18] = -2;
    cells_across_face[19] = 6;

    // cell 6
    cells_across_face[20] = 4;
    cells_across_face[21] = 8;
    cells_across_face[22] = 5;
    cells_across_face[23] = 0;

    // cell 7
    cells_across_face[24] = 5;
    cells_across_face[25] = -3;
    cells_across_face[26] = -2;
    cells_across_face[27] = 8;

    // cell 8
    cells_across_face[28] = 6;
    cells_across_face[29] = -3;
    cells_across_face[30] = 7;
    cells_across_face[31] = 0;

    r.resize(16);
    z.resize(16);

    // cell 1
    z[0] = 2.0;
    z[1] = 2.5;

    // cell 2
    z[2] = 2.5;
    z[3] = 3.0;

    // cell 3
    z[4] = 2.0;
    z[5] = 2.5;

    // cell 4
    z[6] = 2.5;
    z[7] = 3.0;

    // cell 5
    z[8] = 2.0;
    z[9] = 2.5;

    // cell 6
    z[10] = 2.5;
    z[11] = 3.0;

    // cell 7
    z[12] = 2.0;
    z[13] = 2.5;

    // cell 8
    z[14] = 2.5;
    z[15] = 3.0;

    // cell 1
    r[0] = 0.0;
    r[1] = 0.5;

    // cell 2
    r[2] = 0.0;
    r[3] = 0.5;

    // cell 3
    r[4] = 0.5;
    r[5] = 1.0;

    // cell 4
    r[6] = 0.5;
    r[7] = 1.0;

    // cell 5
    r[8] = 1.0;
    r[9] = 1.5;

    // cell 6
    r[10] = 1.0;
    r[11] = 1.5;

    // cell 7
    r[12] = 1.5;
    r[13] = 2.0;

    // cell 8
    r[14] = 1.5;
    r[15] = 2.0;
  }

  if (rtt_c4::node() == 3) {
    full_active.resize(7);
    active_full.resize(6);
    active_full_int.resize(6);
    bnd_cells.resize(3);
    bnd_procs.resize(3);

    full_active[0] = -1;
    full_active[1] = 5;
    full_active[2] = 6;
    full_active[3] = 1;
    full_active[4] = 2;
    full_active[5] = 3;
    full_active[6] = 4;

    active_full_int[0] = 4;
    active_full_int[1] = 5;
    active_full_int[2] = 6;
    active_full_int[3] = 7;
    active_full_int[4] = 2;
    active_full_int[5] = 3;

    bnd_cells[0] = 1;
    bnd_cells[1] = 2;
    bnd_cells[2] = 1;

    bnd_procs[0] = 0;
    bnd_procs[1] = 0;
    bnd_procs[2] = 1;

    for (size_t t = 0; t < active_full_int.size(); t++) {
      active_full[t] = static_cast<int64_t>(active_full_int[t]);
    }

    // make the mesh
    num_cells_across_face = std::vector<int>(24, 1);

    // cells 1-4 are not refined
    // cell 5 low x boundary is refined
    num_cells_across_face[4 * 4] = 2;
    // cell 6 is not refined

    cells_across_face.resize(25);

    // cell 1
    cells_across_face[0] = 1;
    cells_across_face[1] = 3;
    cells_across_face[2] = 1;
    cells_across_face[3] = 2;

    // cell 2
    cells_across_face[4] = 2;
    cells_across_face[5] = 4;
    cells_across_face[6] = 1;
    cells_across_face[7] = -1;

    // cell 3
    cells_across_face[8] = 1;
    cells_across_face[9] = 5;
    cells_across_face[10] = 3;
    cells_across_face[11] = 4;

    // cell 4
    cells_across_face[12] = 2;
    cells_across_face[13] = 5;
    cells_across_face[14] = 3;
    cells_across_face[15] = -1;

    // cell 5
    cells_across_face[16] = 3;
    cells_across_face[17] = 4;
    cells_across_face[18] = 6;
    cells_across_face[19] = 5;
    cells_across_face[20] = -2;

    // cell 6
    cells_across_face[21] = 5;
    cells_across_face[22] = 6;
    cells_across_face[23] = 6;
    cells_across_face[24] = -3;

    r.resize(12);
    z.resize(12);

    // cell 1
    r[0] = 0.0;
    r[1] = 0.5;

    // cell 2
    r[2] = 0.0;
    r[3] = 0.5;

    // cell 3
    r[4] = 0.5;
    r[5] = 1.0;

    // cell 4
    r[6] = 0.5;
    r[7] = 1.0;

    // cell 5
    r[8] = 1.0;
    r[9] = 2.0;

    // cell 6
    r[10] = 2.0;
    r[11] = 3.0;

    // cell 1
    z[0] = 0.0;
    z[1] = 0.5;

    // cell 2
    z[2] = 0.5;
    z[3] = 1.0;

    // cell 3
    z[4] = 0.0;
    z[5] = 0.5;

    // cell 4
    z[6] = 0.5;
    z[7] = 1.0;

    // cell 5
    z[8] = 0.0;
    z[9] = 1.0;

    // cell 6
    z[10] = 0.0;
    z[11] = 1.0;
  }
}

//------------------------------------------------------------------------------------------------//

int Four_Proc::decomposition_scheme() const { return rtt_mesh::Topology::DOMAIN_DECOMPOSITION; }

//------------------------------------------------------------------------------------------------//

int Four_Proc::get_num_bnd_cells() const {
  if (rtt_c4::node() == 0)
    return 8;
  if (rtt_c4::node() == 1)
    return 4;
  if (rtt_c4::node() == 2)
    return 3;
  if (rtt_c4::node() == 3)
    return 3;

  return 0;
}

//------------------------------------------------------------------------------------------------//

int Four_Proc::get_num_cells() const {
  if (rtt_c4::node() == 0)
    return 2;
  if (rtt_c4::node() == 1)
    return 2;
  if (rtt_c4::node() == 2)
    return 8;
  if (rtt_c4::node() == 3)
    return 6;

  return 0;
}

//------------------------------------------------------------------------------------------------//

int Four_Proc::get_num_full_cells() const {
  if (rtt_c4::node() == 0)
    return 2;
  if (rtt_c4::node() == 1)
    return 2;
  if (rtt_c4::node() == 2)
    return 10;
  if (rtt_c4::node() == 3)
    return 7;

  return 0;
}

//------------------------------------------------------------------------------------------------//

const int *Four_Proc::get_num_cells_across_face() const { return &num_cells_across_face[0]; }

//------------------------------------------------------------------------------------------------//

const int64_t *Four_Proc::get_cells_across_face() const { return &cells_across_face[0]; }

//------------------------------------------------------------------------------------------------//

const double *Four_Proc::get_r_dimensions() const { return &r[0]; }

//------------------------------------------------------------------------------------------------//

const double *Four_Proc::get_z_dimensions() const { return &z[0]; }

} // namespace rtt_mesh_test

//------------------------------------------------------------------------------------------------//
// end of Test_Topology.cc
//------------------------------------------------------------------------------------------------//
