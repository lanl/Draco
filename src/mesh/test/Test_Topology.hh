//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/test/Test_Topology.hh
 * \author Mathew Cleveland
 * \date   Nov. 2023
 * \brief  Build a DD Topology for testing.
 * \note   Copyright (C) 2023-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_test_Test_Topology_hh
#define rtt_test_Test_Topology_hh

#include "mesh/Topology.hh"
#include <memory>

// Forward declaraddtion.
namespace rtt_mesh {
class Topology;
} // namespace rtt_mesh

namespace rtt_mesh_test {

//================================================================================================//
// BUILD DUMMY INTERFACES FOR TOPOLOGY BUILDER TESTS
//================================================================================================//
// One processor interface; 6 cell mesh replicated.

class One_Proc {
public:
  One_Proc() = default;
  [[nodiscard]] int decomposition_scheme() const;
  [[nodiscard]] int get_num_global_cells() const { return 6; }
  [[nodiscard]] int get_num_bnd_cells() const { return 0; }
  [[nodiscard]] const int *get_bnd_cells() const { return nullptr; }
  [[nodiscard]] const int *get_bnd_procs() const { return nullptr; }
  [[nodiscard]] const int *get_full_active() const { return nullptr; }
  [[nodiscard]] const int *get_active_full() const { return nullptr; }
  [[nodiscard]] int get_num_cells() const { return 0; }
  [[nodiscard]] int get_num_full_cells() const { return 0; }
  [[nodiscard]] bool has_ghost_cells() const { return get_num_full_cells() != 0; }
};

//------------------------------------------------------------------------------------------------//
// Two processor interface; see autodoc/2_proc_DD.fig

class Two_Proc {
private:
  std::vector<int> full_active;
  std::vector<int64_t> active_full;
  std::vector<int> active_full_int;
  std::vector<int> bnd_cells;
  std::vector<int> bnd_procs;

  // Mesh data
  std::vector<int> num_cells_across_face;
  std::vector<int64_t> cells_across_face;
  std::vector<double> r;
  std::vector<double> z;

public:
  Two_Proc();
  [[nodiscard]] int decomposition_scheme() const;
  [[nodiscard]] int get_num_global_cells() const { return 15; }
  [[nodiscard]] int get_num_bnd_cells() const;
  [[nodiscard]] const int *get_bnd_cells() const { return &bnd_cells[0]; }
  [[nodiscard]] const int *get_bnd_procs() const { return &bnd_procs[0]; }
  [[nodiscard]] const int *get_full_active() const { return &full_active[0]; }
  [[nodiscard]] const int64_t *get_active_full() const { return &active_full[0]; }
  [[nodiscard]] const int *get_active_full_int() const { return &active_full_int[0]; }
  [[nodiscard]] int get_num_cells() const;
  [[nodiscard]] int get_num_full_cells() const;
  [[nodiscard]] bool has_ghost_cells() const { return get_num_full_cells() != 0; }

  [[nodiscard]] const int *get_num_cells_across_face() const;
  [[nodiscard]] const int64_t *get_cells_across_face() const;
  [[nodiscard]] const double *get_r_dimensions() const;
  [[nodiscard]] const double *get_z_dimensions() const;
};

//------------------------------------------------------------------------------------------------//
// Two processor interface no ghost cells; see autodoc/2_proc_DD.fig

class Two_Proc_NGC {
private:
  std::vector<int> bnd_cells;
  std::vector<int> bnd_procs;

public:
  Two_Proc_NGC();
  [[nodiscard]] int decomposition_scheme() const;
  [[nodiscard]] int get_num_global_cells() const { return 15; }
  [[nodiscard]] int get_num_bnd_cells() const;
  [[nodiscard]] const int *get_bnd_cells() const { return &bnd_cells[0]; }
  [[nodiscard]] const int *get_bnd_procs() const { return &bnd_procs[0]; }
  [[nodiscard]] const int *get_full_active() const { return nullptr; }
  [[nodiscard]] int get_num_cells() const;
  [[nodiscard]] int get_num_full_cells() const { return 0; }
  [[nodiscard]] bool has_ghost_cells() const { return get_num_full_cells() != 0; }
};

//------------------------------------------------------------------------------------------------//
// Four processor interface; see autodoc/4_proc_DD.fig

class Four_Proc {
private:
  std::vector<int> full_active;
  std::vector<int64_t> active_full;
  std::vector<int> active_full_int;
  std::vector<int> bnd_cells;
  std::vector<int> bnd_procs;

  // Mesh data
  std::vector<int> num_cells_across_face;
  std::vector<int64_t> cells_across_face;
  std::vector<double> r;
  std::vector<double> z;

public:
  Four_Proc();

  // Topology builder interface.
  [[nodiscard]] int decomposition_scheme() const;
  [[nodiscard]] int get_num_global_cells() const { return 18; }
  [[nodiscard]] int get_num_bnd_cells() const;
  [[nodiscard]] const int *get_bnd_cells() const { return &bnd_cells[0]; }
  [[nodiscard]] const int *get_bnd_procs() const { return &bnd_procs[0]; }
  [[nodiscard]] const int *get_full_active() const { return &full_active[0]; }
  [[nodiscard]] const int64_t *get_active_full() const { return &active_full[0]; }
  [[nodiscard]] const int *get_active_full_int() const { return &active_full_int[0]; }
  [[nodiscard]] int get_num_full_cells() const;
  [[nodiscard]] bool has_ghost_cells() const { return get_num_full_cells() != 0; }

  [[nodiscard]] int get_num_cells() const;
  [[nodiscard]] const int *get_num_cells_across_face() const;
  [[nodiscard]] const int64_t *get_cells_across_face() const;
  [[nodiscard]] const double *get_r_dimensions() const;
  [[nodiscard]] const double *get_z_dimensions() const;
};

//================================================================================================//
// BUILD DD TOPOLOGIES
//================================================================================================//
// build a General Topology for 9 cell mesh described above
std::shared_ptr<rtt_mesh::Topology> build_DD_Topology();

} // namespace rtt_mesh_test

#endif // rtt_mesh_test_Test_Topology_hh

//------------------------------------------------------------------------------------------------//
// end of test/Test_Topology.hh
//------------------------------------------------------------------------------------------------//
