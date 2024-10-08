//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/DD_Topology.i.hh
 * \author H. Park
 * \date   Aug. 1 2023
 * \brief  Member definitions of class DD_Topology
 * \note   Copyright (C) 2023-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef mesh_DD_Topology_i_hh
#define mesh_DD_Topology_i_hh

#include "c4/global.hh"
#include <limits>

namespace rtt_mesh {

//------------------------------------------------------------------------------------------------//
// INLINE FUNCTION DEFINITIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the global cell index given a local cell on-processor.
 *
 * \param local_cell local cell index in range [1, N_local]
 * \return global cell index in range [1, N_global]
 */
int64_t DD_Topology::global_cell(size_t local_cell) const {
  size_t const node = rtt_c4::node();
  Require(local_cell <= num_local_d);

  int64_t const g = local_cell + global_offset_d[node];

  Ensure(g > global_offset_d[node]);
  Ensure(g <= global_offset_d[node + 1]);
  return g;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the local cell index given a global cell on-processor.
 *
 * If the requested global cell does not live on the current processor a value of zero is returned.
 * Otherwise, the local cell index for the particular global cell on the current processor is
 * returned.
 *
 * \param global_cell requested global cell index in range [1, N_global]
 * \return local cell index of global cell on processor or 0 if the global cell does not live on the
 *        current processor
 */
int DD_Topology::local_cell(int64_t global_cell) const {
  size_t const node = rtt_c4::node();
  Require(Base::num_global_cells_d < std::numeric_limits<int>::max());
  Require(global_cell <= Base::num_global_cells_d);

  // convert to local
  auto local_cell = global_cell - global_offset_d[node];

  // check if it is in range and return
  if (local_cell > 0 && local_cell < static_cast<int64_t>(num_local_d + 1)) {
    Check(local_cell < INT32_MAX);
    return static_cast<int>(local_cell);
  }

  // if the global_cell does not live on this node return 0
  return 0;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the remote local cell index for a boundary cell.
 *
 * \param ghost_cell_id boundary cell index in range [1, num_boundary_cells()]
 * \return local cell index on the node returned by remote_node(boundary_cell)
 */
int DD_Topology::remote_local_cell(size_t ghost_cell_id) const {
  Require(ghost_cell_id > 0);
  Require(ghost_cell_id <= ghost_cells_d.size());
  return ghost_cells_d[ghost_cell_id - 1];
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the remote processor index for a boundary cell.
 *
 * \param ghost_cell_id ghost cell index in range [1, num_ghost_cells()]
 * \return processor id that indicates where the boundary cell lives
 */
int DD_Topology::remote_proc(size_t ghost_cell_id) const {
  Require(ghost_cell_id > 0);
  Require(ghost_cell_id <= ghost_cells_d.size());

  // the boundary cell must live on another processor
  Ensure(ghost_procs_d[ghost_cell_id - 1] != rtt_c4::node());
  return ghost_procs_d[ghost_cell_id - 1];
}

} // end namespace rtt_mesh

#endif // mesh_DD_Topology_i_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/DD_Topology.i.hh
//------------------------------------------------------------------------------------------------//
