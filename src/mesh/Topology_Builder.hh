//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Topology_Builder.hh
 * \author H. Park
 * \brief  Topology_Builder class definition.
 * \note   Copyright (C) 2023-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef mesh_Topology_Builder_hh
#define mesh_Topology_Builder_hh

#include "mesh/DD_Topology.hh"
#include <memory>

namespace rtt_mesh {

//================================================================================================//
/*!
 * \class Topology_Builder
 * \brief Builds a Rep_Topology or DD_Topology instance.
 *
 * The Topology_Builder builds a Topology object.  The types of topologies that can be built are:
 * - rtt_mesh::Rep_Topology
 * - rtt_mesh::DD_Topology
 * .
 * The choice of topology is determined by the interface function:
 * \code
 *     int decomposition_scheme() const;
 * \endcode
 * This function returns an rtt_mesh::Topology::Decomposition enumeration.  Single processor runs
 * require rtt_mesh::Rep_Topology.  The interface has several required functions.  These are
 * documented in Topology_Builder().
 *
 * The Topology_Builder can build topologies on meshes with ghost, or inactive, cells.  When the
 * mesh has inactive cells all of the boundary cell data fields are defined on the full mesh (the
 * mesh with the inactive cells).  A full-to-active mesh map must be provided by the interface.  The
 * builder does the requisite communication to map the boundary fields to the active mesh. The
 * resulting topology classes are \b always defined on the active mesh.  See Topology_Builder() for
 * more details.  See the tests for examples.
 *
 * The rtt_mesh::DD_Topology is processor-ordered. The constructed global-processor ordering is for
 * internal use only. It \b cannot be used to map back to global ids used by the client; however, it
 * will be consistent within API.  The algorithm used to calculate the local-to-global cell ids is
 * based on a processor-dependent offset. The offset is calculated on each processor \e p using:
 *
 * \f[
 * O_{p} = \sum_{i=0}^{p-1}N_{a}(i)
 * \f]
 *
 * Then, the global ids on processor \e p are calculated using:
 *
 * \f[
 *  g_{p}(l) = O_{p} + l
 * \f]
 *
 * where \e l is the local cell index in the range \f$[1,N_{a}]\f$.
 *
 * \sa rtt_mesh::DD_Topology
 */
//================================================================================================//

class Topology_Builder {
public:
  // Typedefs.
  using SP_Topology = std::shared_ptr<rtt_mesh::Topology>;

private:
  // >>> DATA

  // The parallel scheme.
  int decomposition_;

  // Global number of active cells in the problem.
  int64_t num_global_cells_;

  // >>> DD DATA

  // Number of boundary cells on processor.
  uint32_t num_bnd_cells_;

  // Remote, local, full cell indices for boundary cells on this processor.
  const int *bnd_cells_;

  // Remote processor indices for bnd_cells.
  const int *bnd_procs_;

  // Number of full cells on processor.
  uint32_t num_full_cells_;

  // Number of active cells on processor.
  uint32_t num_active_cells_;

  // Full to active map on processor.
  const int *full2active_;

public:
  // Constructor.
  Topology_Builder(int decomposition, int64_t num_global_cells, uint32_t num_bnd_cells,
                   const int *bnd_cells, const int *bnd_procs, uint32_t num_full_cells,
                   uint32_t num_active_cells, const int *full2active)
      : decomposition_(decomposition), num_global_cells_(num_global_cells),
        num_bnd_cells_(num_bnd_cells), bnd_cells_(bnd_cells), bnd_procs_(bnd_procs),
        num_full_cells_(num_full_cells), num_active_cells_(num_active_cells),
        full2active_(full2active) {
    Require(decomposition_ == rtt_mesh::Topology::DOMAIN_REPLICATION ||
            decomposition_ == rtt_mesh::Topology::DOMAIN_DECOMPOSITION);

    // assign pointers to DD-required fields if the requested topology is DD
    if (decomposition == rtt_mesh::Topology::DOMAIN_DECOMPOSITION) {

#ifdef CHECK_ON
      uint32_t total_cells = num_active_cells;
      rtt_c4::global_sum(total_cells);
      Check(std::cmp_equal(total_cells, num_global_cells));
#endif

      Ensure(num_active_cells_);
      Ensure(std::cmp_less_equal(num_active_cells_, num_global_cells_));
      Ensure(bnd_cells_ != nullptr);
      Ensure(bnd_procs_ != nullptr);
      Ensure(num_full_cells_ ? full2active_ != nullptr : true);
    }

    Ensure(num_global_cells_);
  }
  // Constructor.

  // Build the topology.
  SP_Topology build_Topology();

private:
  // >>> IMPLEMENTATION

  // Typedefs.
  using Vec_Int = rtt_mesh::DD_Topology::Vec_Int;

  // Build local number of cells per processor on each processor.
  void build_cell_numbers(Vec_Int &cells_per_pe);

  // Build DD boundary cell data.
  void build_bnd_data(Vec_Int &bcells, Vec_Int &bprocs);

  // Remap boundary cells on a mesh with ghost or inactive cells.
  void remap_bnd_data(Vec_Int &bcells, Vec_Int &bprocs);
};

} // namespace rtt_mesh

#endif // mesh_Topology_Builder_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/Topology_Builder.hh
//------------------------------------------------------------------------------------------------//
