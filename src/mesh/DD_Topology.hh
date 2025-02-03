//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/DD_Topology.hh
 * \author H. Park
 * \date   Aug. 1 2023
 * \brief  DD_Topology class definition.
 * \note   Copyright (C) 2023-2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef mesh_DD_Topology_hh
#define mesh_DD_Topology_hh

#include "Topology.hh"

namespace rtt_mesh {

//================================================================================================//
/*!
 * \class DD_Topology
 * \brief Specialization of the topology class for Domain Decomposition (DD).
 *
 * The DD_Topology requires processor-ordered global cell indices, e.g:
 *
 * \verbatim
 * | 1 2 3 4 | 5 6 7 | 8 9 10 | 11 12 |  global cell index space
 * | 1 2 3 4 | 1 2 3 | 1 2  3 |  1  2 |  local cell index space
 * |    0    |   1   |    2   |   3   |  processor id space
 * \endverbatim
 *
 * Thus, to convert from local-to-global cell ids on any processor all that is required are the
 * global cell offsets.  In the example above, we define the following global_offset array:
 *
 * \code
 * global_offset[0] = 0;
 * global_offset[1] = 4;
 * global_offset[2] = 7;
 * global_offset[3] = 10;
 * global_offset[4] = 12;
 * \endcode
 *
 * The range of the global_offset array is [0,Np] where Np is the number of processors.  Note that
 * global_offset[Np] = N_global.  Thus, on any pe, we can define the following relationship between
 * local and global cell id:
 *
 * \code
 * global_cell = local_cell + global_offset[p]
 * \endcode
 *
 * where global_cell has range [1,N_global], local_cell has range [1,N_local[p]], and p is the
 * processor id.  The N_local[p] is easily calculated from global_offset:
 *
 * \code
 * N_local[p] = global_offset[p+1] - global_offset[p]
 * \endcode
 *
 * The number of processors must be greater than 1 to use this class.  For single processor (domain
 * replication) jobs use Rep_Topology.
 *
 * \sa rtt_mesh::Topology
 */
//================================================================================================//

class DD_Topology : public Topology {
public:
  // Typedefs.
  using Base = Topology;
  using packed_type = Base::packed_type;
  using Vec_Int = std::vector<int>;
  using Vec_Int64 = std::vector<int64_t>;

private:
  // >>> DATA

  //! Global cell offset for each processor.
  Vec_Int64 global_offset_d;

  //! Boundary cells.
  Vec_Int ghost_cells_d;

  //! Processors that boundary cells live on.
  Vec_Int ghost_procs_d;

  //! Number of local cells.
  size_t num_local_d;

public:
  // Constructor.
  DD_Topology(const Vec_Int &local_cell_num, Vec_Int ghost_cells, Vec_Int ghost_procs);

  // Unpacking constructor.
  explicit DD_Topology(const packed_type &pt);

  // >>> DERIVED INTERFACE

  //! Number of cells on this processor.
  [[nodiscard]] int64_t num_local_cells() const override { return num_local_d; }

  [[nodiscard]] bool valid_local_cell(int cell) const override {
    return (cell >= 1) && (cell <= static_cast<int>(num_local_d));
  }

  // Get the global cell index.
  [[nodiscard]] inline int64_t global_cell(size_t local_cell) const override;

  // Get the local cell index.
  [[nodiscard]] inline int local_cell(int64_t global_cell) const override;

  //! Get the number of boundary cells on a given processor.
  [[nodiscard]] size_t num_ghost_cells() const override { return ghost_cells_d.size(); }

  // Get the remote local cell index for a boundary cell.
  [[nodiscard]] inline int remote_local_cell(size_t ghost_cell) const override;

  // Get the remote processor id for a boundary cell.
  [[nodiscard]] inline int remote_proc(size_t ghost_cell) const override;

  //! Get the global offset at a rank
  [[nodiscard]] int64_t get_global_offset(int rank) const override { return global_offset_d[rank]; }

  // Packing function.
  [[nodiscard]] packed_type pack() const override;

  // Diagnostic printing.
  void print(std::ostream &out) const override;
};

} // end namespace rtt_mesh

//------------------------------------------------------------------------------------------------//
// INLINE FUNCTION DEFINITIONS
//------------------------------------------------------------------------------------------------//

#include "DD_Topology.i.hh"

#endif // mesh_DD_Topology_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/DD_Topology.hh
//------------------------------------------------------------------------------------------------//
