//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Rep_Topology.hh
 * \author H. Park
 * \date   Aug. 1 2023
 * \brief  Rep_Topology header file.
 * \note   Copyright (C) 2023-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_mesh_Rep_Topology_hh
#define rtt_mesh_Rep_Topology_hh

#include "dsxx/Assert.hh"
#include "mesh/Topology.hh"
#include <limits>

namespace rtt_mesh {

//================================================================================================//
/*!
 * \class Rep_Topology
 *
 * \brief Specialization of the Topology class for Domain Replication (DR).
 *
 * Rep_Topology is a full replication topology class.  Because full replication topologies are
 * trivial, the Rep_Topology class provides the full functionality of rtt_mesh::Topology with
 * minimal performance and memory overhead.  It must be used instead of rtt_mesh::DD_Topology.
 *
 * The Rep_Topology class assumes simple replication.  That is, a local cell and global cell have
 * the same index on each processor and all processors hold all cells. No boundary cells exist.
 *
 * \sa rtt_mesh::Topology
 */
//================================================================================================//

class Rep_Topology : public Topology {
public:
  // Typedefs.
  using Base = Topology;
  using packed_type = Base::packed_type;

public:
  /*!
  * \brief Rep_Topology constructor.
  *
  * The Rep_Topology constructor only requires the number of global cells in the problem to build an
  * instance of Rep_Topology.  All other functions are trivially derived based on the simplicity of
  * a full replication topology.
  *
  * \param[in] num_global_cells number of global cells in the problem
  */
  explicit Rep_Topology(int64_t num_global_cells);

  // Unpacking constructor.
  explicit Rep_Topology(const packed_type &packed);

  // >>> DERIVED INTERFACE

  //! Return the number of cells on a processor.
  [[nodiscard]] int64_t num_local_cells() const override { return Base::num_global_cells_d; }

  [[nodiscard]] bool valid_local_cell(int cell) const override {
    Require(Base::num_global_cells_d < std::numeric_limits<int>::max());
    return (cell >= 1) && (cell <= static_cast<int>(Base::num_global_cells_d));
  }

  //! Get the global cell index on processor.
  [[nodiscard]] int64_t global_cell(size_t lcell) const override {
    validate(lcell);
    return lcell;
  }

  //! Get the local cell index on processor.
  [[nodiscard]] int local_cell(int64_t gcell) const override {
    validate(gcell);
    Check(gcell < INT32_MAX);
    return static_cast<int>(gcell);
  }

  //! Get the number of boundary cells on a processor.
  [[nodiscard]] size_t num_ghost_cells() const override { return 0; }

  //! Get the remote local cell index for a boundary cell.
  [[nodiscard]] int remote_local_cell(size_t /*c*/) const override {
    Insist(0, "Illegal"); // LCOV_EXCL_LINE
    return 0;             // LCOV_EXCL_LINE
  }

  //! Get the remote processor id for a boundary cell.
  [[nodiscard]] int remote_proc(size_t /*c*/) const override {
    Insist(0, "Illegal"); // LCOV_EXCL_LINE
    return 0;             // LCOV_EXCL_LINE
  }

  //! Get the global offset at a rank
  [[nodiscard]] int64_t get_global_offset(int /*rank*/) const override {
    Insist(0, "Illegal"); // LCOV_EXCL_LINE
    return 0;             // LCOV_EXCL_LINE
  }

  // Pack function.
  [[nodiscard]] packed_type pack() const override;

  // Diagnostics.
  void print(std::ostream &out) const override;

private:
  // >>> IMPLEMENTATION

  // Validate a cell index.
  inline void validate(size_t Remember(cell)) const {
    Require(static_cast<int64_t>(cell) <= Base::num_global_cells_d);
  }
};

} // end namespace rtt_mesh

#endif // rtt_mesh_Rep_Topology_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/Rep_Topology.hh
//------------------------------------------------------------------------------------------------//
