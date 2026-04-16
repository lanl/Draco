//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/DD_Topology.cc
 * \author H Park
 * \date   Aug. 1 2023
 * \brief  DD_Topology member definitions.
 * \note   Copyright (C) 2023-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "mesh/DD_Topology.hh"
#include "dsxx/Assert.hh"
#include "dsxx/Packing_Utils.hh"
#include "c4/C4_Functions.hh"
#include "mesh/Topology.hh"
#include <cstddef> // size_t
#include <cstdint> // int64_t
#include <ostream>
#include <string>
#include <utility> // move
#include <vector>

namespace rtt_mesh {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//------------------------------------------------------------------------------------------------//

//! Constructor
DD_Topology::DD_Topology(const Vec_Int &local_cell_num, Vec_Int ghost_cells, Vec_Int ghost_procs)
    : Topology(std::string("DD")), global_offset_d(rtt_c4::nodes() + 1, 0),
      ghost_cells_d(std::move(ghost_cells)), ghost_procs_d(std::move(ghost_procs)), num_local_d(0) {
  size_t const nodes = rtt_c4::nodes();
  size_t const node = rtt_c4::node();
  Require(local_cell_num.size() == nodes);
  Require(nodes > 1);
  Require(global_offset_d.size() == nodes + 1);
  Require(ghost_cells_d.size() == ghost_procs_d.size());

  // build the global_offset array
  for (size_t p = 1, s = nodes + 1; p < s; ++p)
    global_offset_d[p] = local_cell_num[p - 1] + global_offset_d[p - 1];

  // set the global number of cells in the base class
  Base::num_global_cells_d = global_offset_d[nodes];

  // set the number of local cells
  num_local_d = static_cast<size_t>(local_cell_num[node]);

  Ensure(global_offset_d[0] == 0);
  Ensure(global_offset_d[nodes] == Base::num_global_cells_d);
  Ensure(Base::get_parallel_scheme() == "DD");
}

//------------------------------------------------------------------------------------------------//
//! Unpacking constructor.
DD_Topology::DD_Topology(const packed_type &packed)
    : Topology(std::string("DD")), global_offset_d(std::vector<int64_t>(rtt_c4::nodes() + 1, 0)),
      ghost_cells_d(std::vector<int>()), ghost_procs_d(std::vector<int>()), num_local_d(0) {
  size_t const nodes = rtt_c4::nodes();
  size_t const node = rtt_c4::node();
  Require(nodes > 1);
  Require(packed.size() >= 4 * sizeof(int));

  // unpack the data
  rtt_dsxx::Unpacker u;

  // set buffer
  u.set_buffer(packed.size(), &packed[0]);

  // unpack vector sizes
  int ghost_cells_size = 0;
  u >> ghost_cells_size;
  Check(ghost_cells_size > 0);

  // resize vectors
  ghost_cells_d.resize(ghost_cells_size);
  ghost_procs_d.resize(ghost_cells_size);

  // unpack the vectors
  Remember(int64_t last = -1);
  for (size_t i = 0, s = nodes + 1; i < s; ++i) {
    u >> global_offset_d[i];

    Ensure(global_offset_d[i] > last);
    Remember(last = global_offset_d[i]);
  }
  Check(global_offset_d[0] == 0);

  for (int i = 0; i < ghost_cells_size; ++i) {
    u >> ghost_cells_d[i];
    u >> ghost_procs_d[i];
  }

  // set the global number of cells in the base class
  Base::num_global_cells_d = global_offset_d[nodes];

  // set the number of local cells
  num_local_d = static_cast<size_t>(global_offset_d[node + 1] - global_offset_d[node]);

  Ensure(Base::get_parallel_scheme() == "DD");
  Ensure(ghost_cells_d.size() == ghost_procs_d.size());
}

//------------------------------------------------------------------------------------------------//
// PUBLIC INTERFACE FUNCTIONS
//------------------------------------------------------------------------------------------------//
//! Pack up the DD_Topology.
DD_Topology::packed_type DD_Topology::pack() const {
  size_t const nodes = rtt_c4::nodes();
  Require(global_offset_d.size() == nodes + 1);
  Require(ghost_cells_d.size() == ghost_procs_d.size());

  // calculate the size
  size_t const size = sizeof(int64_t) * global_offset_d.size() +
                      sizeof(int) * (1 + ghost_cells_d.size() + ghost_procs_d.size());
  Check(size >= 3 * sizeof(int) + sizeof(int64_t));

  // make packed data
  packed_type packed(size);

  // make a packer and pack the data
  rtt_dsxx::Packer p;
  p.set_buffer(size, &packed[0]);

  // store the scalars
  p << static_cast<int>(ghost_cells_d.size());

  // store the arrays
  for (size_t i = 0, s = nodes + 1; i < s; ++i)
    p << global_offset_d[i];

  for (size_t i = 0, s = ghost_cells_d.size(); i < s; ++i) {
    p << ghost_cells_d[i];
    p << ghost_procs_d[i];
  }

  // return packed topology
  Ensure(&packed[0] + size == p.get_ptr());
  return packed;
}

//------------------------------------------------------------------------------------------------//
//! Print out data in the topology class for diagnostic purposes.
void DD_Topology::print(std::ostream &out) const {
  out << "\n*** DD_TOPOLOGY DATA ***\n------------------------\n"
      << "\nParallel Topology: " << get_parallel_scheme() << "\n------------------------\n\n";
}

} // end namespace rtt_mesh

//------------------------------------------------------------------------------------------------//
// end of DD_Topology.cc
//------------------------------------------------------------------------------------------------//
