//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Rep_Topology.cc
 * \author H. Park
 * \date   Aug. 1 2023
 * \brief  Rep_Topology implementation file.
 * \note   Copyright (C) 2023-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "mesh/Rep_Topology.hh"
#include "dsxx/Assert.hh"
#include "dsxx/Packing_Utils.hh"
#include "mesh/Topology.hh"
#include <cstdint> // int64_t
#include <ostream>

namespace rtt_mesh {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//------------------------------------------------------------------------------------------------//

Rep_Topology::Rep_Topology(int64_t num_global_cells) : Topology(num_global_cells, "replication") {
  Ensure(Base::get_parallel_scheme() == "replication");
}

//------------------------------------------------------------------------------------------------//
//! Unpacking constructor.
Rep_Topology::Rep_Topology(const packed_type &packed) : Topology(1, "replication") {

  Require(packed.size() == sizeof(int64_t));

  // unpack the data
  rtt_dsxx::Unpacker u;

  // set buffer
  u.set_buffer(packed.size(), &packed[0]);

  // unpack the number of global cells
  int64_t tmp(0);
  u >> tmp;
  Base::num_global_cells_d = tmp;

  Ensure(&packed[0] + packed.size() == u.get_ptr());
  Ensure(Base::get_parallel_scheme() == "replication");
}

//------------------------------------------------------------------------------------------------//
// PUBLIC INTERFACE FUNCTIONS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
//! Pack up the Rep_Topology.
Rep_Topology::packed_type Rep_Topology::pack() const {
  // make a packer
  rtt_dsxx::Packer p;

  // calculate the size of packed data
  int const size = sizeof(int64_t);

  // make a packed data buffer
  packed_type packed(size);

  // set the packer
  p.set_buffer(size, &packed[0]);

  // pack the num global cells
  p << Base::num_global_cells_d;

  Ensure(&packed[0] + size == p.get_ptr());

  // return
  return packed;
}

//------------------------------------------------------------------------------------------------//
//! Print out data in the Rep_Topology class for diagnostic purposes.
void Rep_Topology::print(std::ostream &out) const {
  out << "\n*** REP_TOPOLOGY DATA ***\n-------------------------\n\n"
      << "Parallel Topology: " << Base::get_parallel_scheme() << "\n-------------------------\n";
}

} // namespace rtt_mesh

//------------------------------------------------------------------------------------------------//
// end of Rep_Topology.cc
//------------------------------------------------------------------------------------------------//
