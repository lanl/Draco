//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Topology_Builder.cc
 * \author H. Park
 * \date   Aug. 1 2023
 * \brief  Member definitions of class Topology_Builder
 * \note   Copyright (C) 2023-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "mesh/Topology_Builder.hh"
#include "dsxx/Assert.hh"
#include "c4/C4_Functions.hh"
#include "mesh/DD_Topology.hh"
#include "mesh/Rep_Topology.hh"
#include "mesh/Topology.hh"
#include <algorithm>
#include <cstddef> // size_t
#include <cstdint> // uint32_t
#include <map>
#include <memory> // make_shared
#include <numeric>
#include <vector>

namespace rtt_mesh {

//------------------------------------------------------------------------------------------------//
// PUBLIC INTERFACE
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Build a topology.
 *
 * The topology type that is determined by the interface.  Either domain replication
 * (rtt_mesh::Rep_Topology) or domain decomposition (rtt_mesh::DD_Topology) types can be requested.
 * The number of processors must be greater than 1 for rtt_mesh::DD_Topology.
 */
Topology_Builder::SP_Topology Topology_Builder::build_Topology() {
  using rtt_mesh::DD_Topology;
  using rtt_mesh::Rep_Topology;
  using rtt_mesh::Topology;

  // return topology
  SP_Topology topology;

  // determine which topology to build
  if (decomposition_ == Topology::DOMAIN_REPLICATION) {
    // build a Rep_Topology
    topology = std::make_shared<Rep_Topology>(num_global_cells_);
  } else if (decomposition_ == Topology::DOMAIN_DECOMPOSITION) {
    Require(rtt_c4::nodes() > 1);

    // data needed to build DD topology
    Vec_Int cells_per_pe(rtt_c4::nodes(), 0);
    Vec_Int bcells(num_bnd_cells_);
    Vec_Int bprocs(num_bnd_cells_);

    // accumulate the cells per processor on each processor
    build_cell_numbers(cells_per_pe);

    // build the boundary cell data
    build_bnd_data(bcells, bprocs);

    // build the topology
    topology = std::make_shared<DD_Topology>(cells_per_pe, bcells, bprocs);
  } else {
    throw rtt_dsxx::assertion("Invalid topology-type request."); // LCOV_EXCL_LINE
  }

  Ensure(topology);
  Ensure(topology->get_parallel_scheme() == "replication" ||
         topology->get_parallel_scheme() == "DD");
  return topology;
}

//------------------------------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//------------------------------------------------------------------------------------------------//

//! Accumulate the cells-per-processor for every processor on each processor.
void Topology_Builder::build_cell_numbers(Vec_Int &cells_per_pe) {
  using std::accumulate;

  Require(static_cast<int>(cells_per_pe.size()) == rtt_c4::nodes());
  Require(accumulate(cells_per_pe.begin(), cells_per_pe.end(), 0LL) == 0);

  // node info
  int const nodes = rtt_c4::nodes();
  int const node = rtt_c4::node();

  // assign the local number of cells on this processor to the cells_per_pe array
  cells_per_pe[node] = static_cast<int>(num_active_cells_);

  // get the number of cells on each processor using a global reduction
  rtt_c4::global_sum(&cells_per_pe[0], nodes);

  Ensure(static_cast<uint32_t>(cells_per_pe[node]) == num_active_cells_);
  Ensure(accumulate(cells_per_pe.begin(), cells_per_pe.end(), 0LL) == num_global_cells_);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Build boundary cell data needed for DD_Topology instances.
 *
 * The boundary cell data we get from the interface exists on the full mesh. We need to convert it
 * to the active mesh.
 */
void Topology_Builder::build_bnd_data(Vec_Int &bcells, Vec_Int &bprocs) {
  using std::map;
  using std::vector;

  Require(bcells.size() == num_bnd_cells_);
  Require(bprocs.size() == num_bnd_cells_);

  // if the number of full cells is zero then mapping is not required, simply put the data in
  // vectors
  if (num_full_cells_ == 0) {
    std::copy(bnd_procs_, bnd_procs_ + num_bnd_cells_, bprocs.begin());
    std::copy(bnd_cells_, bnd_cells_ + num_bnd_cells_, bcells.begin());
  }

  // otherwise we need to do communication to build the data
  else {
    remap_bnd_data(bcells, bprocs);
  }

  Ensure(bcells.size() == num_bnd_cells_);
  Ensure(bprocs.size() == num_bnd_cells_);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Build boundary cell data when there is a full-active map because of ghost cells.
 *
 * When the mesh has ghost or inactive cells we need to remap the remote, local ids because these
 * are defined on the full mesh.  To do this remapping we need to do communication.
 */
void Topology_Builder::remap_bnd_data(Vec_Int &bcells, Vec_Int &bprocs) {
  using std::map;
  using std::vector;

  Require(bcells.size() == num_bnd_cells_);
  Require(bprocs.size() == num_bnd_cells_);

  // the set of processors and remote cell ids that we communicate with
  map<int, vector<int>> send_full_cells;

  // loop through the boundary cells and place the remote ids we need remapped with the processor
  // that owns them
  for (size_t i = 0; i < num_bnd_cells_; ++i) {
    send_full_cells[bnd_procs_[i]].push_back(bnd_cells_[i]);
  }
  Check(static_cast<int>(send_full_cells.size()) <= rtt_c4::nodes());

  // the number of processors we communicate with
  size_t const comm_procs = send_full_cells.size();

  // post receives for each processor that we communicate with
  vector<rtt_c4::C4_Req> sends(comm_procs);
  vector<rtt_c4::C4_Req> recvs(comm_procs);

  // a map of processor data we receive and their sizes
  map<int, vector<int>> recv_full_cells;
  vector<uint32_t> recv_sizes(comm_procs);

  // the ordered vector of processor ids
  vector<int> ordered_procs(comm_procs);

  // post receives for the sizes we will receive on each processor
  size_t p = 0;

  for (auto &cellset : send_full_cells) {
    Check(p < comm_procs);

    // make the ordered vector of processors
    ordered_procs[p] = cellset.first;

    // post the receive
    rtt_c4::receive_async(recvs[p], &recv_sizes[p], 1, cellset.first, 420);

    // advance the p counter
    ++p;
  }
  Check(p == comm_procs);

  // send out sizes that indicate the number of cells we need information about from each remote
  // processor we communicate with
  for (size_t i = 0; i < comm_procs; ++i) {
    Check(send_full_cells.find(ordered_procs[i]) != send_full_cells.end());

    // the number of cells the remote processor needs to return information about
    auto const size = static_cast<uint32_t>(send_full_cells[ordered_procs[i]].size());

    // send out the sizes to the remote processor
    rtt_c4::send_async(sends[i], &size, 1, ordered_procs[i], 420);

    // make sure the sends complete
    sends[i].wait();
  }

  // receive the sizes from the remote processor
  for (size_t i = 0; i < comm_procs; ++i) {
    // wait on each receive
    recvs[i].wait();
    Check(recv_sizes[i] > 0);

    // size receive full size
    recv_full_cells[ordered_procs[i]] = vector<int>(recv_sizes[i], 0);

    // post receives for the full-mesh cell indices that we will receive from each processor
    rtt_c4::receive_async(recvs[i], &recv_full_cells[ordered_procs[i]][0], recv_sizes[i],
                          ordered_procs[i], 421);
  }

  // send out the full-mesh cell indices that this processor needs remapped
  for (size_t i = 0; i < comm_procs; ++i) {
    Check(send_full_cells.find(ordered_procs[i]) != send_full_cells.end());

    // the number of cells the remote processor needs to return information
    // about
    int const size = static_cast<int>(send_full_cells[ordered_procs[i]].size());

    // send out the cells
    rtt_c4::send_async(sends[i], &send_full_cells[ordered_procs[i]][0], size, ordered_procs[i],
                       421);

    // wait for send to complete
    sends[i].wait();
  }

  // receive all the remote cell lists that need to be remapped
  for (size_t i = 0; i < comm_procs; ++i) {
    Check(send_full_cells.find(ordered_procs[i]) != send_full_cells.end());

    // wait for each complete message
    recvs[i].wait();

    // the number of cells that have been remapped by remote processors that we will receive back
    int const size = static_cast<int>(send_full_cells[ordered_procs[i]].size());

    // post receives for the remapped cells
    rtt_c4::receive_async(recvs[i], &send_full_cells[ordered_procs[i]][0], size, ordered_procs[i],
                          422);
  }

  // now remap cells in the recv_full_cells map using the full-active map on this processor looping
  // through the ordered processors.
  for (size_t i = 0; i < comm_procs; ++i) {
    Check(recv_full_cells.find(ordered_procs[i]) != recv_full_cells.end());

    // make a reference to the cells for easy access
    vector<int> &cells = recv_full_cells[ordered_procs[i]];

    // number of cells to be remapped
    int const size = static_cast<int>(cells.size());

    // loop through the remote cells that live on this processor and remap them
    for (int cell = 0; cell < size; ++cell) {
      Check(cells[cell] <= static_cast<int>(num_full_cells_));

      // remap the cells
      cells[cell] = full2active_[cells[cell] - 1];
      Check(cells[cell] <= static_cast<int>(num_active_cells_));
    }

    // send the remapped cells back to the processor from which they came
    rtt_c4::send_async(sends[i], &recv_full_cells[ordered_procs[i]][0], size, ordered_procs[i],
                       422);

    // wait for send to complete
    sends[i].wait();
  }

  // now receive all the remapped cells
  for (size_t i = 0; i < comm_procs; ++i) {
    // wait for receives to finish
    recvs[i].wait();
  }

  // now make the bcells and bprocs fields from the remapped send_full_cells data; this is the
  // inverse of the first loop where we made the send_full_cells data
  map<int, int> counters;
  for (size_t i = 0; i < num_bnd_cells_; ++i) {
    Check(send_full_cells.find(bnd_procs_[i]) != send_full_cells.end());

    // knowing the processor of the boundary cell we can add the remote-local-active cell ids into
    // the bcell field incrementally
    bcells[i] = send_full_cells[bnd_procs_[i]][counters[bnd_procs_[i]]];

    // increment the counter for this processor
    counters[bnd_procs_[i]]++;

    // map over the processor into bprocs
    bprocs[i] = bnd_procs_[i];
  }

#ifdef CHECK_ON
  int check = 0;
  for (size_t i = 0; i < comm_procs; ++i) {
    check += counters[ordered_procs[i]];
  }
  Check(check == static_cast<int>(num_bnd_cells_));

  for (size_t i = 0; i < comm_procs; ++i) {
    Check(!recvs[i].inuse());
    Check(!sends[i].inuse());
  }
#endif

  // get all processors together before leaving
  rtt_c4::global_barrier();
}

} // namespace rtt_mesh

//------------------------------------------------------------------------------------------------//
// end of mesh/Topology_Builder.cc
//------------------------------------------------------------------------------------------------//
