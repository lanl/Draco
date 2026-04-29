//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/Invert_Comm_Map.cc
 * \author Rob Lowrie
 * \date   Mon Nov 19 10:09:11 2007
 * \brief  Implementation of Invert_Comm_Map
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/Invert_Comm_Map.hh"
#include "dsxx/Assert.hh"
#include "c4/C4_Functions.hh"
#include "c4/C4_Status.hh"
#include "c4/MPI_Traits.hh"
#include <cstddef> // size_t
#include <vector>

namespace {
//------------------------------------------------------------------------------------------------//
void invert_comm_map_body(rtt_c4::Invert_Comm_Map_t const &to_map,
                          rtt_c4::Invert_Comm_Map_t &from_map, const int num_recv) {
  const int my_proc = rtt_c4::node();
  Remember(const int num_procs = rtt_c4::nodes());

  // request handle for the receives
  std::vector<rtt_c4::C4_Req> recvs(num_recv);

  // the number of data elements to be received from each proc.  This data will ultimately be loaded
  // into from_map, once we know the sending proc ids.
  std::vector<size_t> sizes(num_recv);

  // communication tag for sends/recvs
  const int tag = 201;

  // Posts the receives for the data sizes.  We don't yet know the proc numbers sending the data, so
  // use any_source.
  for (int i = 0; i < num_recv; ++i) {
    receive_async(recvs[i], &sizes[i], 1, rtt_c4::any_source, tag);
  }

  from_map.clear(); // empty whatever came in

  // Send the data sizes and take care of on-proc map.
  for (auto it = to_map.begin(); it != to_map.end(); ++it) {
    Require(it->first >= 0);
    Require(it->first < num_procs);
    Require(it->second > 0);
    if (it->first == my_proc) {
      Check(from_map.find(my_proc) == from_map.end());
      // on-proc map
      from_map[my_proc] = it->second;
    } else {
      // we can ignore the request returned, because our send buffers are not shared, and we'll wait
      // on the receives below.
      rtt_c4::send_async(&(it->second), 1, it->first, tag);
    }
  }

  // Wait on the receives and populate the map
  rtt_c4::C4_Status status{};
  for (int i = 0; i < num_recv; ++i) {
    recvs[i].wait(&status);
    Check(status.get_message_size() == sizeof(size_t));
    const int proc = status.get_source();
    Check(proc >= 0);
    Check(proc < num_procs);

    // proc should not yet exist in map
    Check(from_map.find(proc) == from_map.end());

    from_map[proc] = sizes[i];
  }

  return;
}

} // namespace

namespace rtt_c4 {

#ifdef C4_MPI
int get_num_recv(const Invert_Comm_Map_t &to_map) {
  // number of procs we will receive data from MPI all reduce version of get_num_recv function
  const int num_procs = rtt_c4::nodes();
  const int my_proc = rtt_c4::node();

  // ARL: note number of ranks * ranks per node * 4 bytes memory on node, old method used MPI
  // windows, which uses less memory but seems to hang on cray mpich systems
  std::vector<int> num_recv_global(num_procs, 0);

  for (const auto &it : to_map) {
    Require(it.first >= 0);
    Require(it.first < num_procs);
    // treat only non-local sends
    if (it.first != my_proc) {
      num_recv_global[it.first]++;
    }
  }

  global_sum(num_recv_global.data(), num_recv_global.size()); // NOLINT(misc-include-cleaner)

  Ensure(num_recv_global[my_proc] >= 0 && num_recv_global[my_proc] < num_procs);
  return num_recv_global[my_proc];
}

int get_num_recv_rma(const Invert_Comm_Map_t &to_map) {
  Remember(const int num_procs = rtt_c4::nodes());
  const int my_proc = rtt_c4::node();
  const int one = 1;
  int num_recv(0);

  //NOLINTBEGIN(misc-include-cleaner)

  // Create the RMA memory windows for each value
  MPI_Win win;
  MPI_Win_create(&num_recv, 1 * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

  // Assertion value for fences.  Currently, we effectively don't set anything (zero).
  const int fence_assert = 0;

  // Accumulate the local and remote data values
  MPI_Win_fence(fence_assert, win);
  for (const auto &it : to_map) {
    Require(it.first >= 0);
    Require(it.first < num_procs);
    if (it.first != my_proc) { // treat only non-local sends
      // ...increment the remote number of receives
      MPI_Accumulate(&one, 1, MPI_Traits<int>::element_type(), it.first, 0, 1,
                     MPI_Traits<int>::element_type(), MPI_SUM, win);
    }
  }
  MPI_Win_fence(fence_assert, win);
  MPI_Win_free(&win);

  //NOLINTEND(misc-include-cleaner)

  Ensure(num_recv >= 0 && num_recv < num_procs);
  return num_recv;
}

#elif defined(C4_SCALAR)
// SCALAR version of get_num_recv and get_num_recv_rma
int get_num_recv(const Invert_Comm_Map_t & /*to_map*/) { return 0; }
int get_num_recv_rma(const Invert_Comm_Map_t & /*to_map*/) { return 0; }
#else
// error if C4 mode not recognized
int get_num_recv(const Invert_Comm_Map_t & /*to_map*/) {
  Insist(0, "get_num_recv not implemented for this communication type!";)
}
int get_num_recv_rma(const Invert_Comm_Map_t & /*to_map*/) {
  Insist(0, "get_num_recv not implemented for this communication type!";)
}
#endif // ifdef C4_MPI

//------------------------------------------------------------------------------------------------//
void invert_comm_map(Invert_Comm_Map_t const &to_map, Invert_Comm_Map_t &from_map) {
  const int num_recv = get_num_recv(to_map);
  invert_comm_map_body(to_map, from_map, num_recv);
  return;
}

//------------------------------------------------------------------------------------------------//
void invert_comm_map_rma(Invert_Comm_Map_t const &to_map, Invert_Comm_Map_t &from_map) {
  const int num_recv = get_num_recv_rma(to_map);
  invert_comm_map_body(to_map, from_map, num_recv);
  return;
}

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of Invert_Comm_Map.cc
//------------------------------------------------------------------------------------------------//
