//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/Invert_Comm_Map.cc
 * \author Rob Lowrie
 * \date   Mon Nov 19 10:09:11 2007
 * \brief  Implementation of Invert_Comm_Map
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Invert_Comm_Map.hh"
#include "MPI_Traits.hh"
#include "ds++/Assert.hh"
#include <vector>

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
// MPI version of get_num_recv()
#ifdef C4_MPI
int get_num_recv(Invert_Comm_Map_t::const_iterator first, Invert_Comm_Map_t::const_iterator last) {
  const int my_proc = rtt_c4::node();
  Remember(const int num_procs = rtt_c4::nodes());
  const int one = 1;
  int num_recv(0); // return value

  // Create the RMA memory windows for each value
  MPI_Win win;
  MPI_Win_create(&num_recv, 1 * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

  // Assertion value for fences.  Currently, we effectively don't set anything (zero).
  const int fence_assert = 0;

  // Accumulate the local and remote data values
  MPI_Win_fence(fence_assert, win);
  for (auto it = first; it != last; ++it) {
    Require(it->first >= 0);
    Require(it->first < num_procs);
    if (it->first != my_proc) { // treat only non-local sends
      // ...increment the remote number of receives
      MPI_Accumulate(&one, 1, MPI_Traits<int>::element_type(), it->first, 0, 1,
                     MPI_Traits<int>::element_type(), MPI_SUM, win);
    }
  }
  MPI_Win_fence(fence_assert, win);
  MPI_Win_free(&win);

  Ensure(num_recv >= 0 && num_recv < num_procs);
  return num_recv;
}
//------------------------------------------------------------------------------------------------//
// SCALAR version of get_num_recv
#elif defined(C4_SCALAR)
int get_num_recv(Invert_Comm_Map_t::const_iterator /*first*/,
                 Invert_Comm_Map_t::const_iterator /*last*/) {
  return 0;
}
#else
//------------------------------------------------------------------------------------------------//
// Default version of get_num_recv, which throws an error.
int get_num_recv(Invert_Comm_Map_t::const_iterator /*first*/,
                 Invert_Comm_Map_t::const_iterator /*last*/) {
  Insist(0, "get_num_recv not implemented for this communication type!");
}
#endif // ifdef C4_MPI

//------------------------------------------------------------------------------------------------//
void invert_comm_map(Invert_Comm_Map_t const &to_map, Invert_Comm_Map_t &from_map) {
  const int my_proc = rtt_c4::node();
  Remember(const int num_procs = rtt_c4::nodes());

  // number of procs we will receive data from
  const int num_recv = get_num_recv(to_map.begin(), to_map.end());

  // request handle for the receives
  std::vector<C4_Req> recvs(num_recv);

  // the number of data elements to be received from each proc.  This data will ultimately be loaded
  // into from_map, once we know the sending proc ids.
  std::vector<size_t> sizes(num_recv);

  // communication tag for sends/recvs
  const int tag = 201;

  // Posts the receives for the data sizes.  We don't yet know the proc numbers sending the data, so
  // use any_source.
  for (int i = 0; i < num_recv; ++i) {
    receive_async(recvs[i], &sizes[i], 1, any_source, tag);
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
      send_async(&(it->second), 1, it->first, tag);
    }
  }

  // Wait on the receives and populate the map
  C4_Status status{};
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

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of Invert_Comm_Map.cc
//------------------------------------------------------------------------------------------------//
