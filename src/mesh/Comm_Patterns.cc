//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Comm_Patterns.cc
 * \author H. Park
 * \date   Aug. 1 2023
 * \brief  Comm_Patterns implementation file.
 * \note   Copyright (C) 2023-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "mesh/Comm_Patterns.hh"
#include "dsxx/Assert.hh"
#include "c4/C4_Functions.hh"
#include <cstddef> // size_t
#include <cstdint> // INT32_MAX

namespace rtt_mesh {

//------------------------------------------------------------------------------------------------//
// Build The Communication Patterns
//------------------------------------------------------------------------------------------------//
/*
 * \brief Build communication patterns between processor's boundary cells.
 *
 * This function uses the Topology class to determine what processors a processor interacts with
 * across its processor boundary. It calculates two information maps:
 *
 * \arg map of processor/boundary cell data that the processor receives from other processors.
 *
 * \arg map of processor/local cell data that the processor sends out to other processors.
 *
 * The former data map, or field, is a list of boundary cells that live on other processors.  For a
 * processor to perform work, it needs to know where it can get data for its boundary cells.  This
 * field tells it from which processors data will come and into which boundary cells that data is to
 * be placed.  The latter field tells the current processor what data it owns that needs to be sent
 * out to other processors so that they may do their boundary work.  This field requires
 * communication to calculate.  The former field can be calculated without communication.
 *
 * We support two basic topology types: domain replication (replication) and domain decomposition
 * (DD).  If the topology is replication then no data is set because there are no boundary cells.
 */
void Comm_Patterns::calc_patterns(SP_Topology const &topology) {
  Require(topology);

  // empty the query maps
  empty_query_maps();

  // check the topology and launch the correct calculator
  if (topology->get_parallel_scheme() == "replication") {
    // there are no boundary cells so we don't do anything here
    Ensure(send_query_map.empty());
    Ensure(recv_query_map.empty());

    return;
  } else if (topology->get_parallel_scheme() == "DD") {
    // we can use a simple function to calculate the boundary cell maps because we have exact
    // one-to-one communication
    calc_recv_map_DD(topology);

    // do the requisite communication to calculate the send maps after the recv_query_map has been
    // calculated
    calc_send_map_DD(topology);

    Ensure(recv_query_map.size() == send_query_map.size());
    ;
  } else if (topology->get_parallel_scheme() == "DD/replication") {
    // we don't support this yet!
    Insist(0, "We don't support general DD schemes yet!"); // LCOV_EXCL_LINE
  } else {
    Insist(0, "Don't know your topology!"); // LCOV_EXCL_LINE
  }
}

//------------------------------------------------------------------------------------------------//
// Comm_Patterns Implementation
//------------------------------------------------------------------------------------------------//

//! Calculate recv_query_map patterns for DD topologies.
void Comm_Patterns::calc_recv_map_DD(SP_Topology const &topology) {
  using std::pair;

  Require(topology->get_parallel_scheme() == "DD");
  Require(recv_query_map.empty());

  // get number of boundary cells
  size_t const num_ghost_cells = topology->num_ghost_cells();

  // insertion notification for recv_query_map
  pair<iterator, bool> insert_indicator;

  // entries into recv_query_map
  proc_pair entry;

  // global cell and processor indices
  int processor;

  // loop through boundary cells for this processor and determine the processor they live on.
  // Insert the boundary cell index into the recv_query_map.
  for (size_t bc = 1; bc <= num_ghost_cells; ++bc) {
    // find out what processor it's on; for full DD each cell only resides on a single processor
    processor = topology->remote_proc(bc);

    // insert the processor/cell pair into the map
    entry.first = processor;
    Check(bc < INT32_MAX);
    entry.second = sf_int(1, static_cast<int>(bc));
    insert_indicator = recv_query_map.insert(entry);

    // if the data wasn't inserted because this processor key already exists in the map then add the
    // cell to the processor cell data listing
    if (!insert_indicator.second)
      recv_query_map[processor].push_back(static_cast<int>(bc));

    Check(recv_query_map[processor].size() <= num_ghost_cells);
  }

  // the recv_query_map should not be empty because there must be at least one boundary cell
  Ensure(!recv_query_map.empty());

#ifdef ENSURE_ON
  size_t count = 0;
  for (auto const &rq : recv_query_map) {
    count += rq.second.size();
  }
  if (count != num_ghost_cells)
    throw rtt_dsxx::assertion("Improper receive map generated."); // LCOV_EXCL_LINE
#endif
}

//------------------------------------------------------------------------------------------------//
//! Calculate send_query_map patterns for DD topologies.
void Comm_Patterns::calc_send_map_DD(SP_Topology const &topology) {
  using std::pair;
  using std::vector;

  Insist(topology->num_local_cells() > 0,
         "Cannot construct Comm_Patterns with zero cells in a DD mesh");

  Require(topology->get_parallel_scheme() == "DD");
  Require(send_query_map.empty());
  Require(!recv_query_map.empty());

  // number of boundary processors we interact with, it is the number of processors we receive
  // messages from because full DD is one-to-one communication
  size_t const num_bproc = recv_query_map.size();

  // C4 Requests for the number of cells that other processors will query us about
  vector<rtt_c4::C4_Req> size_recv(num_bproc);
  vector<rtt_c4::C4_Req> cells_recv(num_bproc);

  // information that we receive from the processors we communicate with telling us the number of
  // cells that are on this processor that they need information about and what (local) cells those
  // are
  vector<int> number_of_cells(num_bproc);
  vector<int *> list_of_cells(num_bproc);

  // iterators to beginning and end of recv_query_map
  auto recv_begin = recv_query_map.begin();
  auto recv_end = recv_query_map.end();

  // post receives to other processors for the number of cells that they require information about
  // and send out to those processors the number of cells that we need information from them about
  int index = 0;
  int processor = 0;
  //int query_size = 0;
  for (auto i = recv_begin; i != recv_end; ++i) {
    // determine the processor that we are communicating with and the number of our boundary cells
    // that reside on that processor
    processor = i->first;
    int const query_size = static_cast<int>(i->second.size());

    // post receive to that processor
    rtt_c4::receive_async(size_recv[index], &number_of_cells[index], 1, processor, 600);

    // advance index for query to a boundary processor
    index++;

    // send the number of cells that the processor needs to send us
    // information about
    rtt_c4::send(&query_size, 1, processor, 600);
  }
  Check(index == static_cast<int>(recv_query_map.size()));

  // receive the size data from sent from other processors->we only receive one message from each of
  // those processors
  size_t finished = 0;
  while (finished < num_bproc) {
    // loop through size C4 requests and see if the message has completed
    for (size_t i = 0; i < num_bproc; ++i) {
      // check to see if the message has been completed
      if (size_recv[i].complete()) {
        finished++;
        Check(!size_recv[i].inuse());
        Check(!size_recv[i].complete());
      }
    }
  }
  Check(finished == size_recv.size());

  // post receives for the list of cells that we will receive from each processor; send out the list
  // of cells that we require from each processor
  index = 0;
  int ghost_cell = 0;
  for (auto i = recv_begin; i != recv_end; ++i) {
    // determine the processor that we are communicating with and the number of cells that reside on
    // that processor
    processor = i->first;
    int const query_size = static_cast<int>(i->second.size());

    // allocate list of cell data arrays
    list_of_cells[index] = new int[number_of_cells[index]];

    // post receives for the cell lists that we will receive from other processors
    rtt_c4::receive_async(cells_recv[index], list_of_cells[index], number_of_cells[index],
                          processor, 601);

    // advance index for the number of boundary processors we receive from
    index++;

    // send the list of remote local cell indices that this processor needs information about to the
    // communicating processor
    int *query_data = new int[query_size];

    // populate the list with the remote local cell indices
    for (int cell = 0; cell < query_size; ++cell) {
      // figure out the boundary cell that this data will be sent to
      ghost_cell = i->second[cell];
      Check(ghost_cell > 0);
      Check(ghost_cell <= static_cast<int>(topology->num_ghost_cells()));

      // add the remote local cell index to the list
      query_data[cell] = topology->remote_local_cell(ghost_cell);
      Check(query_data[cell] > 0);
      Check(topology->remote_proc(ghost_cell) == processor);
    }

    // send the query list to the processor that will be sending data back
    rtt_c4::send(query_data, query_size, processor, 601);
    delete[] query_data;
  }
  Check(index == static_cast<int>(recv_query_map.size()));

  // information we will need to make the send_query_map
  const_iterator iter;
  proc_pair entry;
  pair<iterator, bool> insert_indicator;

  // receive the cell lists from the processors that we communicate with and build the
  // send_query_map
  finished = 0;
  while (finished < cells_recv.size()) {
    // assign the iterator to the recv_query_map to the beginning location
    iter = recv_begin;

    // loop through cell list C4 requests and see if the message is complete
    for (size_t i = 0; i < cells_recv.size(); ++i) {
      // determine the processor that we are communicating with and make it the key entry for the
      // send_query_map; because this is a DD topology the send_query_map and recv_query_map have
      // the same processor keys
      entry.first = iter->first;

      // see if the message is complete
      if (cells_recv[i].complete()) {
        // indicate that the message is complete
        finished++;
        Check(!cells_recv[i].inuse());
        Check(!cells_recv[i].complete());

        // write the list of cells to the entry
        entry.second.resize(number_of_cells[i]);
        for (size_t cell = 0; cell < entry.second.size(); cell++)
          entry.second[cell] = list_of_cells[i][cell];
        delete[] list_of_cells[i];

        // add the entry to the send_query_map
        insert_indicator = send_query_map.insert(entry);
        Check(insert_indicator.second);
      }

      // advance processor map
      iter++;
    }
    Check(iter == recv_end);
  }
  Check(finished == cells_recv.size());

  // the send_query_map should not be empty because there must be at least one boundary cell
  Ensure(!send_query_map.empty());

  // get all processors to the same place before proceeding
  rtt_c4::global_barrier();
}

//------------------------------------------------------------------------------------------------//
//! Empty the query maps.
void Comm_Patterns::empty_query_maps() {
  // empty send_query_map
  if (!send_query_map.empty())
    send_query_map.erase(send_query_map.begin(), send_query_map.end());

  // empty recv_query_map
  if (!recv_query_map.empty())
    recv_query_map.erase(recv_query_map.begin(), recv_query_map.end());

  Ensure(send_query_map.empty());
  Ensure(recv_query_map.empty());
}

//------------------------------------------------------------------------------------------------//
// Comm_Patterns Interface
//------------------------------------------------------------------------------------------------//

/*!
 * \brief Determine if the Comm_Patterns are set.
 *
 * The overloaded boolean operator determines if the Comm_Patterns are set or not.  It has no way of
 * determining if the data is correct.
 */
Comm_Patterns::operator bool() const { return !send_query_map.empty() && !recv_query_map.empty(); }

} // end namespace rtt_mesh

//------------------------------------------------------------------------------------------------//
// end of Comm_Patterns.cc
//------------------------------------------------------------------------------------------------//
