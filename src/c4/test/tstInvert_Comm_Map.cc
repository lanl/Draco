//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstInvert_Comm_Map.cc
 * \author Mike Buksas, Rob Lowrie
 * \date   Mon Nov 19 16:33:08 2007
 * \brief  Tests Invert_Comm_Map
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/Invert_Comm_Map.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"

using namespace std;
using namespace rtt_c4;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void test2(rtt_c4::ParallelUnitTest &ut) {
  size_t const node = rtt_c4::node();
  Invert_Comm_Map_t to_map;

  if (node == 0) {
    to_map[1] = 15;
  } else if (node != 1) {
    // No communication from node 1.
    FAILMSG("Incorrect node for test2.");
  }

  Invert_Comm_Map_t from_map;
  invert_comm_map(to_map, from_map);

  if (node == 0) {
    if (from_map.size() != 0U)
      FAILMSG("Incorrect map size on node 0.");
  } else if (node == 1) {
    if (from_map.size() != 1U)
      FAILMSG("Incorrect size of map on node 1.");
    if (from_map[0] != 15)
      FAILMSG("Incorrect map contents on node 1.");
  }

  if (ut.numFails == 0)
    ut.passes("test2 passes");
  else
    ut.failure("test2 failed");

  return;
}

//------------------------------------------------------------------------------------------------//
void test4(rtt_c4::ParallelUnitTest &ut) {
  size_t const node = rtt_c4::node();
  Invert_Comm_Map_t to_map;

  if (node == 0) {
    to_map[1] = 10;
    to_map[2] = 20;
    to_map[3] = 30;
  } else if (node == 1)
    to_map[0] = 110;
  else if (node == 2)
    to_map[0] = 120;
  else if (node == 3)
    to_map[0] = 130;
  else
    FAILMSG("Incorrect node for test4.");

  Invert_Comm_Map_t from_map;

  invert_comm_map(to_map, from_map);

  if (node == 0) {
    if (from_map.size() != 3U)
      FAILMSG("Incorrect map size on node 0");
    for (int i = 0; i < 3; ++i) {
      if (from_map[i + 1] != 100U + 10U * (i + 1))
        FAILMSG("Incorrent map contents on node 0");
    }
  } else {
    if (from_map.size() != 1U)
      FAILMSG("Incorrect map size.");
    if (from_map[0] != 10U * node)
      FAILMSG("Incorrect map contents.");
  }

  if (ut.numFails == 0)
    ut.passes("test4 passes");
  else
    ut.failure("test4 failed");

  return;
}

//------------------------------------------------------------------------------------------------//
void test_n_to_n(rtt_c4::ParallelUnitTest &ut) {

  const int node = rtt_c4::node();
  const int nodes = rtt_c4::nodes();

  Invert_Comm_Map_t to_map;
  for (int i = 0; i < nodes; ++i)
    to_map[i] = 10U * node + 1U;

  Invert_Comm_Map_t from_map;
  invert_comm_map(to_map, from_map);

  if (static_cast<int>(from_map.size()) != nodes)
    FAILMSG("Incorrect from_map size.");

  for (int i = 0; i < nodes; ++i) {
    if (from_map[i] != 10U * i + 1U)
      FAILMSG("Incorrect data in map.");
  }

  if (ut.numFails == 0)
    ut.passes("test_n_to_n passes");
  else
    ut.failure("test_n_to_n failed");

  return;
}

//------------------------------------------------------------------------------------------------//
void test_cyclic(rtt_c4::ParallelUnitTest &ut) {

  const int node = rtt_c4::node();
  const int nodes = rtt_c4::nodes();

  Invert_Comm_Map_t to_map;
  to_map[(node + 1) % nodes] = 10U;

  Invert_Comm_Map_t from_map;
  invert_comm_map(to_map, from_map);

  if (from_map.size() != 1U)
    FAILMSG("Incorrect map size.");
  if (from_map[(node + nodes - 1) % nodes] != 10U)
    FAILMSG("Incorrect map contents in cyclc test.");

  if (ut.numFails == 0)
    ut.passes("test_cyclic passes");
  else
    ut.failure("test_cyclic failed");

  return;
}

//------------------------------------------------------------------------------------------------//
void test_empty(rtt_c4::ParallelUnitTest &ut) {

  Invert_Comm_Map_t to_map;
  Invert_Comm_Map_t from_map;

  invert_comm_map(to_map, from_map);

  if (from_map.size() != 0U)
    FAILMSG("Incorrect map size in empty test.");

  if (ut.numFails == 0)
    ut.passes("test_empty passes");
  else
    ut.failure("test_empty failed");

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    if (nodes() == 2)
      test2(ut);
    if (nodes() == 4)
      test4(ut);
    test_n_to_n(ut);
    test_cyclic(ut);
    test_empty(ut);
    if (ut.numFails == 0)
      ut.passes("All passed");
    else
      ut.failure("Failed");
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstInvert_Comm_Map.cc
//------------------------------------------------------------------------------------------------//
