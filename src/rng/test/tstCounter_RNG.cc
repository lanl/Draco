//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   rng/test/tstCounter_RNG.cc
 * \author Peter Ahrens
 * \date   Fri Aug 3 16:53:23 2012
 * \brief  Counter_RNG tests.
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "rng/Counter_RNG.hh"
#include <set>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_rng;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_equality(UnitTest &ut) {
  // Create a Counter_RNG by specifying a seed and stream number.
  uint32_t seed = 1;
  uint64_t streamnum = 2;
  Counter_RNG rng(seed, streamnum);

  FAIL_IF_NOT(rng.get_num() == streamnum);
  FAIL_IF_NOT(rng.size() == CBRNG_DATA_SIZE);
  FAIL_IF_NOT(rng.size_bytes() == CBRNG_DATA_SIZE * sizeof(uint64_t));

  // Create another Counter_RNG with a different seed.
  seed = 2;
  Counter_RNG rng2(seed, streamnum);

  // rng2's stream number should match rng's, but the two generators should not be identical.
  FAIL_IF_NOT(rng2.get_num() == streamnum);
  FAIL_IF_NOT(rng2.get_num() == rng.get_num());
  FAIL_IF(rng2.get_unique_num() == rng.get_unique_num());
  FAIL_IF(rng2 == rng);

  // Create another Counter_RNG with a different stream number.
  seed = 1;
  streamnum = 3;
  Counter_RNG rng3(seed, streamnum);

  // rng3 should be different from the previous two generators.
  FAIL_IF_NOT(rng3.get_num() == streamnum);
  FAIL_IF(rng3.get_unique_num() == rng.get_unique_num());
  FAIL_IF(rng3.get_unique_num() == rng2.get_unique_num());
  FAIL_IF(rng3 == rng);
  FAIL_IF(rng3 == rng2);

  // Create another Counter_RNG with the original seed and stream number.
  streamnum = 2;
  Counter_RNG rng4(seed, streamnum);

  // rng4 should be equal to rng but different from rng2 and rng3.
  FAIL_IF_NOT(rng4.get_num() == streamnum);
  FAIL_IF_NOT(rng4.get_unique_num() == rng.get_unique_num());
  FAIL_IF(rng4.get_unique_num() == rng2.get_unique_num());
  FAIL_IF(rng4.get_unique_num() == rng3.get_unique_num());
  FAIL_IF_NOT(rng4 == rng);
  FAIL_IF(rng4 == rng2);
  FAIL_IF(rng4 == rng3);

  // Create a Counter_RNG from a data array.
  vector<uint64_t> data(CBRNG_DATA_SIZE);
  data[0] = 1234;
  data[1] = 5678;
  data[2] = 9012;
  data[3] = 3456;
  Counter_RNG rng5(&data[0], &data[0] + CBRNG_DATA_SIZE);

  streamnum = data[2];
  FAIL_IF_NOT(rng5.get_num() == streamnum);
  FAIL_IF(rng5.get_unique_num() == rng.get_unique_num());
  FAIL_IF(rng5.get_unique_num() == rng2.get_unique_num());
  FAIL_IF(rng5.get_unique_num() == rng3.get_unique_num());
  FAIL_IF(rng5.get_unique_num() == rng4.get_unique_num());
  FAIL_IF(rng5 == rng);
  FAIL_IF(rng5 == rng2);
  FAIL_IF(rng5 == rng3);
  FAIL_IF(rng5 == rng4);

  // Create a Counter_RNG from a data array that should match rng and rng4.
  data[0] = 0;
  data[1] = static_cast<uint64_t>(1) << 32U;
  data[2] = 2;
  data[3] = 0;
  Counter_RNG rng6(&data[0], &data[0] + CBRNG_DATA_SIZE);

  streamnum = data[2];
  FAIL_IF_NOT(rng6.get_num() == streamnum);
  FAIL_IF_NOT(rng6.get_unique_num() == rng.get_unique_num());
  FAIL_IF(rng6.get_unique_num() == rng2.get_unique_num());
  FAIL_IF(rng6.get_unique_num() == rng3.get_unique_num());
  FAIL_IF_NOT(rng6.get_unique_num() == rng4.get_unique_num());
  FAIL_IF(rng6.get_unique_num() == rng5.get_unique_num());
  FAIL_IF_NOT(rng6 == rng);
  FAIL_IF(rng6 == rng2);
  FAIL_IF(rng6 == rng3);
  FAIL_IF_NOT(rng6 == rng4);
  FAIL_IF(rng6 == rng5);

// Try to create a Counter_RNG from a data array that's too short.
// 1. Only test exceptions if DbC is enabled.
// 2. However, do not run these tests if no-throw DbC is enabled (DBC & 8)
#ifdef REQUIRE_ON
#if !(DBC & 8)
  bool caught = false;
  try {
    Counter_RNG rng7(&data[0], &data[0] + CBRNG_DATA_SIZE - 1);
  } catch (rtt_dsxx::assertion &err) {
    cout << "Good, caught assertion: " << err.what() << endl;
    caught = true;
  }
  FAIL_IF_NOT(caught);
#endif
#endif

  // Test for equality using iterators.
  FAIL_IF_NOT(std::equal(rng6.begin(), rng6.end(), rng.begin()));
  FAIL_IF(std::equal(rng6.begin(), rng6.end(), rng2.begin()));
  FAIL_IF(std::equal(rng6.begin(), rng6.end(), rng3.begin()));
  FAIL_IF(!std::equal(rng6.begin(), rng6.end(), rng4.begin()));
  FAIL_IF(std::equal(rng6.begin(), rng6.end(), rng5.begin()));

  if (ut.numFails == 0)
    PASSMSG("test_equality passed");
}

//------------------------------------------------------------------------------------------------//
void test_stream(UnitTest &ut) {
  // Create two identical Counter_RNGs.
  uint32_t seed = 0x12121212;
  uint64_t streamnum = 1234;
  Counter_RNG rng(seed, streamnum);
  Counter_RNG rng2(seed, streamnum);

  FAIL_IF_NOT(rng == rng2);

  // Generate a random double (and advance the stream) from rng.
  double x = rng.ran();

  // rng and rng2 should no longer match, but their stream numbers and unique identifiers should be
  // the same.
  FAIL_IF(rng == rng2);
  FAIL_IF_NOT(rng.get_num() == streamnum);
  FAIL_IF_NOT(rng.get_num() == rng2.get_num());
  FAIL_IF_NOT(rng.get_unique_num() == rng2.get_unique_num());

  // Generate a random double (and advance the stream) from rng2.
  double y = rng2.ran();

  // Now rng and rng2 should match again, and the two generated doubles should be identical.
  FAIL_IF_NOT(rng == rng2);
  FAIL_IF_NOT(soft_equiv(x, y));

  // Generate another random double from rng.
  double z = rng.ran();

  // Now they should differ again.
  FAIL_IF(rng == rng2);
  FAIL_IF_NOT(rng.get_num() == streamnum);
  FAIL_IF_NOT(rng.get_num() == rng2.get_num());
  FAIL_IF_NOT(rng.get_unique_num() == rng2.get_unique_num());
  FAIL_IF(soft_equiv(x, z));

  // Create a Counter_RNG from a data array.
  vector<uint64_t> data(CBRNG_DATA_SIZE);
  data[0] = 0;
  data[1] = static_cast<uint64_t>(seed) << 32U;
  data[2] = streamnum;
  data[3] = 0;
  Counter_RNG rng3(&data[0], &data[0] + CBRNG_DATA_SIZE);

  // Initially, rng3 should exactly match neither rng nor rng2, but all three should have the same
  // stream number and "unique" identifier.
  FAIL_IF_NOT(std::equal(rng3.begin(), rng3.end(), data.begin()));
  FAIL_IF(rng3 == rng);
  FAIL_IF(rng3 == rng2);
  FAIL_IF_NOT(rng3.get_num() == streamnum);
  FAIL_IF_NOT(rng3.get_unique_num() == rng.get_unique_num());
  FAIL_IF_NOT(rng3.get_unique_num() == rng2.get_unique_num());

  // Generate a random double from rng3; it should match rng2 but not data afterward.
  double w = rng3.ran();
  FAIL_IF(rng3 == rng);
  FAIL_IF_NOT(rng3 == rng2);
  FAIL_IF(std::equal(rng3.begin(), rng3.end(), data.begin()));
  FAIL_IF_NOT(rng3.get_num() == streamnum);
  FAIL_IF_NOT(rng3.get_unique_num() == rng.get_unique_num());
  FAIL_IF_NOT(rng3.get_unique_num() == rng2.get_unique_num());
  FAIL_IF_NOT(soft_equiv(w, y));

  if (ut.numFails == 0)
    PASSMSG("test_stream passed");
}

//------------------------------------------------------------------------------------------------//
void test_alias(UnitTest &ut) {
  // Create four Counter_RNGs; rng and rng2 are identical, and rng, rng2, and rng3 have the same
  // stream number.
  uint64_t streamnum = 0x20202020;
  Counter_RNG rng(0x1111, streamnum);
  Counter_RNG rng2(0x1111, streamnum);
  Counter_RNG rng3(0x2222, streamnum);
  ++streamnum;
  Counter_RNG rng4(0x3333, streamnum);

  FAIL_IF_NOT(rng.get_num() == rng2.get_num());
  FAIL_IF_NOT(rng.get_num() == rng3.get_num());
  FAIL_IF(rng.get_num() == rng4.get_num());
  FAIL_IF_NOT(rng2.get_num() == rng3.get_num());
  FAIL_IF(rng2.get_num() == rng4.get_num());
  FAIL_IF(rng3.get_num() == rng4.get_num());
  FAIL_IF_NOT(rng.get_unique_num() == rng2.get_unique_num());
  FAIL_IF(rng.get_unique_num() == rng3.get_unique_num());
  FAIL_IF(rng.get_unique_num() == rng4.get_unique_num());
  FAIL_IF(rng2.get_unique_num() == rng3.get_unique_num());
  FAIL_IF(rng2.get_unique_num() == rng4.get_unique_num());
  FAIL_IF(rng3.get_unique_num() == rng4.get_unique_num());
  FAIL_IF_NOT(rng == rng2);
  FAIL_IF(rng == rng3);
  FAIL_IF(rng == rng4);
  FAIL_IF(rng2 == rng3);
  FAIL_IF(rng2 == rng4);
  FAIL_IF(rng3 == rng4);

  // Create a Counter_RNG_Ref from rng.
  Counter_RNG_Ref ref(rng.ref());

  FAIL_IF_NOT(ref.get_num() == rng.get_num());
  FAIL_IF_NOT(ref.get_num() == rng2.get_num());
  FAIL_IF_NOT(ref.get_num() == rng3.get_num());
  FAIL_IF(ref.get_num() == rng4.get_num());
  FAIL_IF_NOT(ref.get_unique_num() == rng.get_unique_num());
  FAIL_IF_NOT(ref.get_unique_num() == rng2.get_unique_num());
  FAIL_IF(ref.get_unique_num() == rng3.get_unique_num());
  FAIL_IF(ref.get_unique_num() == rng4.get_unique_num());
  FAIL_IF_NOT(ref.is_alias_for(rng));
  FAIL_IF(ref.is_alias_for(rng2));
  FAIL_IF(ref.is_alias_for(rng3));
  FAIL_IF(ref.is_alias_for(rng4));

  // Generate a random double (and advance the stream) from rng via ref.
  double x = ref.ran();

  FAIL_IF_NOT(ref.get_num() == rng.get_num());
  FAIL_IF_NOT(ref.get_num() == rng2.get_num());
  FAIL_IF_NOT(ref.get_num() == rng3.get_num());
  FAIL_IF(ref.get_num() == rng4.get_num());
  FAIL_IF_NOT(ref.get_unique_num() == rng.get_unique_num());
  FAIL_IF_NOT(ref.get_unique_num() == rng2.get_unique_num());
  FAIL_IF(ref.get_unique_num() == rng3.get_unique_num());
  FAIL_IF(ref.get_unique_num() == rng4.get_unique_num());
  FAIL_IF_NOT(ref.is_alias_for(rng));
  FAIL_IF(ref.is_alias_for(rng2));
  FAIL_IF(ref.is_alias_for(rng3));
  FAIL_IF(ref.is_alias_for(rng4));

  // Invoking ref.ran should have altered rng; it should still have the same stream number as rng2
  // and rng3, but it should be identical to none of them.
  FAIL_IF_NOT(rng.get_num() == rng2.get_num());
  FAIL_IF_NOT(rng.get_num() == rng3.get_num());
  FAIL_IF(rng.get_num() == rng4.get_num());
  FAIL_IF_NOT(rng.get_unique_num() == rng2.get_unique_num());
  FAIL_IF(rng.get_unique_num() == rng3.get_unique_num());
  FAIL_IF(rng.get_unique_num() == rng4.get_unique_num());
  FAIL_IF(rng == rng2);
  FAIL_IF(rng == rng3);
  FAIL_IF(rng == rng4);

  // Create a bare data array that should match rng.
  vector<uint64_t> data(CBRNG_DATA_SIZE);
  data[0] = 1;
  data[1] = static_cast<uint64_t>(0x1111) << 32U;
  data[2] = 0x20202020;
  data[3] = 0;

  FAIL_IF_NOT(std::equal(rng.begin(), rng.end(), data.begin()));

  // Create a Counter_RNG_Ref from a bare data array.
  data[0] = 0;
  data[1] = static_cast<uint64_t>(0x2222) << 32U;
  data[2] = 0x20202020;
  data[3] = 0;
  Counter_RNG_Ref ref2(&data[0], &data[0] + CBRNG_DATA_SIZE);

  // ref2 should have the same stream number as rng, rng2, and rng3 but shouldn't be an alias for
  // any of them.
  FAIL_IF_NOT(ref2.get_num() == rng.get_num());
  FAIL_IF_NOT(ref2.get_num() == rng2.get_num());
  FAIL_IF_NOT(ref2.get_num() == rng3.get_num());
  FAIL_IF(ref2.get_num() == rng4.get_num());
  FAIL_IF(ref2.get_unique_num() == rng.get_unique_num());
  FAIL_IF(ref2.get_unique_num() == rng2.get_unique_num());
  FAIL_IF_NOT(ref2.get_unique_num() == rng3.get_unique_num());
  FAIL_IF(ref2.get_unique_num() == rng4.get_unique_num());
  FAIL_IF(ref2.is_alias_for(rng));
  FAIL_IF(ref2.is_alias_for(rng2));
  FAIL_IF(ref2.is_alias_for(rng3));
  FAIL_IF(ref2.is_alias_for(rng4));

  // Generate a random double from ref.
  double y = ref2.ran();

  // The underlying data array should have changed.
  FAIL_IF_NOT(data[0] == 1);
  FAIL_IF_NOT(data[1] == static_cast<uint64_t>(0x2222) << 32U);
  FAIL_IF_NOT(data[2] == 0x20202020);
  FAIL_IF_NOT(data[3] == 0);
  FAIL_IF(soft_equiv(y, x));

  // Generate a random double from rng3; it should match the one from ref2.
  double z = rng3.ran();

  FAIL_IF_NOT(soft_equiv(z, y));
  FAIL_IF(soft_equiv(z, x));

// Try to create a Counter_RNG_Ref with a data array that's too short.
// 1. Only test exceptions if DbC is enabled.
// 2. However, do not run these tests if no-throw DbC is enabled (DBC & 8)
#ifdef REQUIRE_ON
#if !(DBC & 8)
  bool caught = false;
  try {
    Counter_RNG_Ref ref3(&data[0], &data[0] + CBRNG_DATA_SIZE - 1);
  } catch (rtt_dsxx::assertion &err) {
    cout << "Good, caught assertion: " << err.what() << endl;
    caught = true;
  }
  FAIL_IF_NOT(caught);
#endif
#endif

  if (ut.numFails == 0)
    PASSMSG("test_alias passed");
}

//------------------------------------------------------------------------------------------------//
void test_rollover(UnitTest &ut) {
  // Create a Counter_RNG with a large counter value.
  vector<uint64_t> data(CBRNG_DATA_SIZE);
  data[0] = 0xfffffffffffffffd;
  data[1] = 1;
  data[2] = 0xabcd;
  data[3] = 0xef00;
  Counter_RNG rng(&data[0], &data[0] + CBRNG_DATA_SIZE);

  // Increment data[0], generate a random double, and compare.
  ++data[0];
  double x = rng.ran();
  FAIL_IF_NOT(std::equal(rng.begin(), rng.end(), data.begin()));

  // ... and again.
  ++data[0];
  double y = rng.ran();
  FAIL_IF(soft_equiv(y, x));
  FAIL_IF_NOT(std::equal(rng.begin(), rng.end(), data.begin()));

  // Generate another random double and verify that the counter has incremented correctly.
  data[0] = 0;
  data[1] = 2;
  double z = rng.ran();
  FAIL_IF(soft_equiv(z, x));
  FAIL_IF(soft_equiv(z, y));
  FAIL_IF_NOT(std::equal(rng.begin(), rng.end(), data.begin()));

  // Repeat the test with a Counter_RNG_Ref.
  data[0] = 0xfffffffffffffffe;
  data[1] = 1;
  Counter_RNG_Ref ref(&data[0], &data[0] + CBRNG_DATA_SIZE);

  double y2 = ref.ran();
  FAIL_IF_NOT(soft_equiv(y2, y));
  FAIL_IF_NOT(data[0] == 0xffffffffffffffff);
  FAIL_IF_NOT(data[1] == 1);
  FAIL_IF_NOT(data[2] == 0xabcd);
  FAIL_IF_NOT(data[3] == 0xef00);

  double z2 = ref.ran();
  FAIL_IF_NOT(soft_equiv(z2, z));
  FAIL_IF_NOT(data[0] == 0);
  FAIL_IF_NOT(data[1] == 2);
  FAIL_IF_NOT(data[2] == 0xabcd);
  FAIL_IF_NOT(data[3] == 0xef00);

  if (ut.numFails == 0)
    PASSMSG("test_rollover passed");
}

//------------------------------------------------------------------------------------------------//
void test_spawn(UnitTest &ut) {
  // Create a generator.
  uint32_t seed = 0xabcdef;
  uint64_t streamnum = 0;
  Counter_RNG rng(seed, streamnum);

  // Spawn a new generator.
  Counter_RNG rng_child1;
  rng.spawn(rng_child1);

  // rng_child1 should have the same stream number as rng but should not be identical to rng.
  FAIL_IF_NOT(rng_child1.get_num() == streamnum);
  FAIL_IF_NOT(rng.get_num() == streamnum);
  FAIL_IF(rng_child1 == rng);

  // Create a reference to rng, and spawn from the reference.
  Counter_RNG_Ref rng_ref(rng.ref());
  FAIL_IF_NOT(rng_ref.is_alias_for(rng));

  Counter_RNG rng_child2;
  rng_ref.spawn(rng_child2);

  // rng_child2 should have the same stream number as rng and rng_child1 but should not be identical
  // to either previous generator.
  FAIL_IF_NOT(rng_child2.get_num() == streamnum);
  FAIL_IF_NOT(rng_child1.get_num() == streamnum);
  FAIL_IF_NOT(rng.get_num() == streamnum);
  FAIL_IF(rng_child2 == rng_child1);
  FAIL_IF(rng_child2 == rng);
  FAIL_IF(rng_child1 == rng);

  // Spawn a generator from rng_child1.
  Counter_RNG rng_grandchild1;
  rng_child1.spawn(rng_grandchild1);

  FAIL_IF_NOT(rng_grandchild1.get_num() == streamnum);
  FAIL_IF_NOT(rng_child2.get_num() == streamnum);
  FAIL_IF_NOT(rng_child1.get_num() == streamnum);
  FAIL_IF_NOT(rng.get_num() == streamnum);
  FAIL_IF(rng_grandchild1 == rng_child2);
  FAIL_IF(rng_grandchild1 == rng_child1);
  FAIL_IF(rng_grandchild1 == rng);
  FAIL_IF(rng_child2 == rng_child1);
  FAIL_IF(rng_child2 == rng);
  FAIL_IF(rng_child1 == rng);

  // Spawn a generator from rng_child2.
  Counter_RNG rng_grandchild2;
  rng_child2.spawn(rng_grandchild2);

  FAIL_IF_NOT(rng_grandchild2.get_num() == streamnum);
  FAIL_IF_NOT(rng_grandchild1.get_num() == streamnum);
  FAIL_IF_NOT(rng_child2.get_num() == streamnum);
  FAIL_IF_NOT(rng_child1.get_num() == streamnum);
  FAIL_IF_NOT(rng.get_num() == streamnum);
  FAIL_IF(rng_grandchild2 == rng_grandchild1);
  FAIL_IF(rng_grandchild2 == rng_child2);
  FAIL_IF(rng_grandchild2 == rng_child1);
  FAIL_IF(rng_grandchild2 == rng);
  FAIL_IF(rng_grandchild1 == rng_child2);
  FAIL_IF(rng_grandchild1 == rng_child1);
  FAIL_IF(rng_grandchild1 == rng);
  FAIL_IF(rng_child2 == rng_child1);
  FAIL_IF(rng_child2 == rng);
  FAIL_IF(rng_child1 == rng);

  // Create two identical generators.
  Counter_RNG original(seed, streamnum);
  Counter_RNG parent(seed, streamnum);

  // Repeatedly spawn from parent.  See how long it takes to create a duplicate of a previous
  // generator.
  set<uint64_t> spawn_id;
  spawn_id.insert(original.begin()[3]);

  unsigned int expected_period = 0;
  for (unsigned int i = 0; i < 8 * sizeof(uint64_t); ++i)
    expected_period += i;

  unsigned int timeout = 16384;
  unsigned int gen = 0;
  for (gen = 0; gen < timeout; ++gen) {
    Counter_RNG child;
    parent.spawn(child);

    // The child generator should always have the same stream number as the parent and the original
    // but should never be identical to either.
    FAIL_IF_NOT(child.get_num() == streamnum);
    FAIL_IF_NOT(parent.get_num() == streamnum);
    FAIL_IF(child == parent);
    FAIL_IF(child == original);
    FAIL_IF(parent == original);

    // Look for both parent and child identifiers in the spawn_id set.
    uint64_t parent_id = parent.begin()[3];
    uint64_t child_id = child.begin()[3];

    if (spawn_id.find(parent_id) != spawn_id.end()) {
      std::cout << "Found parent among previous generators after " << gen << " generations."
                << std::endl;
      break;
    }

    if (spawn_id.find(child_id) != spawn_id.end()) {
      std::cout << "Found child among previous generators after " << gen << " generations."
                << std::endl;
      break;
    }

    // Insert the parent and child identifiers and continue.
    spawn_id.insert(parent_id);
    spawn_id.insert(child_id);
  }
  FAIL_IF_NOT(gen == expected_period);

  // Go again from that parent, this time through a reference.  How long until it repeats this time?
  Counter_RNG_Ref parent_ref(parent.ref());
  FAIL_IF_NOT(parent_ref.is_alias_for(parent));

  spawn_id.clear();
  spawn_id.insert(original.begin()[3]);
  spawn_id.insert(parent.begin()[3]);

  // Start from generation 1, to include the generational difference in this experiment between
  // original and parent.
  for (gen = 1; gen < timeout; ++gen) {
    Counter_RNG child;
    parent_ref.spawn(child);
    FAIL_IF_NOT(parent_ref.is_alias_for(parent));

    // The child generator should always have the same stream number as the parent and the original
    // but should never be identical to either.
    FAIL_IF_NOT(child.get_num() == streamnum);
    FAIL_IF_NOT(parent.get_num() == streamnum);
    FAIL_IF(child == parent);
    FAIL_IF(child == original);
    FAIL_IF(parent == original);

    // Look for both parent and child identifiers in the spawn_id set.
    uint64_t parent_id = parent.begin()[3];
    uint64_t child_id = child.begin()[3];

    if (spawn_id.find(parent_id) != spawn_id.end()) {
      std::cout << "Found parent among previous generators after " << gen << " generations."
                << std::endl;
      break;
    }

    if (spawn_id.find(child_id) != spawn_id.end()) {
      std::cout << "Found child among previous generators after " << gen << " generations."
                << std::endl;
      break;
    }

    // Insert the parent and child identifiers and continue.
    spawn_id.insert(parent_id);
    spawn_id.insert(child_id);
  }
  FAIL_IF_NOT(gen == expected_period);

  // Repeat the experiment, but this time use the first child as the starting parent.
  Counter_RNG child;
  parent.spawn(child);

  FAIL_IF_NOT(child.get_num() == streamnum);
  FAIL_IF_NOT(parent.get_num() == streamnum);
  FAIL_IF(child == parent);
  FAIL_IF(child == original);
  FAIL_IF(parent == original);

  spawn_id.clear();
  spawn_id.insert(original.begin()[3]);
  spawn_id.insert(parent.begin()[3]);
  spawn_id.insert(child.begin()[3]);

  expected_period = 1;
  for (unsigned int i = 0; i < 8 * sizeof(uint64_t) - 2; ++i)
    expected_period += i;

  for (gen = 0; gen < timeout; ++gen) {
    Counter_RNG grandchild;
    child.spawn(grandchild);

    // The grandchild generator should always have the same stream number as its parent,
    // grandparent, and the original generator but should never be identical to any of them.
    FAIL_IF_NOT(grandchild.get_num() == streamnum);
    FAIL_IF_NOT(child.get_num() == streamnum);
    FAIL_IF_NOT(parent.get_num() == streamnum);
    FAIL_IF(grandchild == child);
    FAIL_IF(grandchild == parent);
    FAIL_IF(grandchild == original);

    // Look for both child and grandchild identifiers in the spawn_id set.
    uint64_t child_id = child.begin()[3];
    uint64_t grandchild_id = grandchild.begin()[3];

    if (spawn_id.find(child_id) != spawn_id.end()) {
      std::cout << "Found child among previous generators after " << gen << " generations."
                << std::endl;
      break;
    }

    if (spawn_id.find(grandchild_id) != spawn_id.end()) {
      std::cout << "Found grandchild among previous generators after " << gen << " generations."
                << std::endl;
      break;
    }

    // Insert the child and grandchild identifiers and continue.
    spawn_id.insert(child_id);
    spawn_id.insert(grandchild_id);
  }
  FAIL_IF_NOT(gen == expected_period);

  if (ut.numFails == 0)
    PASSMSG("test_spawn passed");
}

//------------------------------------------------------------------------------------------------//
void test_unique(UnitTest &ut) {
  // Create three identical generators.
  uint32_t seed = 332211;
  uint64_t streamnum = 2468;
  Counter_RNG rng(seed, streamnum);
  Counter_RNG rng2(seed, streamnum);
  Counter_RNG rng3(seed, streamnum);

  Counter_RNG_Ref rng_ref(rng.ref());
  Counter_RNG_Ref rng2_ref(rng2.ref());
  Counter_RNG_Ref rng3_ref(rng3.ref());

  FAIL_IF_NOT(rng == rng2);
  FAIL_IF_NOT(rng == rng3);
  FAIL_IF_NOT(rng.get_unique_num() == rng2.get_unique_num());
  FAIL_IF_NOT(rng.get_unique_num() == rng3.get_unique_num());

  FAIL_IF_NOT(rng_ref.is_alias_for(rng));
  FAIL_IF_NOT(rng2_ref.is_alias_for(rng2));
  FAIL_IF_NOT(rng3_ref.is_alias_for(rng3));

  // Generate some random numbers from rng2.  The stream number and unique number of rng2 should
  // remain the same during this process.
  set<uint64_t> ids;
  ids.insert(rng.get_unique_num());

  for (int i = 0; i < 1000000; ++i) {
    rng2.ran();

    FAIL_IF_NOT(rng2.get_num() == rng.get_num());
    FAIL_IF_NOT(rng2_ref.get_unique_num() == rng2.get_unique_num());
    FAIL_IF(ids.find(rng2.get_unique_num()) == ids.end());
  }

  // Spawn some generators from rng3.  While all spawned generators should have the same stream
  // number as the original parent, their unique numbers should differ.
  unsigned int expected_period = 0;
  for (unsigned int i = 0; i < 8 * sizeof(uint64_t); ++i)
    expected_period += i;

  for (unsigned int i = 0; i < expected_period; ++i) {
    Counter_RNG rng3_child;
    rng3.spawn(rng3_child);

    FAIL_IF_NOT(rng3.get_num() == rng.get_num());
    FAIL_IF_NOT(rng3_child.get_num() == rng.get_num());

    FAIL_IF_NOT(rng_ref.get_unique_num() == rng.get_unique_num());
    FAIL_IF_NOT(rng3_ref.get_unique_num() == rng3.get_unique_num());

    FAIL_IF_NOT(ids.find(rng3.get_unique_num()) == ids.end());
    FAIL_IF_NOT(ids.find(rng3_child.get_unique_num()) == ids.end());

    // Insert new unique identifiers and continue.
    ids.insert(rng3.get_unique_num());
    ids.insert(rng3_child.get_unique_num());
  }

  if (ut.numFails == 0)
    PASSMSG("test_unique passed");
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    test_equality(ut);
    test_stream(ut);
    test_alias(ut);
    test_rollover(ut);
    test_spawn(ut);
    test_unique(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstCounter_RNG.cc
//------------------------------------------------------------------------------------------------//
