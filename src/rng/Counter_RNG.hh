//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   rng/Counter_RNG.hh
 * \author Peter Ahrens
 * \date   Fri Aug 3 16:53:23 2012
 * \brief  Declaration of class Counter_RNG.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC., All rights reserved */
//------------------------------------------------------------------------------------------------//

#ifndef Counter_RNG_hh
#define Counter_RNG_hh

#include "device/config.h"
#include "rng/config.h"

#ifdef _MSC_FULL_VER
// - 4521: Engines have multiple copy constructors, quite legal C++, disable MSVC complaint.
// - 4244: possible loss of data when converting between int types.
// - 4204: nonstandard extension used - non-constant aggregate initializer
// - 4127: conditional expression is constant
#pragma warning(push)
#pragma warning(disable : 4521 4244 4204 4127)
#endif

#if defined(__ICC)
// Suppress Intel's "unrecognized preprocessor directive" warning, triggered by use of #warning in
// Random123/features/sse.h.
#pragma warning disable 11
#endif

#if defined(__GNUC__) && !defined(__clang__)

#pragma GCC diagnostic push
#if (DBS_GNUC_VERSION >= 70000)
#pragma GCC diagnostic ignored "-Wexpansion-to-defined"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

#if defined(__clang__) && !defined(__ibmxl__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wexpansion-to-defined"
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#pragma clang diagnostic ignored "-Wimplicit-int-conversion"
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#pragma clang diagnostic ignored "-Wextra-semi"
#endif

#if defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-identifier"
#endif

#include "Random123/threefry.h"
#include "uniform.hpp"

#if defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif

#if defined(__clang__) && !defined(__ibmxl__)
// Restore clang diagnostics to previous state.
#pragma clang diagnostic pop
#endif

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

#ifdef _MSC_FULL_VER
#pragma warning(pop)
#endif

#include "ds++/Data_Table.hh"
#include <algorithm>

namespace rtt_rng {

// Forward declaration.
class Counter_RNG;

// Select a particular counter-based random number generator from Random123.
using CBRNG = r123::Threefry2x64;

// Counter and key types.
using ctr_type = CBRNG::ctr_type;
using key_type = CBRNG::key_type;

#define CBRNG_DATA_SIZE 4

namespace // anonymous
{

//------------------------------------------------------------------------------------------------//
/*! \brief Generate a nearly-unique identifier.
 *
 * Given a pointer to RNG state data, this function generates a 64-bit identifier unique to this
 * generator but not to the specific position of its RNG stream.  In other words, the identifier
 * associated with a given generator will not change as random numbers are generated from it.
 * However, this insensitivity to the specific stream position also means that repeated spawning
 * will eventually produce two generators with the same identifier.
 *
 * This function simply applies the chosen counter-based RNG to a shuffled version of the RNG seed,
 * stream number, and spawn indicator and then returns the lower 64 bits of the result.
 */
GPU_HOST_DEVICE
inline uint64_t _get_unique_num(const ctr_type::value_type *const data) {
  CBRNG hash;
  const ctr_type ctr = {{data[3], data[2]}};
  const key_type key = {{data[1] >> 32, 0}};
  const ctr_type result = hash(ctr, key);
  return result[0];
}

//------------------------------------------------------------------------------------------------//
/*! \brief Generate a random double.
 *
 * Given a pointer to RNG state data, this function returns a random double in the open interval
 * (0, 1).  i.e., excluding the endpoints.
 */
GPU_HOST_DEVICE inline double _ran(ctr_type::value_type *const data) {
  CBRNG rng;

  // Assemble a counter from the first two elements in data.
  ctr_type ctr = {{data[0], data[1]}};

  // Assemble a key from the last two elements in data.
  const key_type key = {{data[2], data[3]}};

  // Invoke the counter-based rng.
  const ctr_type result = rng(ctr, key);

  // Increment the counter.
  ctr.incr();

  // Copy the updated counter back into data.
  data[0] = ctr[0];
  data[1] = ctr[1];

  // Convert the first 64 bits of the RNG output into a double-precision value in the open interval
  // (0, 1) and return it.
  return r123::u01fixedpt<double, ctr_type::value_type>(result[0]);
}

} // namespace

//================================================================================================//
/*!
 * \class Counter_RNG_Ref
 * \brief A reference to a Counter_RNG.
 *
 * Counter_RNG_Ref provides an interface to a counter-based random number generator from the
 * Random123 library from D. E. Shaw Research
 * (http://www.deshawresearch.com/resources_random123.html).  Unlike Counter_RNG, Counter_RNG_Ref
 * doesn't own its own RNG state (i.e., key and counter); instead, it operates using a data block
 * specified during construction.
 */
//================================================================================================//
class Counter_RNG_Ref {
public:
  //! Constructor.  db and de specify the extents of an RNG state block.
  GPU_HOST_DEVICE
  Counter_RNG_Ref(ctr_type::value_type *const db, ctr_type::value_type *const de) : data(db, de) {
    Require(std::distance(db, de) * sizeof(ctr_type::value_type) ==
            sizeof(ctr_type) + sizeof(key_type));
  }

  //! Return a random double in the open interval (0, 1).
  GPU_HOST_DEVICE
  double ran() const { return _ran(data.access()); }

  //! Spawn a new, independent generator from this reference.
  inline void spawn(Counter_RNG &new_gen) const;

  //! Return the stream number.
  GPU_HOST_DEVICE
  uint64_t get_num() const { return data[2]; }

  //! Return a unique identifier for this generator.
  GPU_HOST_DEVICE
  uint64_t get_unique_num() const { return _get_unique_num(data.access()); }

  //! Is this Counter_RNG_Ref a reference to rng?
  GPU_HOST_DEVICE
  inline bool is_alias_for(Counter_RNG const &rng) const;

private:
  mutable rtt_dsxx::Data_Table<ctr_type::value_type> data;
};

//================================================================================================//
/*!
 * \class Counter_RNG
 * \brief A counter-based random-number generator.
 *
 * Counter_RNG provides an interface to a counter-based random number generator from the Random123
 * library from D. E. Shaw Research (http://www.deshawresearch.com/resources_random123.html).
 *
 * Counter_RNG_Ref is a friend of Counter_RNG because spawning a new generator modifies both the
 * parent and the child generator in ways that should not be exposed through the public interface of
 * Counter_RNG.
 *
 * Similarly, Rnd_Control is a friend of Counter_RNG because initializing a generator requires
 * access to private data that should not be exposed through the public interface.  Rnd_Control
 * takes no responsibility for instantiating Counter_RNGs itself, and since copying Counter_RNGs is
 * disabled (via a private copy constructor), an Rnd_Control must be able to initialize a generator
 * that was instantiated outside of its control.
 */
//================================================================================================//
class Counter_RNG {
  friend class Counter_RNG_Ref;
  friend class Rnd_Control;

public:
  using const_iterator = ctr_type::const_iterator;

  /*! \brief Default constructor.
   *
   * This default constructor is invoked when a client wants to create a Counter_RNG but delegate
   * its initialization to an Rnd_Control object.
   */
  Counter_RNG() {
    Remember(constexpr bool is_data_ok = sizeof(data) == sizeof(ctr_type) + sizeof(key_type));
    Require(is_data_ok);
  }

  //! Construct a Counter_RNG using a seed and stream number.
  GPU_HOST_DEVICE
  Counter_RNG(const uint32_t seed, const uint64_t streamnum) { initialize(seed, streamnum); }

  //! Create a new Counter_RNG from data.
  Counter_RNG(const ctr_type::value_type *const begin, const ctr_type::value_type *const end) {
    Require(std::distance(begin, end) * sizeof(ctr_type::value_type) ==
            sizeof(ctr_type) + sizeof(key_type));

    std::copy(begin, end, data);
  }

  //! Disable copy and move constructor.
  Counter_RNG(const Counter_RNG &) = delete;
  Counter_RNG(const Counter_RNG &&) = delete;

  //! Disable copy- and move- assignment operator.
  Counter_RNG &operator=(const Counter_RNG &) = delete;
  Counter_RNG &operator=(const Counter_RNG &&) = delete;

  //! Return a random double in the interval (0, 1).
  GPU_HOST_DEVICE
  double ran() const { return _ran(data); }

  //! Spawn a new, independent generator from this one.
  void spawn(Counter_RNG &new_gen) const { new_gen._spawn(data); }

  //! Return the stream number.
  GPU_HOST_DEVICE
  uint64_t get_num() const { return data[2]; }

  //! Return a unique identifier for this generator.
  GPU_HOST_DEVICE
  uint64_t get_unique_num() const { return _get_unique_num(data); }

  //! Return an iterator to the beginning of the state block.
  GPU_HOST_DEVICE
  const_iterator begin() const { return data; }

  //! Return an iterator to the end of the state block.
  GPU_HOST_DEVICE
  const_iterator end() const { return data + size(); }

  //! Test for equality.
  bool operator==(Counter_RNG const &rhs) const { return std::equal(begin(), end(), rhs.begin()); }

  //! Test for inequality.
  GPU_HOST_DEVICE
  bool operator!=(Counter_RNG const &rhs) const { return !std::equal(begin(), end(), rhs.begin()); }

  //! Return a Counter_RNG_Ref corresponding to this Counter_RNG.
  GPU_HOST_DEVICE
  Counter_RNG_Ref ref() const { return Counter_RNG_Ref(data, data + size()); } // NOLINT

  //! Return the size of this Counter_RNG.
  GPU_HOST_DEVICE
  size_t size() const { return size_bytes() / sizeof(ctr_type::value_type); }

  //! Return the size of this Counter_RNG in bytes.
  GPU_HOST_DEVICE
  size_t size_bytes() const { return sizeof(data); }

private:
  //! \bug Can this be changed to std::array<ctr_type::value_type,CBRNG_DATA_SIZE> ?
  mutable ctr_type::value_type data[CBRNG_DATA_SIZE]; // NOLINT

  //! Initialize internal state from a seed and stream number.
  GPU_HOST_DEVICE
  inline void initialize(const uint32_t seed, const uint64_t streamnum);

  //! Spawn a new, independent generator from the provided state block.
  inline void _spawn(ctr_type::value_type *const parent_data);
};

//------------------------------------------------------------------------------------------------//
// Implementation
//------------------------------------------------------------------------------------------------//

//! Spawn a new, independent generator from this reference.
inline void Counter_RNG_Ref::spawn(Counter_RNG &new_gen) const { new_gen._spawn(data.access()); }

//------------------------------------------------------------------------------------------------//
//! Is this Counter_RNG_Ref a reference to rng?
inline bool Counter_RNG_Ref::is_alias_for(Counter_RNG const &rng) const {
  return rng.begin() == data.access();
}

//------------------------------------------------------------------------------------------------//
//! \brief Initialize internal state from a seed and stream number.
GPU_HOST_DEVICE
inline void Counter_RNG::initialize(const uint32_t seed, const uint64_t streamnum) {
  // Low bits of the counter.
  data[0] = 0;

  // High bits of the counter; used for the seed.
  data[1] = static_cast<uint64_t>(seed) << 32;

  // Low bits of the key; used for the stream number.
  data[2] = streamnum;

  // High bits of the key; used as a spawn counter.
  data[3] = 0;
}

//------------------------------------------------------------------------------------------------//
/*! \brief Spawn a new, independent generator from the provided state block.
 *
 * To provide parallel reproducibility independent of the number of ranks or threads, the set of
 * generators used in a calculation must be the same regardless of rank or thread identifier or
 * count.  To provide that level of reproducibility, the SPRNG library of RNGs implemented a
 * binary-tree algorithm for subdividing the set of available generators and creating a new
 * generator from any existing generator without communication.  Counter_RNG adopts the same
 * approach.
 *
 * The current Counter_RNG (Threefry2x64) uses 128-bit keys and therefore provides 2^128 possible
 * generators.  Using a 64-bit stream number to subdivide the key space produces 2^64 families of
 * generators, each with 2^64 members.
 *
 * Given 2^M possible generators, arranging them in a binary tree produces a tree of depth M.
 *
 * \verbatim
 *                                      0
 *                                    /   \
 *                                  1       2
 *                                /   \   /   \
 *                                3   4   5   6
 *                               / \ / \ / \ / \
 *
 *                              [...]
 *                                    \
 *                                      N
 *                                    /   \
 *                                 2N+1   2N+2
 * \endverbatim
 *
 * If every root generator has a different stream number, the generators spawned from that root will
 * be independent of the generators spawned from any other root.  With 2^64 possible generators per
 * stream number, each root generator can support 63 spawned generations before any repetition might
 * occur.
 *
 * In addition to providing a fixed number of guaranteed-independent generations from spawning as
 * described above, this implementation tries to maximize the number of independent generators that
 * can be spawned in a row from a single parent by shifting that parent to an unused portion of the
 * key space when it reaches the bottom of the tree.
 *
 * When generator N spawns, this implementation creates a new generator at 2N+2 and shifts the
 * parent generator from N to 2N+1.  Spawning repeatedly from the same parent results in a
 * progression down the left side of the tree rooted at N.  When this process runs out of bits (and
 * would lead to overflow, which would lead to generator reuse), the parent and new generators are
 * instead shifted to the first level in the unused subtree below the first spawned child in the
 * previous descent.  This process repeats, each time shifting to subtrees rooted at the first
 * spawned child in the previous descent, until it has iterated through all available subtrees and
 * must wrap back to 0, the original root of the tree.  Starting from node 0, this process provides
 * \f$\sum_{i=1}^{M-1} i = 2016\f$ generators for \f$M = 64\f$.
 */
inline void Counter_RNG::_spawn(ctr_type::value_type *const parent_data) {
  // Initialize this generator with the seed and stream number from the parent.
  auto seed = static_cast<uint32_t>(parent_data[1] >> 32);
  uint64_t streamnum = parent_data[2];
  initialize(seed, streamnum);

  ctr_type::value_type next_id = parent_data[3];

  // If the child generator would overflow the key...
  if (2 * parent_data[3] + 2 < parent_data[3]) {
    // ... look back up the tree for the parent of the first spawned child; it will be the first
    // even-numbered node...
    while (next_id % 2)
      next_id = (next_id - 1) / 2;

    // ... shift to the right subtree of that original parent...
    next_id = 2 * next_id + 2;

    // ... and wrap back to 0 if we've run out of subtrees.
    if (next_id > parent_data[3])
      next_id = 0;
  }

  // Shift the parent to the left child.
  parent_data[3] = 2 * next_id + 1;

  // Shift this generator to the right child.
  data[3] = parent_data[3] + 1;
}

} // end namespace rtt_rng

#endif

//------------------------------------------------------------------------------------------------//
// end Counter_RNG.hh
//------------------------------------------------------------------------------------------------//
