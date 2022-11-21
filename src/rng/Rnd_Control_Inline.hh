//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   rng/Rnd_Control_Inline.hh
 * \author Paul Henning
 * \brief  Rnd_Control header file.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_rng_Rnd_Control_Inline_hh
#define rtt_rng_Rnd_Control_Inline_hh

#include "Counter_RNG.hh"
#include <limits>

namespace rtt_rng {

//================================================================================================//
/*!
 * \class Rnd_Control
 * \brief A stream number manager that can initialize RNGs.
 *
 * Rnd_Control manages application-facing RNG information---a seed used by all generators, and the
 * next available stream number to be used when constructing a new generator.
 */
//================================================================================================//
class Rnd_Control {
private:
  //! Seed for initialization of random number streams.
  const uint32_t d_seed;

  //! Next available stream number.
  uint64_t d_streamnum;

  //! Total number of streams supported.
  const uint64_t d_max_streams;

public:
  //! Constructor.
  explicit Rnd_Control(const uint32_t seed, const uint64_t streamnum = 0,
                       const uint64_t max_streams = std::numeric_limits<uint64_t>::max())
      : d_seed(seed), d_streamnum(streamnum), d_max_streams(max_streams) {
    Require(max_streams > 0);
    Require(streamnum < max_streams);
  }

  //! Return the next available stream number.
  uint64_t get_num() const { return d_streamnum; }

  //! Reset the stream number.
  void set_num(const uint64_t num) {
    Require(num < d_max_streams);

    d_streamnum = num;
  }

  //! Return the seed value.
  uint32_t get_seed() const { return d_seed; }

  //! Return the maximum number of streams allowed.
  uint64_t get_max_streams() const { return d_max_streams; }

  GPU_HOST_DEVICE
  inline void initialize(const uint64_t snum, Counter_RNG &cbrng);
  GPU_HOST_DEVICE
  inline void initialize(Counter_RNG &cbrng);
};

//------------------------------------------------------------------------------------------------//
//! Update the stream number and initialize the Counter_RNG.
GPU_HOST_DEVICE
inline void Rnd_Control::initialize(const uint64_t snum, Counter_RNG &cbrng) {
  Require(snum < d_max_streams);

  // Reset the stream number.
  set_num(snum);

  // Continue initialization.
  initialize(cbrng);
}

//------------------------------------------------------------------------------------------------//
//! Initialize the Counter_RNG with the next available stream number.
GPU_HOST_DEVICE
inline void Rnd_Control::initialize(Counter_RNG &cbrng) {
  Require(d_streamnum < d_max_streams);

  // Initialize the counter-based RNG.
  cbrng.initialize(d_seed, d_streamnum);

  // Advance to the next stream number.
  ++d_streamnum;
}

} // end namespace rtt_rng

#endif // rtt_rng_Rnd_Control_hh

//------------------------------------------------------------------------------------------------//
// end of rng/Rnd_Control.hh
//------------------------------------------------------------------------------------------------//
