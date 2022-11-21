//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   rng/Subrandom_Generator.hh
 * \author Kent Budge
 * \brief  Definition of class Subrandom_Generator
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rng_Subrandom_Generator_hh
#define rng_Subrandom_Generator_hh

#include "ds++/config.h"

namespace rtt_rng {

//================================================================================================//
/*!
 * \class Subrandom_Generator
 *
 * \brief Abstract class representing various generators for a sequence of subrandom vectors.
 *
 * The name is something of a misnomer, since only the Sobol_Sequence child class actually generates
 * a subrandom sequence. The others generate either pseudorandom sequences or nonrandom sequences.
 *
 * Conceptually, the class computes a sequence of subrandom vectors of the desired dimension. The
 * current vector is accessed an element at a time, which is appropriate for Markov chains where the
 * chain might terminate early. Thus each Markov chain is represented by a different subrandom
 * vector and each event in the chain by a different element of the subrandom vector.  Note that
 * this implies that we have defined a maximum permissible length for the Markov chain.
 */
//================================================================================================//

class Subrandom_Generator {
public:
  // CREATORS

  virtual ~Subrandom_Generator() = default;
  Subrandom_Generator(Subrandom_Generator const &rhs) = delete;
  Subrandom_Generator(Subrandom_Generator &&rhs) noexcept = delete;
  Subrandom_Generator &operator=(Subrandom_Generator const &rhs) = delete;
  Subrandom_Generator &operator=(Subrandom_Generator &&rhs) noexcept = delete;

  //! Advance sequence to the next vector.
  virtual void shift_vector() = 0;

  //! Get the next element in the current vector.
  virtual double shift() = 0;

  // ACCESSORS

  unsigned count() const { return count_; }

protected:
  /*!
   * \brief Constructor for Subrandom_Generator.
   * \param[in] count Dimension of the vector calculated by this generator.
   */
  explicit Subrandom_Generator(unsigned count) : count_(count), element_(0) { /* empty */
  }

  // MANIPULATORS

  // DATA

  unsigned count_;   // Current vector of the sequence
  unsigned element_; // Current element of the current vector
};

} // end namespace rtt_rng

#endif // rng_Subrandom_Generator_hh

//------------------------------------------------------------------------------------------------//
// end of rng/Subrandom_Generator.hh
//------------------------------------------------------------------------------------------------//
