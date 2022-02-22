//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   rng/Halton_Sequence.hh
 * \author Kent Budge
 * \date   Thu Dec 22 13:38:35 2005
 * \brief  Definition of class Halton_Sequence
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rng_Halton_Sequence_hh
#define rng_Halton_Sequence_hh

#include "ds++/config.h"
#include <vector>

namespace rtt_rng {

//================================================================================================//
/*!
 * \class Halton_Sequence
 * \brief Generator for a sequence of subrandom (pseudorandom) numbers.
 *
 * The Halton sequence is similar to the Sobol' sequence, but unlike the Sobol' sequence, the Halton
 * sequence does not decorrelate successive iterates by using a gray code. Thus, it is not really
 * subrandom: The approximation of an integral using this sequence to generate quadrature point
 * converges as 1/N like the Sobol' sequence, but it does not do so smoothly for a strongly peaked
 * function. On the other hand, it is easy to generate a Halton sequence of arbitrarily high
 * dimension and to restart a Halton sequence.
 *
 * A Halton sequence exists on the interval (0,1) for every prime integer \f$p\f$. The sequence
 * takes the form
 *
 * \f$1/p, 2/p, ..., (p-1)/p,\f$<br>
 * \f$1/p^2, (1+p)/p^2, (1+2p)/p^2, ... (1+(p-1)p)/p^2,\f$<br>
 * \f$ 2/p^2, (2+p)/p^2, (2+2p)/p^2, ... (2+(p-1)p)/p^2, ... \f$<br>
 * \f$(p-1)/p^2, (p-1+p)/p^2, (p-1+2p)/p^2, ... (p-1+(p-1)p)/p^2, ...\f$<br>
 * \f$1/p^3, (1+p^2)/p^3, (1+2p^2)/p^3, ... (1+(p^2-1)p^2)/p^3,\f$<br>
 * \f$ 2/p^3, (2+p^2)/p^3, (2+2p^2)/p^3, ... (2+(p^2-1)p^2)/p^3, ... \f$<br>
 * \f$(p^2-1)/p^3, (p^2-1+p^2)/p^3, (p^2-1+2p^2)/p^3, ... (p^2-1+(p^2-1)p^2)/p^3, ...\f$
 *
 * For the simplest case \f$p=2\f$ this is
 *
 * \f$1/2, \f$<br>
 * \f$1/4, 3/4,\f$<br>
 * \f$1/8, 3/8, 5/8, 7/8,\f$<br>
 * \f$...\f$
 *
 * while for \f$p=3\f$ this is
 *
 * \f$1/3, 2/3,\f$<br>
 * \f$1/9, 4/9, 7/9,\f$ <br>
 * \f$2/9, 5/9, 8/9,\f$<br>
 * \f$1/81, 10/81, 19/81, 28/81, 37/81, 46/81, 55/81, 64/81, 73/81,\f$<br>
 * \f$2/81, 11/81, 20/81, 29/81, 38/81, 47/81, 56/81, 65/81, 74/81,\f$
 * \f$...\f$
 *
 * These are the same points mapped by the corresponding Sobol' sequence, but they are not scrambled
 * at the lowest level using a gray code. Hence the sequence is easier to generate but not truly
 * subrandom. */
//================================================================================================//

class Halton_Sequence {
public:
  // CREATORS

  //! Default constructor for declaring arrays et al.
  Halton_Sequence() : n_() { /*empty*/
  }

  //! Normal constructor.
  Halton_Sequence(unsigned const base_index, unsigned const count = 1);

  // MANIPULATORS

  //! Get next element of sequence and advance sequence.
  double shift();

  //! Look at next element of sequence without advancing sequence.
  double lookahead() const { return value_; }

  // ACCESSORS

  unsigned base() const { return base_; }
  unsigned count() const { return count_; }

  bool check_class_invariants() const;

private:
  // DATA

  unsigned base_{0};  // Base of the sequence
  unsigned count_{0}; // Current element of the sequence
  double value_{0.0}; // Value of current element of the sequence

  std::vector<unsigned> n_; // Digits in base of current element of sequence
};

} // end namespace rtt_rng

#endif // rng_Halton_Sequence_hh

//------------------------------------------------------------------------------------------------//
// end of rng/Halton_Sequence.hh
//------------------------------------------------------------------------------------------------//
