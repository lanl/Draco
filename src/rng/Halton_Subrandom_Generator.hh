//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   rng/Halton_Subrandom_Generator.hh
 * \author Kent Budge
 * \brief  Definition of class Halton_Subrandom_Generator
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rng_Halton_Subrandom_Generator_hh
#define rng_Halton_Subrandom_Generator_hh

#include "Halton_Sequence.hh"
#include "Subrandom_Generator.hh"

namespace rtt_rng {

//================================================================================================//
/*!
 * \class Halton_Subrandom_Generator
 * \brief Generator for a sequence of subrandom (pseudorandom) vectors.
 *
 * Actually a bit of a misnomer, since the vectors are not really subrandom. See the documentation
 * for Halton_Sequence for further information and examples.
 *
 * This class generates the components of the vector sequence as individual Halton sequences based
 * on different primes. Thus, a subrandom 3-vector sequence would based on the Halton sequences for
 * 2, 3, and 5. */
//================================================================================================//

class Halton_Subrandom_Generator : public Subrandom_Generator {
public:
  // CREATORS

  //! Normal constructor.
  explicit Halton_Subrandom_Generator(unsigned const count = 1);

  // MANIPULATORS

  //! Advance sequence.
  void shift_vector() override;

  //! Get the next element in the current vector.
  double shift() override;

  // ACCESSORS

private:
  // DATA

  std::vector<Halton_Sequence> sequences_;
};

} // end namespace rtt_rng

#endif // rng_Halton_Subrandom_Generator_hh

//------------------------------------------------------------------------------------------------//
// end of rng/Halton_Subrandom_Generator.hh
//------------------------------------------------------------------------------------------------//
