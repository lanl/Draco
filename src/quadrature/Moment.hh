//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Moment.hh
 * \author Kent Budge
 * \brief  Declaration file for the class rtt_quadrature::Moment.
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef quadrature_Moment_hh
#define quadrature_Moment_hh

#include "ds++/DracoMath.hh"

namespace rtt_quadrature {

//================================================================================================//
/*!
 * \class Moment
 * \brief Class containing description of a moment.
 *
 * Each moment is characterized by its degree \f$l\f$ and order \f$m\f$ of the corresponding
 * spherical harmonic. We could represent this as a simple std::pair<unsigned, int>, but these must
 * satisfy certain constraints that are protected by the class representation. We use uppercase L
 * and M for these numbers to avoid the unpleasantness of l (which is nearly indistinguishable from
 * 1) as a symbol name.
 */
//================================================================================================//

class Moment {
public:
  // CREATORS

  //! Create an uninitialized Moment.  This is required by the constructor for vector<Moment>.
  Moment() = default;

  //! Construct a Moment with specified degree and order.
  Moment(unsigned const L, int const M) : L_(L), M_(M) {
    // Chief constraint satisfied by the Moment.
    Require(static_cast<unsigned>(std::abs(M)) <= L);
  }

  // Accessors

  unsigned L() const { return L_; }
  int M() const { return M_; }

private:
  // DATA

  // The data must be kept private in order to protect the invariant.

  unsigned L_{0};
  int M_{0};
};

//------------------------------------------------------------------------------------------------//
//! Test moments for equality
inline bool operator==(Moment const &a, Moment const &b) {
  return a.L() == b.L() && a.M() == b.M();
}

} // end namespace rtt_quadrature

#endif // quadrature_Moment_hh

//------------------------------------------------------------------------------------------------//
// end of quadrature/Moment.hh
//------------------------------------------------------------------------------------------------//
