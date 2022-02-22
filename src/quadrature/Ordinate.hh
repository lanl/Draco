//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Ordinate.hh
 * \author Kent Budge
 * \date   Tue Dec 21 14:20:03 2004
 * \brief  Declaration file for the class rtt_quadrature::Ordinate.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_quadrature_Ordinate_hh
#define rtt_quadrature_Ordinate_hh

#include "ds++/Soft_Equivalence.hh"
#include "ds++/config.h"

namespace rtt_quadrature {

//================================================================================================//
/*!
 * \class Ordinate
 * \brief Class containing angle cosines and weights for an element of an ordinate set.
 *
 * Provides a container that represents \f$ \mathbf\Omega_m = \mu_m \mathbf e_x + \eta_m \mathbf e_y
 * + \xi_m \mathbf e_z \f$ plus the associated point weight, \f$ w_m \f$. We could represent this as
 * a simple 4-tuple of doubles, but the ordinates must satisfy certain invariants that are protected
 * by the class representation.
 */
//================================================================================================//

class Ordinate {
public:
  // CREATORS

  //! Create an uninitialized Ordinate.  This is required by the constructor for vector<Ordinate>.
  Ordinate() = default;

  //! Construct an Ordinate from the specified vector and weight.
  Ordinate(double const mu, double const eta, double const xi, double const wt)
      : mu_(mu), eta_(eta), xi_(xi), wt_(wt) {
    Require(rtt_dsxx::soft_equiv(mu * mu + eta * eta + xi * xi, 1.0));
  }

  //! Construct a 1D Ordinate from the specified angle and weight.
  inline Ordinate(double const mu, double const wt) : mu_(mu), wt_(wt) {
    Require(mu >= -1.0 && mu <= 1.0);
  }

  // Accessors

  double mu() const { return mu_; }
  double eta() const { return eta_; }
  double xi() const { return xi_; }
  double wt() const { return wt_; }

  void set_wt(double const wt) { wt_ = wt; }

  double const *cosines() const {
    // This is a little crufty, but guaranteed to work according to C++ object layout rules.
    return &mu_;
  }

private:
  // DATA

  // The data must be kept private in order to protect the norm invariant.

  /*!
   * \brief Angle cosines for the ordinate.
   *
   * Do not change the layout of these members! They must be declared in this sequence for cosines()
   * to work as expected! */
  double mu_{0}, eta_{0}, xi_{0};

  //! Quadrature weight for the ordinate.
  double wt_{0};
};

//------------------------------------------------------------------------------------------------//
//! Test ordinates for equality
inline bool operator==(Ordinate const &a, Ordinate const &b) {
  return rtt_dsxx::soft_equiv(a.mu(), b.mu()) && rtt_dsxx::soft_equiv(a.eta(), b.eta()) &&
         rtt_dsxx::soft_equiv(a.xi(), b.xi()) && rtt_dsxx::soft_equiv(a.wt(), b.wt());
}

} // end namespace rtt_quadrature

#endif // rtt_quadrature_Ordinate_hh

//------------------------------------------------------------------------------------------------//
// end of quadrature/Ordinate.hh
//------------------------------------------------------------------------------------------------//
