//--------------------------------------------*-C++-*---------------------------------------------//
/*!
  \file   Norms_Base.hh
  \author Rob Lowrie
  \date   Fri Jan 14 13:01:19 2005
  \brief  Header file for Norms_Base.
  \note   Copyright (C) 2010-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_norms_Norms_Base_hh
#define rtt_norms_Norms_Base_hh

#include "ds++/Assert.hh"
#include <cmath>

namespace rtt_norms {

//================================================================================================//
/*!
  \class Norms_Base

  \brief Base class for all (template) Norms classes.

  This class is not intended to be used by clients.  It allows non-template dependent functionality
  of Norms to be compiled in a single place.
*/
//================================================================================================//
class Norms_Base {
protected:
  // DATA

  // sum of absolute values
  double d_sum_L1 = {0.0};

  // sum of squares of values
  double d_sum_L2 = {0.0};

  // max norm.
  double d_Linf = {0.0};

  // sum of weights
  double d_sum_weights = {0.0};

public:
  // CREATORS

  Norms_Base();

  // Use default copy ctor, dtor, and assignment.

  /// Destructor for Norms_Base.
  virtual ~Norms_Base() = default;

  //! Copy assignment
  Norms_Base &operator=(Norms_Base const &rhs) = default;
  Norms_Base &operator=(Norms_Base &&rhs) = default;

  //! Copy constructor
  Norms_Base(Norms_Base const &rhs) = default;
  Norms_Base(Norms_Base &&rhs) = default;

  // MANIPULATORS

  // Re-initializes the norm values.
  virtual void reset();

  // Equality operator.
  bool operator==(const Norms_Base &n) const;

  // ACCESSORS

  /// Returns the current \f$ L_1 \f$ norm.
  double L1() const {
    Require(d_sum_weights > 0.0);
    return d_sum_L1 / d_sum_weights;
  }

  /// Returns the current \f$ L_2 \f$ norm.
  double L2() const {
    Require(d_sum_weights > 0.0);
    return std::sqrt(d_sum_L2 / d_sum_weights);
  }

  /// Returns the current \f$ L_{\infty} \f$ norm.
  double Linf() const { return d_Linf; }
};

} // namespace rtt_norms

#endif // rtt_norms_Norms_Base_hh

//------------------------------------------------------------------------------------------------//
// end of norms/Norms_Base.hh
//------------------------------------------------------------------------------------------------//
