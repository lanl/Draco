//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/fnorm.hh
 * \author Kent Budge
 * \date   Wed Aug 11 15:21:38 2004
 * \brief  Find the norm of a set of functions.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_fnorm_hh
#define linear_fnorm_hh

#include "ds++/DracoMath.hh"

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
/*! Find the norm of a set of functions.
 *
 * \arg \a RandomContainer A random access container
 * \arg \a Funtion_N_to_N A type representing a set of N functions of N variables.  This type must
 *         support <code>operator()( RandomContainer const&, RandomContainer &)</code>.
 * \param x       Point at which to evaluate norm.
 * \param fvec    On return, contains the values of the functions at the point.
 * \param vecfunc Functor whose norm is to be computed.
 *
 * \pre \c fvec.size()==x.size()
 */
template <typename RandomContainer, typename Function_N_to_N>
typename RandomContainer::value_type fnorm(const RandomContainer &x, RandomContainer &fvec,
                                           const Function_N_to_N &vecfunc) {
  using rtt_dsxx::conj;

  const size_t n = x.size();

  fvec.resize(n);
  vecfunc(x, fvec);

  typename RandomContainer::value_type sum = 0;
  for (size_t i = 0; i < n; ++i) {
    sum += fvec[i] * conj(fvec[i]);
  }
  Ensure(fvec.size() == n);
  return 0.5 * sum;
}

} // end namespace rtt_linear

#endif // linear_fnorm_hh

//------------------------------------------------------------------------------------------------//
// end of linear/fnorm.hh
//------------------------------------------------------------------------------------------------//
