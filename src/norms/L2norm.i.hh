//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   norms/L2norm.i.hh
 * \author Kent Budge
 * \date   Tue Sep 18 08:22:09 2007
 * \brief  Member definitions of class L2norm
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef norms_L2norm_i_hh
#define norms_L2norm_i_hh

#include "c4/C4_Functions.hh"
#include "norms/norm.hh"
#include <cmath>
#include <iostream>
#include <numeric>

namespace rtt_norms {

//------------------------------------------------------------------------------------------------//
/*! Helper type for L2norm.
 *
 * \arg \a Field A real type such as float or double.
 */
template <typename Field> double accumulate_norm_(double const init, Field const &x) {
  return init + norm<Field>(x);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \arg \a In An input container type whose elements are real, such as \c vector<double> or \c
 *            list<float>.
 *
 * \param x Container representing a real vector whose norm is desired.
 */
template <typename In> double L2norm(In const &x) {
  double norm = std::accumulate(x.begin(), x.end(), 0.0, accumulate_norm_<typename In::value_type>);

  rtt_c4::global_sum(norm);

  Check(x.size() < UINT_MAX);
  auto xlength(static_cast<unsigned>(x.size()));

  rtt_c4::global_sum(xlength);
  Require(xlength > 0);

  norm = sqrt(norm / xlength);

  Ensure(norm >= 0.0);
  return norm;
}

//------------------------------------------------------------------------------------------------//
/*!
 * This function computes the norm of the difference between two vectors. We have found that this is
 * a surprisingly common operation, and there is advantage to not having to compute the difference
 * vector if all we want is its norm.
 *
 * \arg \a In1 An input container type whose elements are real, such as /c vector<double> or /c
 *             list<float>.
 *
 * \arg \a In2 An input container type whose elements are real, such as /c vector<double> or /c
 *             list<float>.
 *
 * \param x Container representing a real vector.
 *
 * \param y Container representing a real vector.
 */
template <typename In1, typename In2> double L2norm_diff(In1 const &x, In2 const &y) {
  Require(x.size() == y.size());

  auto xi = x.begin();
  auto yi = y.begin();
  // Looping this way avoids restriction to random access containers.
  double norm = 0.0;
  for (; xi != x.end(); ++xi, ++yi) {
    norm += norm_diff<typename In1::value_type, typename In2::value_type>(*xi, *yi);
  }

  rtt_c4::global_sum(norm);

  Check(x.size() < UINT_MAX);
  auto xlength(static_cast<unsigned>(x.size()));

  rtt_c4::global_sum(xlength);
  Require(xlength > 0);

  norm = sqrt(norm / xlength);

  Ensure(norm >= 0.0);
  return norm;
}

} // end namespace rtt_norms

#endif // norms_L2norm_i_hh

//------------------------------------------------------------------------------------------------//
// end of norms/L2norm.i.hh
//------------------------------------------------------------------------------------------------//
