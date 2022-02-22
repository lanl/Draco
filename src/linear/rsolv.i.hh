//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/rsolv.i.hh
 * \author Kent Budge
 * \date   Tue Aug 10 13:01:02 2004
 * \brief  Solve an upper triangular system of equations
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_rsolv_i_hh
#define linear_rsolv_i_hh

#include "ds++/Assert.hh"

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Solve an upper triangular system of equations.
 *
 * \arg \a RandomContainer A random access container.
 *
 * \param R Upper triangular matrix
 * \param n Rank of the matrix
 * \param b Right-hand side of the system of equations.  On exit, contains the solution of the
 *          system.
 *
 * \pre \c R.size()==n*n
 * \pre \c b.size()==n
 * \pre The diagonal elements of R are nonzero.
 */

template <class RandomContainer>
void rsolv(const RandomContainer &R, const unsigned n, RandomContainer &b) {
  Require(R.size() == n * n);
  Require(b.size() == n);

  for (int i = n - 1; i >= 0; --i) {
    b[i] /= R[i + n * i];
    for (int j = i - 1; j >= 0; --j) {
      b[j] -= R[j + n * i] * b[i];
    }
  }
}

} // end namespace rtt_linear

#endif // linear_rsolv_i_hh

//------------------------------------------------------------------------------------------------//
// end of linear/rsolv.i.hh
//------------------------------------------------------------------------------------------------//
