//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   linear/rsolv.hh
 * \author Kent Budge
 * \brief  Solve an upper triangular system of equations
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef linear_rsolv_hh
#define linear_rsolv_hh

#include "dsxx/Assert.hh"
#include <ranges>

namespace rtt_linear {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Solve an upper triangular system of equations.
 *
 * \arg \a Vector A random access container like `std::vector`. Must provide `operator[]`
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
template <std::ranges::random_access_range Vector>
  requires std::ranges::sized_range<Vector>
void rsolv(const Vector &R, const std::size_t n, Vector &b) {
  Require(std::ranges::size(R) == n * n);
  Require(std::ranges::size(b) == n);

  for (size_t i = n; i-- > 0;) {
    b[i] /= R[i + n * i];
    for (size_t j = i; j-- > 0;) {
      b[j] -= R[j + n * i] * b[i];
    }
  }
}

} // namespace rtt_linear

#endif // linear_rsolv_hh

//------------------------------------------------------------------------------------------------//
// end of linear/rsolv.hh
//------------------------------------------------------------------------------------------------//
