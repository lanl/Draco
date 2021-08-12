//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/dbc.i.hh
 * \author Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  Template implementation for dbc
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved.
 *
 * This header defines several function templates that perform common numerical operations not
 * standardized in the STL algorithm header. It also defines some useful STL-style predicates. These
 * predicates are particularly useful for writing Design by Contract assertions.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_dbc_i_hh
#define rtt_dsxx_dbc_i_hh

#include "Soft_Equivalence.hh"
#include <algorithm>
#include <functional>

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Check whether a sequence is monotonically increasing.
 *
 * Checks whether every element in a sequence is less than or equal to the next element of the
 * sequence.  This is particularly useful for Design by Contract assertions that check that a
 * sequence is sorted.
 *
 * \tparam Forward_Iterator A forward iterator whose value type supports \c operator<.
 * \param[in,out] first Points to the first element of the sequence.
 * \param[in]     last Points one element past the end of the sequence.
 *
 * \return \c true if \f$a_i<=a_{i+1}\f$ for all \f$a_i\f$ in the sequence;
 * \c false otherwise.
 */
template <typename Forward_Iterator>
bool is_monotonic_increasing(Forward_Iterator first, Forward_Iterator const last) {
  Forward_Iterator prev = first;
  while (++first != last) {
    if (*first < *prev)
      return false;
    prev = first;
  }
  Ensure(prev != first);
  return true;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Check whether a sequence is strictly monotonically increasing.
 *
 * Checks whether every element in a sequence is less than the next element of the sequence.  This
 * is particularly useful for Design by Contract assertions that check the validity of a table of
 * data.
 *
 * \tparam Forward_Iterator A forward iterator whose value type supports \c operator<.
 * \param[in,out] first Points to the first element of the sequence.
 * \param[in]     last Points one element past the end of the sequence.
 *
 * \return \c true if \f$a_i<a_{i+1}\f$ for all \f$a_i\f$ in the sequence; \c false otherwise.
 */
template <typename Forward_Iterator>
bool is_strict_monotonic_increasing(Forward_Iterator first, Forward_Iterator const last) {
  Forward_Iterator prev = first;
  while (++first != last) {
    if (!(*prev < *first))
      return false;
    prev = first;
  }
  return true;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Check whether a sequence is strictly monotonically decreasing.
 *
 * Checks whether every element in a sequence is greater than the next element of the sequence.
 *
 * \tparam Forward_Iterator A forward iterator whose value type supports \c operator<.
 * \param[in,out] first Points to the first element of the sequence.
 * \param[in]     last Points one element past the end of the sequence.
 *
 * \pre \c last>first
 *
 * \return \c true if \f$a_{i+1}<a_i\f$ for all \f$a_i\f$ in the sequence;
 * \c false otherwise.
 */
template <typename Forward_Iterator>
bool is_strict_monotonic_decreasing(Forward_Iterator first, Forward_Iterator const last) {
  Require(first < last);
  Forward_Iterator prev = first;
  while (++first != last) {
    if (*prev < *first)
      return false;
    prev = first;
  }
  return true;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Check whether a matrix is symmetric.
 *
 * \tparam Random_Access_Container A random access container type.
 * \param[in] A Matrix that is supposed to be symmetric.
 * \param[in] n Rank of the matrix.
 * \param[in] tolerance Tolerance for comparing matrix elements.
 *
 * \pre \c A.size()==n*n
 * \pre \c tolerance>=0.0
 *
 * \return \c true if <code>A[i+n*j]==A[j+n*i]</code> for all \c i and \c j; \c false otherwise.
 */
template <typename Random_Access_Container>
bool is_symmetric_matrix(Random_Access_Container const &A, unsigned const n,
                         double const tolerance) {
  Require(A.size() == n * n);
  Require(tolerance >= 0.0);

  for (unsigned i = 1; i < n; ++i) {
    for (unsigned j = 0; j < i; ++j) {
      if (!rtt_dsxx::soft_equiv(A[i + n * j], A[j + n * i], tolerance)) {
        return false;
      }
    }
  }
  return true;
}

} // namespace rtt_dsxx

#endif // rtt_dsxx_dbc_i_hh
//------------------------------------------------------------------------------------------------//
// end of dbc.i.hh
//------------------------------------------------------------------------------------------------//
