//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/dbc.hh
 * \author Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  Extensions to the STL algorithm library
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved.
 *
 * This header defines several function templates that perform common numerical operations not
 * standardized in the STL algorithm header. It also defines some useful STL-style predicates. These
 * predicates are particularly useful for writing Design by Contract assertions.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_dbc_hh
#define rtt_dsxx_dbc_hh

namespace rtt_dsxx {

//! Check whether a sequence is monotonically increasing.
template <typename Forward_Iterator>
bool is_monotonic_increasing(Forward_Iterator first, Forward_Iterator const last);

//! Check whether a sequence is strictly monotonically increasing.
template <typename Forward_Iterator>
bool is_strict_monotonic_increasing(Forward_Iterator first, Forward_Iterator const last);

//! Check whether a sequence is strictly monotonically decreasing.
template <typename Forward_Iterator>
bool is_strict_monotonic_decreasing(Forward_Iterator first, Forward_Iterator const last);

//! Check whether a matrix is symmetric.
template <typename Random_Access_Container>
bool is_symmetric_matrix(Random_Access_Container const &A, unsigned const n,
                         double const tolerance = 1.0e-12);

//! Return the positive difference of the arguments.
template <typename Ordered_Group_Element>
inline Ordered_Group_Element dim(Ordered_Group_Element a, Ordered_Group_Element b);

} // namespace rtt_dsxx

// Use implicit instantiation for these templatized functions
#include "dbc.i.hh"

#endif // rtt_dsxx_dbc_hh
//------------------------------------------------------------------------------------------------//
// end of dbc.hh
//------------------------------------------------------------------------------------------------//
