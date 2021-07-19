//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Index_Set.i.hh
 * \author Mike Buksas
 * \date   Thu Feb  2 10:01:46 2006
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef dsxx_Index_Set_i_hh
#define dsxx_Index_Set_i_hh

#include <functional>
#include <numeric>

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/**
 * \brief Set the size of the Index_Set. Discards old size information
 * \arg sizes Pointer to unsigned integers for the index set sizes.
 */
template <unsigned D, int OFFSET>
void Index_Set<D, OFFSET>::set_size(unsigned const *const dimensions_) {
  std::copy(dimensions_, dimensions_ + D, m_dimensions.begin());
  Require(sizes_okay());
  compute_size();
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Set the size of the Index_Set to a uniform dimension. Discards old size information
 * \arg dimension The uniform dimension of the index set.
 */
template <unsigned D, int OFFSET> void Index_Set<D, OFFSET>::set_size(const unsigned dimension) {
  for (unsigned dim = 0; dim < D; ++dim)
    m_dimensions[dim] = dimension;
  compute_size();
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Comparison routine
 * \arg The Index_Set object to compare to.
 */
template <unsigned D, int OFFSET>
inline bool Index_Set<D, OFFSET>::operator==(const Index_Set &rhs) const {
  if (m_array_size != rhs.m_array_size)
    return false;
  return std::equal(m_dimensions.begin(), m_dimensions.begin() + D, rhs.m_dimensions.begin());
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Make sure the indices are with the range for each dimension
 * \arg iterator An iterator to a range of indices.
 */
template <unsigned D, int OFFSET>
template <typename IT>
bool Index_Set<D, OFFSET>::indices_in_range(IT indices) const {

  int dimension = 0;
  for (unsigned index = 0; index < D; ++index, ++dimension)
    if (!index_in_range(indices[index], dimension))
      return false;

  return true;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return true iff the given index is within the range for the given dimension
 *
 * \arg index     The index value
 * \arg dimension The dimension of the index
 */
template <unsigned D, int OFFSET>
inline bool Index_Set<D, OFFSET>::index_in_range(int index, unsigned dimension) const {
  Check(dimension_okay(dimension));

  return ((index >= OFFSET) && (index < static_cast<int>(m_dimensions[dimension]) + OFFSET));
}

//------------------------------------------------------------------------------------------------//
// IMPLEMENTAION
//------------------------------------------------------------------------------------------------//
template <unsigned D, int OFFSET> inline void Index_Set<D, OFFSET>::compute_size() {

  m_array_size = std::accumulate(m_dimensions.begin(), m_dimensions.end(), 1u, std::multiplies<>());
  Ensure(m_array_size > 0);
}

} // end namespace rtt_dsxx

#endif // dsxx_Index_Set_i_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/Index_Set.i.hh
//------------------------------------------------------------------------------------------------//
