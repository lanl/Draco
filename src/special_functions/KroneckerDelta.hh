//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/KroneckerDelta.hh
 * \author Kelly Thompson
 * \date   Mon Nov 8 11:17:12 2004
 * \brief  Provide declaration of templatized KroneckerDelta function.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef sf_KroneckerDelta_hh
#define sf_KroneckerDelta_hh

#include "ds++/config.h"
#include <type_traits>

namespace rtt_sf {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief kronecker_delta
 *
 * \param[in] test_value
 * \param[in] offset
 * \return 1 if test_value == offset, otherwise return 0;
 */
template <typename T>
unsigned int
kronecker_delta(T const test_value, T const offset,
                typename std::enable_if<std::is_integral<T>::value>::type * /*unused*/ = nullptr) {
  return (test_value == offset) ? 1 : 0;
}

template <typename T>
unsigned int kronecker_delta(
    T const test_value, T const offset,
    typename std::enable_if<std::is_floating_point<T>::value>::type * /*unused*/ = nullptr) {
  T const eps = std::numeric_limits<T>::epsilon();
  return rtt_dsxx::soft_equiv(test_value, offset, eps) ? 1 : 0;
}

} // end namespace rtt_sf

#endif // sf_KroneckerDelta_hh

//------------------------------------------------------------------------------------------------//
// end of sf/KroneckerDelta.hh
//------------------------------------------------------------------------------------------------//
