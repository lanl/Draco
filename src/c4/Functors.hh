//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/Functors.hh
 * \author Thomas M. Evans
 * \date   Fri Sep 30 12:54:09 2005
 * \brief  Functors used in rtt_c4.
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved.
 *
 * Functors for use with TiggsTrace (scatter). */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_c4_Functors_hh
#define rtt_c4_Functors_hh

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief A binary functor to return the maximum of two items
 * \param[in] x, the first operand
 * \param[in] y, the second operand
 * \return the greater value of x or y
 */
template <typename T> struct max {
  T operator()(const T &x, const T &y) const { return (x > y) ? x : y; }
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief A binary functor to find the minimum of two items
 * \param[in] x, the first operand
 * \param[in] y, the second operand
 * \return the lesser value of x or y
 */
template <typename T> struct min {
  T operator()(const T &x, const T &y) const { return (x < y) ? x : y; }
};

} // end namespace rtt_c4

#endif // rtt_c4_Functors_hh

//------------------------------------------------------------------------------------------------//
// end of rtt_c4/Functors.hh
//------------------------------------------------------------------------------------------------//
