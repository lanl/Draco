//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/atomics.hh
 * \author Tim Kelley
 * \date   Thursday, Sept. 6, 2018, 10:50 am
 * \brief  Header file for atomic functions on floatint-point (until C++20)
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_Atomics_hh
#define rtt_dsxx_Atomics_hh

#include <atomic>
#include <type_traits> // std::is_floating_point_v

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/**\brief atomically add arg to a.
 * \tparam FpT a floating point type (integer types in std lib since C++11).
 * \param a    an atomic of type FpT that will be updated.
 * \param arg  quantity to add to a
 * \param m_o  std memory order, default is memory_order_relaxed
 * \return     value of a after update
 *
 * By default, uses memory_order_relaxed, meaning (I think) that other atomic operations on 'a' can
 * be moved before or after this one.
 */
template <typename FpT>
FpT fetch_add(std::atomic<FpT> &a, FpT arg, std::memory_order m_o = std::memory_order_relaxed) {
  static_assert(std::is_floating_point<FpT>::value,
                "Template parameter ought to be floating point, use C++11 std "
                "for integral types");
  FpT expected = a.load();
  FpT to_store = expected + arg;
  while (!a.compare_exchange_weak(expected, to_store, m_o)) {
    expected = a.load();
    to_store = arg + expected;
  }
  return to_store;
}

//------------------------------------------------------------------------------------------------//
/**\brief atomically subtract a from arg.
 * \tparam FpT a floating point type (integer types in std lib since C++11).
 * \param a    an atomic of type FpT that will be updated.
 * \param arg  quantity to subtract from a
 * \param m_o  std memory order, default is memory_order_relaxed
 * \return value of a after update
 *
 * By default, uses memory_order_relaxed, meaning (I think) that other atomic operations on 'a' can
 * be moved before or after this one.
 */
template <class FpT>
FpT fetch_sub(std::atomic<FpT> &a, FpT arg, std::memory_order m_o = std::memory_order_relaxed) {
  static_assert(std::is_floating_point<FpT>::value,
                "Template parameter ought to be floating point, use C++11 std "
                "for integral types");
  FpT expected = a.load();
  FpT to_store = expected - arg;
  while (!a.compare_exchange_weak(expected, to_store, m_o)) {
    expected = a.load();
    to_store = expected - arg;
  }
  return to_store;
}

} // namespace rtt_dsxx

#endif // rtt_dsxx_Atomics_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/atomics.hh
//------------------------------------------------------------------------------------------------//
