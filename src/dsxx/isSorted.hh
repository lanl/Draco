//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/isSorted.hh
 * \author Randy M. Roberts
 * \brief  Functions that checks to see if a container is sorted.
 * \note   Copyright (C) 2024-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_isSorted_hh
#define rtt_dsxx_isSorted_hh

namespace rtt_dsxx {

//================================================================================================//
/*!
 * \fn    isSorted
 * \brief Checks to see if elements of [first, last) are sorted, via "<".
 */
//================================================================================================//

template <typename ForwardIterator> bool isSorted(ForwardIterator first, ForwardIterator last) {
  if (first == last)
    return true;

  bool isSorted = true;
  ForwardIterator prev = first;
  ++first;
  while (isSorted && first != last) {
    isSorted = !(*first < *prev);
    prev = first;
    ++first;
  }
  return isSorted;
}

//================================================================================================//
/*!
 * \fn isSorted
 * \brief Checks to see if elements of [first, last) are sorted, via "comp".
 */
template <typename ForwardIterator, typename StrictWeakOrdering>
bool isSorted(ForwardIterator first, ForwardIterator last, StrictWeakOrdering comp) {
  if (first == last)
    return true;

  bool isSorted = true;
  ForwardIterator prev = first;
  ++first;
  while (isSorted && first != last) {
    isSorted = !comp(*first, *prev);
    prev = first;
    ++first;
  }
  return isSorted;
}

} // end namespace rtt_dsxx

#endif // rtt_dsxx_isSorted_hh

//------------------------------------------------------------------------------------------------//
// end of dsxx/isSorted.hh
//------------------------------------------------------------------------------------------------//
