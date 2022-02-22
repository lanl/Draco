//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Range_Finder.hh
 * \author Mike Buksas
 * \date   Thu Feb  6 12:10:56 2003
 * \brief  Header file for Range_finder
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_Range_finder_hh
#define rtt_dsxx_Range_finder_hh

#include "Assert.hh"
#include <algorithm>
#include <iterator>

namespace rtt_dsxx {

//================================================================================================//
/*!
 * Range_finder
 *
 * These functions locate a value in intervals described by an increasing array.  E.g. Let \c v[i]
 * be an increasing array and \c r a value. These functions look for \c i such that: \code v[i] < r
 * < v[i+1] \endcode. The different versions of the function have different behavior in the event
 * that \code r = v[i] \endcode for some \c i.
 *
 * The 'left' versions assume the intervals are closed on the left, so that \code v[i] <= r < v[i+1]
 * \endcode. Likewise, the 'right' versions assume closure on the right: \code v[i] < r <= v[i+1]
 * \endcode.
 *
 * The 'catch end' versions will allow equality for values at the ends of the range when it would
 * not normally be considered a part of the interval. For example if \c r=v[0], Range_finder_left
 * will consider \c r to be outside the range, but Range_finder_left_catch_end will return \c
 * i=0. Likewise for \c r=v[n-1] (v's last value) Range_finder_right_catch_end will return \c n-1
 * instead of considering the value to be out of range.
 *
 * The functions Range_funder and Range_finder_catch_end take an enumeration parameter \c
 * RANGE_DIRECTION to determine the direction: LEFT or RIGHT.
 *
 * Values which are out of range will fail a DCB check, but pass silently if DBC is off. You have
 * been warned!
 *
 * This function uses the equal_range algorithm, which returns an iterator pair.  Both iterators
 * will point to the first value _after_ \c r in \c v, unless \c v[i]==r for some \c i, then the
 * first iterator will point to \c v[i] and the second to \c v[i+1]. This enables simple
 * verification of out of bound index detection and a simple computation for the index which works
 * whether or not equality is obtained.
 */
//================================================================================================//

enum class RANGE_DIRECTION { LEFT = 0, RIGHT = 1 };

namespace {

// Return false if the iterator pair indicate an out of range result.
template <typename IT> bool validate(const std::pair<IT, IT> &it, IT begin, IT end) {
  return !(it.first == end ||   // v > *(end-1)
           it.second == begin); // v < *begin
}
} // namespace

//------------------------------------------------------------------------------------------------//
template <typename IT>
int Range_finder_left(IT begin, IT end, typename std::iterator_traits<IT>::value_type value) {
  const std::pair<IT, IT> it = std::equal_range(begin, end, value);
  Require(validate(it, begin, end));

  const int index = static_cast<int>(it.second - begin) - 1;
  Ensure(index >= 0);
  Ensure(begin + index < end);

  return index;
}

//------------------------------------------------------------------------------------------------//
template <typename IT>
int Range_finder_right(IT begin, IT end, typename std::iterator_traits<IT>::value_type value) {
  const std::pair<IT, IT> it = std::equal_range(begin, end, value);
  Require(validate(it, begin, end));

  const int index = static_cast<int>(it.first - begin) - 1;
  Ensure(index >= 0);
  Ensure(begin + index < end);

  return index;
}

//------------------------------------------------------------------------------------------------//
template <typename IT>
int Range_finder_left_catch_end(IT begin, IT end,
                                typename std::iterator_traits<IT>::value_type value) {
  const std::pair<IT, IT> it = std::equal_range(begin, end, value);
  Require(validate(it, begin, end));

  const int index = static_cast<int>(it.second - begin) - 1;
  const int range_size = static_cast<int>(end - begin) - 1;

  // If index==range_size, subtract one. This happens when value=*end.
  const int fix_index = std::min(index, range_size - 1);
  Ensure(fix_index >= 0);
  Ensure(begin + fix_index < end);

  return fix_index;
}

//------------------------------------------------------------------------------------------------//
template <typename IT>
int Range_finder_right_catch_end(IT begin, IT end,
                                 typename std::iterator_traits<IT>::value_type value) {
  const std::pair<IT, IT> it = std::equal_range(begin, end, value);
  Require(validate(it, begin, end));

  // Extract the coordinate index (0..n-1)
  const int index = static_cast<int>(it.first - begin) - 1;

  // If we got -1 here, then v=v[0] and we want to catch this end value:
  const int fix_index = std::max(index, 0);
  Ensure(fix_index >= 0);
  Ensure(begin + fix_index < end);

  return fix_index;
}

//------------------------------------------------------------------------------------------------//
// Generic versions.
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
template <typename IT>
int Range_finder(IT begin, IT end, typename std::iterator_traits<IT>::value_type value,
                 RANGE_DIRECTION direction_indicator) {
  Check(direction_indicator == RANGE_DIRECTION::LEFT ||
        direction_indicator == RANGE_DIRECTION::RIGHT);

  if (direction_indicator == RANGE_DIRECTION::LEFT)
    return Range_finder_left(begin, end, value);
  else if (direction_indicator == RANGE_DIRECTION::RIGHT)
    return Range_finder_right(begin, end, value);
  else
    Insist(0, "Invalid direction indicator in Range_finder");

  return -1;
}

//------------------------------------------------------------------------------------------------//
template <typename IT>
int Range_finder_catch_end(IT begin, IT end, typename std::iterator_traits<IT>::value_type value,
                           RANGE_DIRECTION direction_indicator) {
  Check(direction_indicator == RANGE_DIRECTION::LEFT ||
        direction_indicator == RANGE_DIRECTION::RIGHT);

  if (direction_indicator == RANGE_DIRECTION::LEFT)
    return Range_finder_left_catch_end(begin, end, value);
  else if (direction_indicator == RANGE_DIRECTION::RIGHT)
    return Range_finder_right_catch_end(begin, end, value);
  else
    Insist(0, "Invalid direction indicator in Range_finder");

  return -1;
}

} // end namespace rtt_dsxx

#endif // rtt_dsxx_Range_finder_hh

//------------------------------------------------------------------------------------------------//
// end of dsxx/Range_finder.hh
//------------------------------------------------------------------------------------------------//
