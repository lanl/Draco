//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/Range_Finder.hh
 * \author Mike Buksas
 * \date   Thu Feb  6 12:10:56 2003
 * \brief  Header file for Range_finder
 * \note   Copyright (C) 2024-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_Range_finder_hh
#define rtt_dsxx_Range_finder_hh

#include "Assert.hh"
#include "dsxx/device_config.h"
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <ranges>
#include <utility>

namespace rtt_dsxx {

//================================================================================================//
/*!
 * Range_finder
 *
 * These functions locate a value in intervals described by an increasing array.  E.g. Let `v[i]`
 * be an increasing array and `r` a value. These functions look for `i` such that:
 * `v[i] < r < v[i+1]`. The different versions of the function have different behavior in the event
 * that `r = v[i]` for some `i`.
 *
 * The 'left' versions assume the intervals are closed on the left, so that `v[i] <= r < v[i+1]`.
 * Likewise, the 'right' versions assume closure on the right: `v[i] < r <= v[i+1]`.
 *
 * The 'catch end' versions will allow equality for values at the ends of the range when it would
 * not normally be considered a part of the interval. For example if `r=v[0]`, Range_finder_left
 * will consider `r` to be outside the range, but Range_finder_left_catch_end will return `i=0`.
 * Likewise for `r=v[n-1]` (`v`'s last value) Range_finder_right_catch_end will return `n-1` instead
 * of considering the value to be out of range.
 *
 * The functions Range_funder and Range_finder_catch_end take an enumeration parameter
 * `RANGE_DIRECTION` to determine the direction: `LEFT` or `RIGHT`.
 *
 * Values which are out of range will fail a DCB check, but pass silently if DBC is off. You have
 * been warned!
 *
 * This function uses the equal_range algorithm, which returns an iterator pair.  Both iterators
 * will point to the first value _after_ `r` in `v`, unless `v[i]==r` for some `i`, then the first
 * iterator will point to `v[i]` and the second to `v[i+1]`. This enables simple verification of out
 * of bound index detection and a simple computation for the index which works whether or not
 * equality is obtained.
 *
 * \todo Replace signatures with range based versions...
 * \code
 * #include <ranges>
 * int Range_finder_left(std::ranges::input_range auto&& r,
 *                       std::ranges::range_value_t<decltype(r)> value);
 * \endcode
 */
//================================================================================================//

enum class RANGE_DIRECTION : uint8_t { LEFT = 0, RIGHT = 1 };

namespace {

//! Return false if the iterator pair indicate an out of range result.
template <std::forward_iterator It>
bool validate(const std::pair<It, It> &it, const It &begin, const It &end) {
  return !(it.first == end ||   // v > *(end-1)
           it.second == begin); // v < *begin
}

} // namespace

//------------------------------------------------------------------------------------------------//
template <typename It>
int Range_finder_left(It const &begin, It const &end, std::iter_value_t<It> value) {
  auto it = std::equal_range(begin, end, value);
  Require(validate(it, begin, end));

  const int index = static_cast<int>(it.second - begin) - 1;
  Ensure(index >= 0);
  Ensure(begin + index < end);

  return index;
}

//------------------------------------------------------------------------------------------------//
template <typename It>
int Range_finder_right(It const &begin, It const &end, std::iter_value_t<It> value) {
  auto it = std::equal_range(begin, end, value);
  Require(validate(it, begin, end));

  const int index = static_cast<int>(it.first - begin) - 1;
  Ensure(index >= 0);
  Ensure(begin + index < end);

  return index;
}

//------------------------------------------------------------------------------------------------//
template <typename It>
GPU_HOST_DEVICE int Range_finder_left_catch_end(It begin, It end, std::iter_value_t<It> value) {
  const auto it = std::equal_range(begin, end, value);
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
template <typename It>
int Range_finder_right_catch_end(It begin, It end, std::iter_value_t<It> value) {
  auto it = std::equal_range(begin, end, value);
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
template <typename It>
int Range_finder(It const &begin, It const &end, std::iter_value_t<It> value,
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
template <typename It>
int Range_finder_catch_end(It begin, It end, std::iter_value_t<It> value,
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
