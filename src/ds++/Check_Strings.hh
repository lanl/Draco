//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Check_Strings.hh
 * \author John McGhee
 * \date   Sun Jan 30 14:57:09 2000 *
 * \brief  Provides some utilities to check containers of strings.
 * \note Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 *
 * Functions are provided to examine a container of strings: 1) for the occurrence of certain
 * characters, 2) for the length of each string outside of high and low limits, and 3) to determine
 * if there are any duplicate strings in the container. Return value is a vector of iterators which
 * point to any strings which are selected.
 *
 * \example ds++/test/tstCheck_Strings.cc
 * The following code provides examples of how to use the Check_Strings utilities.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_Check_Strings_hh
#define rtt_dsxx_Check_Strings_hh

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace rtt_dsxx {

// Private functors for internal use by string checking utilities.

struct char_in_string {
  std::string str2;
  explicit char_in_string(std::string s) : str2(std::move(s)) {}
  bool operator()(const std::string &str1) const {
    size_t out = str1.find_first_of(str2);
    return out != std::string::npos;
  }
};

struct string_length_out_of_range {
  size_t low;
  size_t high;
  string_length_out_of_range(const int l, const int h)
      : low(static_cast<size_t>(l)), high(static_cast<size_t>(h)) { /* empty */
  }
  bool operator()(const std::string &str1) const {
    size_t i = str1.size();
    return (i < low) || (i > high);
  }
};

struct strings_equal {
  std::string str2;
  explicit strings_equal(std::string s) : str2(std::move(s)) {}
  bool operator()(const std::string &str1) const { return str1 == str2; }
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Looks through a container of strings to see if any of the strings in the container use any
 *        of the characters specified in the input parameter "match_chars".
 *
 * For example, if you want to create a set of files or directories from a list of strings, this
 * function can be used to check for characters (like "*" for instance) that really shouldn't be
 * used for file or directory names.
 *
 * \param first iterator for the first string in the container to be checked.
 *
 * \param last iterator for the last string in the container to be to be checked plus one.
 *
 * \param match_chars a string containing the characters for which to search.  i.e. - If you want to
 *        check for the presence of "*" or "^" then set std::string match_chars = "*^".
 *
 * \return a vector of iterators for any strings which contain one or more of the match
 *         characters. If the size of this vector is 0 no strings that contained any of the matching
 *         characters were found.
 *
 * \sa Other string checking utilities are available in rtt_dsxx::check_string_lengths, and
 *         rtt_dsxx::check_strings_unique.
 */
template <class IT>
std::vector<IT> check_string_chars(IT const &first, IT const &last,
                                   std::string const &match_chars) {
  std::vector<IT> result_vector;
  if (first == last)
    return result_vector;

  IT out = std::find_if(first, last, char_in_string(match_chars));
  while (out != last) {
    result_vector.push_back(out);
    ++out;
    out = std::find_if(out, last, char_in_string(match_chars));
  }
  return result_vector;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Looks through a container of strings to see if the size of any of the strings in the
 *        container is outside the specified range.
 *
 * For example, if you want to check a container of strings to make sure that all lengths are in the
 * range 1:19 inclusive, set parameter low=1, and parameter high=19.
 *
 * \param first iterator for the first string in the container to be checked.
 *
 * \param last iterator for the last string in the container to be to be checked plus one.
 *
 * \param low lower limit of acceptable string length. If less than zero, it is assumed to be zero.
 *
 * \param high upper limit of acceptable string length
 *
 * \return a vector of iterators for any strings which are of length less than low or greater than
 *         high.  If the size of this vector is 0 no strings with out-of-range lengths were found.
 *
 * \sa Other string checking utilities are available in rtt_dsxx::check_string_chars, and
 *         rtt_dsxx::check_strings_unique.
 */
template <class IT>
std::vector<IT> check_string_lengths(IT const &first, IT const &last, int const low,
                                     int const high) {
  std::vector<IT> result_vector;
  if (first == last)
    return result_vector;
  int lw = low;
  if (lw < 0)
    lw = 0;
  IT out = std::find_if(first, last, string_length_out_of_range(lw, high));
  while (out != last) {
    result_vector.push_back(out);
    ++out;
    out = std::find_if(out, last, string_length_out_of_range(lw, high));
  }
  return result_vector;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Looks through a container of strings to see if there are any duplicates.
 *
 * If a string is duplicated more than once, it will appear in the output vector more than once.
 *
 * \param first iterator for the first string in the container to be checked.
 *
 * \param last iterator for the last string in the container to be to be checked plus one.
 *
 * \return a vector of iterators for any strings which are found to be duplicated.  If the size of
 *         this vector is 0 no duplicate strings were found, that is, all the strings in the
 *         container in the range first:last-1 are unique.
 *
 * \sa Other string checking utilities are available in rtt_dsxx::check_string_chars, and
 *         rtt_dsxx::check_string_lengths.
 */
template <class IT> std::vector<IT> check_strings_unique(IT first, IT const &last) {
  std::vector<IT> result_vector;
  if (first == last)
    return result_vector;
  IT current = first;
  IT next = ++first;
  while (next != last) {
    IT out = std::find_if(next, last, strings_equal(*current));
    if (out != last)
      result_vector.push_back(out);
    ++current;
    ++next;
  }
  return result_vector;
}

} // end namespace rtt_dsxx

#endif // rtt_dsxx_Check_Strings_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/Check_Strings.hh
//------------------------------------------------------------------------------------------------//
