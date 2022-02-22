//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/DracoStrings.i.hh
 * \author Kelly G. Thompson <kgt@lanl.gov
 * \date   Wednesday, Aug 23, 2017, 12:48 pm
 * \brief  Enscapulates common string manipulations (implicit template implementation).
 * \note   Copyright (C) 2017-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_DracoStrings_i_hh
#define rtt_dsxx_DracoStrings_i_hh

#include "Assert.hh"
#include <algorithm>
#include <iostream>

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
//! Convert a string into a numeric-type type with error checking.
template <typename T> auto parse_number(std::string const &str, bool verbose) -> T {
  T retval(0);
  try {
    retval = parse_number_impl<T>(str);
  } catch (std::invalid_argument &e) {
    // if no conversion could be performed
    if (verbose)
      std::cerr << "\n==ERROR==\nrtt_dsxx::parse_number:: "
                << "No valid conversion from string to a numeric value could be "
                << "found.\n"
                << "\tstring = \"" << str << "\"\n"
                << std::endl;
    throw e;
  } catch (std::out_of_range & /*error*/) {
    // if the converted value would fall out of the range of the result type or if the underlying
    // function (std::strtol or std::strtoull) sets errno to ERANGE.
    if (verbose)
      std::cerr << "\n==ERROR==\nrtt_dsxx::parse_number:: "
                << "Type conversion from string to a numeric value resulted in "
                << "a value that is out of range.\n"
                << "\tstring = \"" << str << "\"\n"
                << std::endl;
    throw;
  } catch (...) {
    // everything else
    Insist(false, "Unknown failures in call to std::sto[x] from parse_number.");
  }
  return retval;
}

//------------------------------------------------------------------------------------------------//
//! Convert a string into a vector of floating-point or an integral values.
template <typename T>
std::vector<T> string_to_numvec(std::string const &str, std::string const &range_symbols,
                                std::string const &delimiters) {

  // for vector data, first and last char might be some delimiter character.  For example "[1,2,3)"
  // or "{3.3, 4.4}"
  Insist(range_symbols.length() == 0 ||
             (range_symbols.length() == 2 && str[0] == range_symbols[0] &&
              str[str.length() - 1] == range_symbols[range_symbols.length() - 1]),
         "String data is malformed. It is missing { or }.");

  // remove beginning and ending braces and any extra whitespace.
  std::string const tstr = trim(str, range_symbols + " \t");

  // Convert the string into a vector<string>
  std::vector<std::string> vsdata = tokenize(tstr, delimiters);

  // convert from vector<string> to vector<T>
  std::vector<T> retval(vsdata.size());
  std::transform(vsdata.begin(), vsdata.end(), retval.begin(),
                 [](const std::string &val) { return parse_number<T>(val); });

  return retval;
}

} // namespace rtt_dsxx

#endif // namespace rtt_dsxx_DracoStrings_i_hh

//------------------------------------------------------------------------------------------------//
// end of DracoStrings.i.hh
//------------------------------------------------------------------------------------------------//
