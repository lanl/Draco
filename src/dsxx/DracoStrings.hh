//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/DracoStrings.hh
 * \author Kelly G. Thompson <kgt@lanl.gov
 * \brief  Encapsulates common string manipulations.
 * \note   Copyright (C) 2024-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_DracoStrings_hh
#define rtt_dsxx_DracoStrings_hh

#include "Assert.hh"
#include "dsxx/config.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <locale>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Convert a string to all lower case
 *
 * \param[in] string_in This string will be converted letter by letter to lowercase.
 * \return A string that contains no uppercase letters.
 *
 * There are many complexities not considered here (e.g.: non-ASCI character sets) and many third
 * party libraries like Boost provide a more complete solution.
 */
std::string string_tolower(std::string const &string_in);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Convert a string to all upper case
 *
 * \param[in] string_in This string will be converted letter by letter to uppercase.
 * \return A string that contains no lowercase letters.
 *
 * \sa string_toupper
 */
std::string string_toupper(std::string const &string_in);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Convert a number into a string
 *
 * \param[in] num The integral or float-point type number that will be converted to a string.
 * \param[in] precision How many digits will be preserved? (default: 23)
 * \return A string representation of the numeric value
 *
 * \tparam T type of number that will be converted to a string. Typically a double or int.
 *
 * \sa http://public.research.att.com/~bs/bs_faq2.html
 */
template <typename T> std::string to_string(T const num, unsigned int const precision = 23) {
  std::ostringstream s;
  s.precision(precision);
  s << num;
  return s.str();
}

//------------------------------------------------------------------------------------------------//
//! trim whitespace (or other characters) from before and after main text.
std::string trim(std::string const &str, std::string const &whitespace = " \t\n");

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Removes all specified characters from a string.
 *
 * \param[in] orig_str The string that will be processed
 * \param[in] chars_to_remove A set of characters (as a std::string) that will be removed.
 * \return A new, possibly shortened, string that does not contain the unwanted characters.
 */
std::string prune(std::string const &orig_str, std::string const &chars_to_remove);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Split a string into a vector<string> using a specified delimiter.
 *
 * \param[in] str The string that will be split
 * \param[in] delimiters Set of character that separate each 'word'.  (default: " ")
 * \param[in] keepEmptyStrings Should empty values be saved in the result?  (default: false)
 * \return    a vector of strings. The delimiter is always removed.
 */
std::vector<std::string> tokenize(std::string const &str, std::string const &delimiters = " ",
                                  bool keepEmptyStrings = false);

//! Concept that represents simple numbers like floats or integers.
template <typename T>
concept arithmetic = std::integral<T> || std::floating_point<T>;

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Convert a string into a floating-point type or an integral type without error checking.
 *
 * \param[in] str The string that contains a number.
 * \return A numeric value.
 *
 * \tparam T The template type will be deduced based on the return type ("-> T").
 */
template <arithmetic T> auto parse_number_impl(std::string const &str) -> T;

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Convert a string into a floating-point type or an integral type with error checking.
 *
 * \param[in] str The string that contains a number.
 * \param[in] verbose Should the function print conversion message warnings (default: true).
 * \return A numeric value.
 *
 * \tparam T The template type will be deduced based on the return type ("-> T").
 *
 * These functions will throw exceptions if the conversion is invalid or if the converted value is
 * out-of-range. This generic function calls unchecked parse_number_impl specialization that is
 * appropriate.
 *
 * Consider catching thrown conversion errors using code similar to this:
 *
 * \code
 try {
    parse_number<T>(str);
  } catch (std::invalid_argument &e) {
    if no conversion could be performed
    std::cerr << "\n==ERROR==\nrtt_dsxx::parse_number:: "
              << "No valid conversion from string to a numeric value could be found.\n"
              << "\tstring = \"" << str << "\"\n" << std::endl;
    throw e;
  } catch (std::out_of_range &e) {
    if the converted value would fall out of the range of the result type
    or if the underlying function (std::strtol or std::strtoull) sets
    errno to ERANGE.
    std::cerr << "\n==ERROR==\nrtt_dsxx::parse_number:: "
              << "Type conversion from string to a numeric value resulted in a value that is out "
              << "of range.\n\tstring = \"" << str << "\"\n" << std::endl;
    throw;
  } catch (...) {
    everything else
    Insist(false, "Unknown failures in call to std::sto[x] from parse_number.");
  }
 * \endcode
 */
template <arithmetic T> auto parse_number(std::string const &str, bool verbose = true) -> T {
  T retval(0);
  try {
    retval = parse_number_impl<T>(str);
  } catch (std::invalid_argument &e) {
    // if no conversion could be performed
    if (verbose)
      std::cerr << "\n==ERROR==\nrtt_dsxx::parse_number:: "
                << "No valid conversion from string to a numeric value could be found.\n"
                << "\tstring = \"" << str << "\"\n\n";
    throw e;
  } catch (std::out_of_range & /*error*/) {
    // if the converted value would fall out of the range of the result type or if the underlying
    // function (std::strtol or std::strtoull) sets errno to ERANGE.
    if (verbose)
      std::cerr << "\n==ERROR==\nrtt_dsxx::parse_number:: Type conversion from string to a numeric"
                << " value resulted in a value that is out of range.\n"
                << "\tstring = \"" << str << "\"\n\n";
    throw;
  } catch (...) {
    // everything else
    Insist(false, "Unknown failures in call to std::sto[x] from parse_number.");
  }
  return retval;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Convert a string into a vector of floating-point or an integral values.
 *
 * \param[in] str The string that contains a number.
 * \param[in] range_symbols Parenthesis or braces that mark the beginning or end of the value
 *                 range. (default: "{}")
 * \param[in] delimiters A character that separates each numeric entry.  (default: ",")

 * \return A vector of numeric values.
 *
 * \tparam T plain-old-data type for storing numeric values, such as double, unsigned, float or int.
 *
 * std::string foo = "{ 1.0, 2.0, 3.0 }" will be converted to
 * std::vector<double> bar = { 1.0, 2.0, 3.0 }
 */
template <arithmetic T>
std::vector<T> string_to_numvec(std::string const &str, std::string const &range_symbols = "{}",
                                std::string const &delimiters = ",") {

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

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parse an ostream and create a map containing [word]:[num_occurances]
 *
 * \param[in] data An ostringstream that will be parsed
 * \param[in] verbose Echo the data map to stdout when this function is run (default: false).
 * \return a map<string,uint> that contains [word]:[num_occurances]
 */
std::map<std::string, unsigned> get_word_count(std::ostringstream const &data,
                                               bool verbose = false);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parse a file and create a map containing [word]:[num_occurances]
 *
 * \param[in] filename The name of the file (full path recommended) that will be parsed.
 * \param[in] verbose Echo the data map to stdout when this function is run (default: false).
 * \return a map<string,uint> that contains [word]:[num_occurances]
 */
std::map<std::string, unsigned> get_word_count(std::string const &filename, bool verbose = false);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Remove control characters from string that might produce color, etc.
 *
 * \param[in] colored_string A string with ANSI escape sequences that enable color output (see
 *               DracoTerminal.hh)
 * \return A string w/o ANSI escape sequences.
 *
 * As mentioned at
 * https://superuser.com/questions/380772/removing-ansi-color-codes-from-text-stream, another way to
 * strip color is to pipe output through sed to remove the color codes...
 * \code
 * sed 's/\x1b\[[0-9;]*m//g'           # Remove color sequences only
 * sed 's/\x1b\[[0-9;]*[a-zA-Z]//g'    # Remove all escape sequences
 * sed 's/\x1b\[[0-9;]*[mGKH]//g'      # Remove color and move sequences
 * sed 's/\x1b\[[0-9;]*[mGKF]//g'      # Remove color and move sequences
 * \endcode
 */
std::string remove_color(std::string const &colored_string);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Extract version of the form M.m.p from a string.
 *
 * \param[in] string_in String that contains a version number as a substring.
 * \param[in] digits number of version digits to print. Default = 3.
 * \return The version as a substring.
 *
 * Ex: Given the string "/ccs/opt/vendors-ec/ndi/2.1.3/share/gendir.all", return "2.1.3".
 *
 * A version includes any integers separated by a period.  A text postfix 'alpha' or 'beta' is also
 * supported for any portion of the version string.
 */
std::string extract_version(std::string const &string_in, size_t digits = 3);

} // namespace rtt_dsxx

#endif // rtt_dsxx_DracoStrings_hh

//------------------------------------------------------------------------------------------------//
// end of DracoStrings.hh
//------------------------------------------------------------------------------------------------//
