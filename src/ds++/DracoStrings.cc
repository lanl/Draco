//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/DracoStrings.cc
 * \author Kelly G. Thompson <kgt@lanl.gov>
 * \date   Wednesday, Aug 23, 2017, 12:48 pm
 * \brief  Encapsulates common string manipulations (implementation).
 * \note   Copyright (C) 2017-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "DracoStrings.hh"
#include "Assert.hh"
#include <fstream>
#include <iostream>
#include <regex>

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
//! Convert a string to all lower case
std::string string_tolower(std::string const &string_in) {
  std::locale loc;
  std::ostringstream string_out;
  for (auto elem : string_in)
    string_out << std::tolower(elem, loc);
  return string_out.str();
}

//------------------------------------------------------------------------------------------------//
//! Convert a string to all upper case
std::string string_toupper(std::string const &string_in) {
  std::locale loc;
  std::ostringstream string_out;
  for (auto elem : string_in)
    string_out << std::toupper(elem, loc);
  return string_out.str();
}

//------------------------------------------------------------------------------------------------//
// Definitions for fully specialized template functions
//------------------------------------------------------------------------------------------------//

template <> auto parse_number_impl<int32_t>(std::string const &str) -> int32_t {
  return std::stoi(str);
}
template <> auto parse_number_impl<int64_t>(std::string const &str) -> int64_t {
  return std::stol(str);
}
template <> auto parse_number_impl<uint32_t>(std::string const &str) -> uint32_t {
  return static_cast<uint32_t>(std::stoul(str));
}
template <> auto parse_number_impl<uint64_t>(std::string const &str) -> uint64_t {
  return std::stoull(str); // use stoull or stul?
}

// See notes in DracoStrings.hh about this CPP block
#if defined(_MSC_VER) || defined(APPLE)

template <> auto parse_number_impl<long>(std::string const &str) -> long {
  return std::stol(str); // use stoull or stul?
}
template <> auto parse_number_impl<unsigned long>(std::string const &str) -> unsigned long {
  return std::stoul(str); // use stoull or stul?
}
#endif

template <> auto parse_number_impl<float>(std::string const &str) -> float {
  return std::stof(str);
}
template <> auto parse_number_impl<double>(std::string const &str) -> double {
  return std::stod(str);
}

//------------------------------------------------------------------------------------------------//
// Definitions for regular functions
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief trim whitespace (or other characters) from before and after main text.
 *
 * \param[in] str The string that will be processed
 * \param[in] whitespace A set of characters that will be removed.  (default: " \t")
 * \return A new, probably shortened, string without unwanted leading/training characters.
 */
std::string trim(std::string const &str, std::string const &whitespace) {
  std::string retval(str.data()); // prune at first '\0' if str is a c-string.
  auto const strBegin = retval.find_first_not_of(whitespace);
  if (strBegin == std::string::npos)
    return ""; // no content
  auto const strEnd = retval.find_last_not_of(whitespace);
  auto const strRange = strEnd - strBegin + 1;
  return retval.substr(strBegin, strRange);
}

//------------------------------------------------------------------------------------------------//
//! Removes all specified characters from a string.
std::string prune(std::string const &orig_str, std::string const &chars_to_remove) {
  std::string str(orig_str);
  for (char c : chars_to_remove)
    str.erase(std::remove(str.begin(), str.end(), c), str.end());
  return str;
}

//------------------------------------------------------------------------------------------------//
//! Split a string into a vector<string> using a specified delimiter.
std::vector<std::string> tokenize(std::string const &str, std::string const &delimiters,
                                  bool keepEmptyStrings) {
  std::vector<std::string> retval; // Storage for the result
  // convert a string into a stream to be processed by getline.

  // Simple implementation if only one delimiter
  if (delimiters.size() == 1) {
    std::istringstream iss(str);
    std::string single_word; // local storage
    while (getline(iss, single_word, delimiters[0]))
      if (single_word.length() > 0 || keepEmptyStrings)
        retval.push_back(trim(single_word));
  } else {
    // Allow multiple delimiter characters.
    size_t prev = 0;
    size_t next = 0;
    while ((next = str.find_first_of(delimiters, prev)) != std::string::npos) {
      if (keepEmptyStrings || (next - prev != 0))
        retval.push_back(str.substr(prev, next - prev));
      prev = next + 1;
    }
    if (prev < str.size())
      retval.push_back(str.substr(prev));
  }

  return retval;
}

//------------------------------------------------------------------------------------------------//
//! Parse msg to provide a list of words and the number of occurrences of each.
std::map<std::string, unsigned> get_word_count(std::ostringstream const &msg, bool verbose) {
  using std::cout;
  using std::endl;
  using std::map;
  using std::string;

  map<string, unsigned> word_list;
  string msgbuf(msg.str());
  string delims(" \n\t:,.;");

  { // Build a list of words found in msgbuf.  Count the number of occurrences.

    // Find the beginning of the first word.
    string::size_type begIdx = msgbuf.find_first_not_of(delims);
    string::size_type endIdx;

    // While beginning of a word found
    while (begIdx != string::npos) {
      // search end of actual word
      endIdx = msgbuf.find_first_of(delims, begIdx);
      if (endIdx == string::npos)
        endIdx = msgbuf.length();

      // the word is we found is...
      string word(msgbuf, begIdx, endIdx - begIdx);

      // add it to the map
      word_list[word]++;

      // search to the beginning of the next word
      begIdx = msgbuf.find_first_not_of(delims, endIdx);
    }
  }

  if (verbose) {
    cout << "The messages from the message stream contained the following "
         << "words/occurrences." << endl;
    // print the word_list
    for (auto it : word_list)
      cout << it.first << ": " << it.second << endl;
  }

  return word_list;
}

//------------------------------------------------------------------------------------------------//
//! Parse text file to provide a list of words and the number of occurrences of each.
std::map<std::string, unsigned> get_word_count(std::string const &filename, bool verbose) {
  // open the file
  std::ifstream infile;
  infile.open(filename.c_str());
  Insist(infile, std::string("Cannot open specified file = \"") + filename + std::string("\"."));

  // read and store the text file contents
  std::ostringstream data;
  std::string line;
  if (infile.is_open())
    while (infile.good()) {
      getline(infile, line);
      data << line << std::endl;
    }

  infile.close();
  return get_word_count(data, verbose);
}

//------------------------------------------------------------------------------------------------//
//! Remove control characters from string that might produce color, etc.
std::string remove_color(std::string const &colored_string) {
  std::regex color_regex("\033["
                         "[[:digit:]]+[m]");
  return std::regex_replace(colored_string, color_regex, "");
}

//------------------------------------------------------------------------------------------------//
//! Extract version of the form M.m.p from a string.
std::string extract_version(std::string const &string_in, size_t digits) {
  Require(digits > 0);
  // build up a regex that will be used for version extraction.
  std::string regexString{"[A-Za-z/\\\\]*([[:digit:]]+[abehlpt]*)"};
  for (size_t i = 1; i < digits; ++i)
    regexString += "[.]([[:digit:]]+[abehlpt]*)";
  regexString += "(.*)";
  // instantiate the regex object and create variable for results.
  std::regex VerNumRegex{regexString.c_str()};
  std::smatch matchResults;
  // Search
  std::regex_match(string_in, matchResults, VerNumRegex);
  // Build the resulting version string
  std::string version;
  if (matchResults.size() > 0) {
    version += matchResults.str(1);
    for (size_t i = 2; i <= digits; ++i)
      version += "." + matchResults.str(i);
  }
  return version;
}

} // namespace rtt_dsxx

//------------------------------------------------------------------------------------------------//
// end of DracoStrings.cc
//------------------------------------------------------------------------------------------------//
