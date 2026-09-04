//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/XGetopt.cc
 * \author Kelly Thompson, Katherine Wang
 * \brief  Command line argument handling similar to getopt.
 * \note   Copyright (C) 2024-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "XGetopt.hh"
#include "dsxx/Assert.hh"
#include <algorithm>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <ranges>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

// cSpell: ignore csmap, mato

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Convert a string into easy-to-parse vectors and strip colons.
 * \param[in] shortopts A string of character and optional colons. Each character represents one
 *        option that is registered. Characters followed by a colon represent options that require
 *        additional information to be read from the command line.
 * \return A vector of chars that represent all support options.
 *
 * - This function also populates the class member vshortopts_hasarg.
 * - Valid shortopts are single characters or a character followed by a colon.
 *
 * \pre shortopts is non‑empty.
 * \pre Each option character is a valid, single‑byte character (no multibyte UTF‑8).
 * \pre The caller owns `vshortopts_hasarg` (e.g. an unordered_map<char,bool>) and expects it to be
 *      populated for every option character that appears in *shortopts*.
 *
 * \post keys contains exactly `shortopts.size()` characters (including any ':' characters).
 * \post For every option character `c` in *shortopts* (including ':'), `vshortopts_hasarg[c]` is
 *       set to true iff the character immediately following `c` in *shortopts* is ':'.
 * \post The returned vector is non‑empty.
 */
std::vector<char> XGetopt::decompose_shortopts(std::string_view shortopts) {
  Require(!shortopts.empty());
  std::vector<char> keys;
  keys.reserve(shortopts.size()); // avoid reallocations

  for (std::size_t i = 0; i < shortopts.size(); ++i) {
    const char opt = shortopts[i];
    keys.push_back(opt);

    // If the next character is ':' the current option takes an argument.
    const bool has_arg = (i + 1 < shortopts.size() && shortopts[i + 1] == ':');
    vshortopts_hasarg[opt] = has_arg;
  }
  Ensure(!keys.empty() && keys.size() <= shortopts.size());
  return keys;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Convert a map<char,string> into easy-to-parse vectors and strip colons.
 * \param[in] longopts_ A map<char,string> data structure that contains all supported options that
 *        are to be registered. Each entry pair represents one option that is registered, a short
 *        option that is used with a single dash (-h) and a long option that is used with a double
 *        dash (--help). Strings followed by a colon represent options that require additional
 *        information to be read from the command line.
 * \return A vector of chars that represent all support options.
 *
 * If the string value ends in a colon, ':', strip it from the map.  This data has already been
 * processed to indicate that an argument following the option is expected.
 */
XGetopt::csmap XGetopt::store_longopts(csmap const &longopts_) {
  Require(!longopts_.empty());
  csmap retValue; // build up the data to be returned in this container.
  // Convert each entry, stripping a trailing ':' if it exists.
  std::ranges::transform(longopts_, std::inserter(retValue, retValue.end()),
                         [](auto const &kv) -> std::pair<char, std::string> {
                           const char key = kv.first;
                           std::string_view val = kv.second; // view, no copy yet
                           if (!val.empty() && val.back() == ':')
                             val.remove_suffix(1);         // drop the colon
                           return {key, std::string{val}}; // copy only once
                         });
  Ensure(retValue.size() == longopts_.size());
  return retValue;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Decompose a map of long‑option characters into a plain vector of the characters while
 *        populating the short‑option‑has‑argument table.
 *
 * \param[in] longopts_  A map whose **key** is a single‑character option identifier (`char`) and
 *          whose **value** is the option description (normally a string that may end with ':' to
 *          indicate that the option expects an argument).
 * \return  A `std::vector<char>` that contains exactly one copy of each key from `longopts_`.  The
 *          order of the elements follows the iteration order of the map (unspecified for
 *          `std::unordered_map`).
 *
 * \pre `!longopts_.empty()` – the input map must contain at least one entry.
 * \pre No duplicate keys are present in `longopts_`.  The function asserts this property; if a
 *      duplicate is found the program aborts in a debug build.
 * \post The returned vector `keys` satisfies:
 *       - `keys.size() == longopts_.size()`.
 *       - Every element of `keys` is a key from `longopts_` and each key appears exactly once.
 *       - For each `key` in `longopts_` the member `vshortopts_hasarg[key]` is set to `true` iff
 *         the associated value ends with a colon (`':'`); otherwise it is set to `false`.
 *
 * \note This function is intended for internal use by `XGetopt`.  It relies on the member container
 *       `vshortopts_hasarg` to store the “has‑argument” flag for each short option.
 */
std::vector<char> XGetopt::decompose_longopts(csmap const &longopts_) {
  Require(!longopts_.empty());
  std::vector<char> keys;
  keys.reserve(longopts_.size()); // avoid repeated allocations
  for (auto const &[key, value] : longopts_) {
    // Guard against duplicate keys.
    Insist(std::ranges::find(keys, key) == keys.end(),
           "You cannot use the same single character command line argument more than once.");
    keys.push_back(key); // Create a container of keys.
    // Record whither this option takes an argument.  A trailing ':' in the description signals "has
    // argument".
    vshortopts_hasarg[key] = !value.empty() && value.back() == ':'; // create container of hasarg.
  }
  Ensure(keys.size() == longopts_.size());
  return keys;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Store the argument value that follows an option which requires one.
 *
 * This helper is called after an option (short or long) that has been identified as needing a
 * value.  It checks that a next token exists, advances the loop index so the caller skips the
 * value, and records the association in matched_arguments_values.
 *
 * \param[in]  opt   The short‑option character that requires a value (e.g. ‘a’ for “‑a”).  The
 *                   character is also the key used in matched_arguments_values.
 * \param[in,out] i Index of the current command‑line argument in cmd_line_args.  On success this
 *                   index is incremented to point at the consumed value.
 */
void XGetopt::mato_handle_value(char opt, std::size_t &i) {
  // We are looking at the *value* token that must follow the option.
  const bool has_next = (i + 1) < cmd_line_args.size();

  // Produce a helpful error message if the value is missing.
  Insist(has_next,
         (std::string{"The command argument '-"} + opt + "' expected a value to be provided."));

  // Advance `i` so the outer loop skips the value token we just consumed.
  ++i;
  matched_arguments_values[opt] = cmd_line_args[i];
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Process a short option of the form “‑x”.
 *
 * The function validates that *x* is a recognised short option, records the match, and, if the
 * option expects an argument, forwards the work to handle_value.
 *
 * \param[in]  arg        The raw command‑line token (e.g. “‑f”).  It is guaranteed to have length 2
 *                        and to start with ‘‑’.
 * \param[in]  short_set  A set containing all allowed short‑option characters.
 * \param[in,out] i       Index of the current token in cmd_line_args. May be incremented if the
 *                        option consumes a following value.
 *
 * \note Unrecognised short options are placed into unmatched_arguments.
 */
void XGetopt::mato_process_short(std::string_view arg, std::unordered_set<char> const &short_set,
                                 std::size_t &i) {
  const char opt = arg[1]; // arg is guaranteed to be of the form "-x"

  if (!short_set.contains(opt)) {
    unmatched_arguments.emplace_back(arg); // Unknown short option → unmatched.
    return;
  }

  // Record the match.
  matched_arguments.push_back(opt);

  // If this option expects a value, pull the next command‑line token.
  if (vshortopts_hasarg[opt]) {
    mato_handle_value(opt, i);
  }
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Process a long option of the form “--name”.
 *
 * The function strips the leading “--”, looks the name up in the longopts map to obtain the
 * associated short‑option character, records the match, and delegates to handle_value if the option
 * requires an argument.
 *
 * \param[in]  arg The raw command‑line token beginning with “--”.
 * \param[in,out] i Index of the current token in cmd_line_args. May be incremented if the option
 *                consumes a following value.
 *
 * \note If the long name cannot be found in longopts, the original token is appended to
 *       unmatched_arguments.
 */
void XGetopt::mato_process_long(std::string_view arg, std::size_t &i) {
  // Strip the leading “--”.
  std::string_view const long_name = arg.substr(2);

  // Find the corresponding short key (if any).
  char short_opt = '\0';
  for (auto const &[key, name] : longopts) {
    if (name == long_name) {
      short_opt = key;
      break;
    }
  }

  if (short_opt == '\0') {
    // No registered long option → unmatched.
    unmatched_arguments.emplace_back(arg);
    return;
  }

  // Record the match.
  matched_arguments.push_back(short_opt);

  // If this option needs a value, fetch it.
  if (vshortopts_hasarg[short_opt]) {
    mato_handle_value(short_opt, i);
  }
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Match the command‑line arguments against the registered short and long options.
 *
 * The function walks the vector `cmd_line_args` (typically populated from `argv`) and
 * classifies each token into one of three categories:
 *
 * - **Matched short options** – tokens of the form “‑x”.  The character *x* is verified against
 *   `vshortopts`.  If the option requires a value (according to `vshortopts_hasarg`) the next token
 *   is taken as that value.
 *
 * - **Matched long options** – tokens that start with “--”.  The text after the two dashes is
 *   looked up in the map `longopts` (which maps a short character to its long name).  The
 *   corresponding short key is stored, and any required argument is consumed exactly as for short
 *   options.
 *
 * - **Unmatched arguments** – any token that does not fit the two categories above, or that fails
 *   the validation (unknown short option or unknown long name) is copied into
 *   `unmatched_arguments`.
 *
 * The special token “--” marks the end of option parsing; all tokens that appear after it are
 * ignored by this routine.
 *
 * \details
 * - The function updates the following members of the class:
 *   - `matched_arguments` – list of short keys that were successfully matched.
 *   - `matched_arguments_values` – map from short key to the associated value (if any).
 *   - `unmatched_arguments` – list of tokens that could not be interpreted as a known option.
 *
 * - If an option that requires a value is encountered at the end of the argument list, an
 *   `option_error` exception is thrown with a descriptive message.
 *
 * - The routine runs in O(N) where N is the number of command‑line arguments; short‑option look‑ups
 *   use an `unordered_set` for O(1) access.
 */
void XGetopt::match_args_to_options() {

  // Fast O(1) lookup for short options.
  std::unordered_set<char> const short_set(vshortopts.begin(), vshortopts.end());

  for (std::size_t i = 0; i < cmd_line_args.size(); ++i) {
    std::string_view const arg = cmd_line_args[i];

    // 1) End‑of‑options marker “--”.
    if (arg == "--")
      break; // stop processing further args

    // 2) Short option (exactly two characters, e.g. "-a").
    if (arg.size() == 2 && arg[0] == '-') {
      mato_process_short(arg, short_set, i);
      continue;
    }

    // 3) Long option (starts with “--”).
    if (arg.size() > 2 && arg.starts_with("--")) {
      mato_process_long(arg, i);
      continue;
    }

    // 4) Anything else is “unmatched”.
    unmatched_arguments.emplace_back(arg);
  }
  return;
}

namespace { // anonymous namespace for helper functions used by display_help()

//------------------------------------------------------------------------------------------------//
//! Helper: indent every line after the first with a “hanging indent”.
std::string indent_multiline(std::string_view txt, std::size_t hangingIndent) {
  if (txt.empty())
    return {};

  const std::string indent(hangingIndent, ' ');
  std::string result;
  result.reserve(txt.size() + (txt.size() / 30) * hangingIndent); // rough reservation

  std::size_t pos = 0;
  while (true) {
    // copy up to the next '\n' (or the end of string)
    std::size_t const next = txt.find('\n', pos);
    result.append(txt.substr(pos, next - pos));
    if (next == std::string_view::npos)
      break;

    // we hit a line‑break – insert it and the extra indent for the next line
    result.push_back('\n');
    result.append(indent);
    pos = next + 1; // start after the '\n'
  }
  return result;
}

} // end anonymous namespace

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return a help string for the program.
 *
 * The format is roughly:
 * ```text
 * Usage:
 *   <appName> [options]
 *
 * Options without arguments:
 *   -a | --alpha            : description
 *   -b | --beta             : description
 *
 * Options requiring arguments:
 *   -c | --gamma <value>    : description that may span several lines
 *   ...
 *   -- : Stop processing any further arguments.
 * ```
 *
 * The function never throws; it only uses the data stored in the object.
 */
std::string XGetopt::display_help(std::string_view appName) const {
  Require(!appName.empty());

  // 1. Gather every option we have to display.
  struct OptionInfo {
    char shortOpt;
    std::string_view longOpt; // empty for “short‑only” options
    bool hasArg;
    std::string_view help; // raw help text (may contain '\n')
  };

  std::vector<OptionInfo> opts;

  if (!longopts.empty()) {
    for (auto const &[s, l] : longopts) {
      bool const hasArg = vshortopts_hasarg.at(s);
      auto const itHelp = helpstrings.find(s);
      std::string_view const help =
          (itHelp != helpstrings.end()) ? std::string_view{itHelp->second} : std::string_view{};
      opts.push_back(
          {.shortOpt = s, .longOpt = std::string_view{l}, .hasArg = hasArg, .help = help});
    }
  } else {
    for (char const s : vshortopts) {
      if (s == ':')
        continue; // skip the placeholder used elsewhere
      bool const hasArg = vshortopts_hasarg.at(s);
      auto itHelp = helpstrings.find(s);
      std::string_view const help =
          (itHelp != helpstrings.end()) ? std::string_view{itHelp->second} : std::string_view{};
      opts.emplace_back(s, std::string_view{}, hasArg, help);
    }
  }

  // 2️. Split into “no‑argument” and “needs‑argument” groups.
  auto noArgRange = opts | std::views::filter([](auto const &o) noexcept { return !o.hasArg; });
  auto argRange = opts | std::views::filter([](auto const &o) noexcept { return o.hasArg; });

  // 3. Determine the column width needed to align the long‑option column.
  auto longest_long_opt = [](auto &&range) -> std::size_t {
    auto it = std::ranges::max_element(range, std::less<>{},
                                       [](auto const &o) noexcept { return o.longOpt.size(); });
    return (it == std::ranges::end(range)) ? 0 : it->longOpt.size();
  };

  const std::size_t maxLongNoArg = longest_long_opt(noArgRange);
  const std::size_t maxLongArg = longest_long_opt(argRange);

  // 4. Build the final help text.
  std::ostringstream out;
  out << "\nUsage:\n  " << appName << " [options]\n\nOptions without arguments:\n";

  // ----- options that do NOT need a value ---------------------------------
  for (auto const &o : noArgRange) {
    // column 1: short option
    out << "   -" << o.shortOpt;

    // column 2: optional long option, padded to maxLongNoArg
    if (!o.longOpt.empty()) {
      out << " | --" << std::left << std::setw(static_cast<int>(maxLongNoArg)) << o.longOpt;
    }
    out << " : " << o.help << '\n';
  }

  // ----- options that need a value ----------------------------------------
  out << "\nOptions requiring arguments:\n";
  for (auto const &o : argRange) {
    // column 1: short option
    out << "   -" << o.shortOpt;

    // column 2: long option (if any) plus the literal “<value>”, padded.
    if (!o.longOpt.empty()) {
      out << " | --" << std::left << std::setw(static_cast<int>(maxLongArg)) << o.longOpt
          << " <value>";
    } else {
      out << " <value>";
    }

    // column 3: help text with hanging indent (for multi‑line strings)
    const std::size_t hangingIndent = 21 // known fixed width portion
                                      + std::max(maxLongArg, o.longOpt.size());
    out << " : " << indent_multiline(o.help, hangingIndent) << '\n';
  }

  out << "   -- : Stop processing any further arguments.\n";
  return out.str();
}

} // end namespace rtt_dsxx

//------------------------------------------------------------------------------------------------//
// end of XGetopt.cc
//------------------------------------------------------------------------------------------------//
