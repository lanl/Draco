//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/XGetopt.hh
 * \author Kelly Thompson, Katherine Wang
 * \date   Tuesday, Oct 27, 2016, 15:17 pm
 * \brief  Command line argument handling similar to getopt.
 * \note   Copyright (C) 2015-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_XGetopt_hh
#define rtt_dsxx_XGetopt_hh

#include "Assert.hh"
#include "ds++/config.h"
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace rtt_dsxx {

//================================================================================================//
/*!
 * \class XGetopt
 *
 * This class provides features related to parsing command line options.
 *
 * In 'main', the variables argc and argv are the argument count and array of strings, respectively,
 * as passed into the application on program invocation.  XGetopt provides tools to store the
 * command line, registered options, processing of options that read a user provided value, and
 * provides tools for presenting a help string and information about unprocessed data.
 *
 * On construction, XGetopt records registered options and help strings.  These are matched against
 * provided command line arguments and accessor functions are provided for retrieving this
 * information.
 *
 * \todo Option letters may be combined, e.g., "-ab" is equivalent to "-a -b".  Option letters are
 * case sensitive. (This features is not yet available as of 2016 October).
 *
 * If the special option "--" is encountered, all futher processing of command line arguments
 * ceases.
 *
 * Sample use:
 * \code
 *       rtt_dsxx::XGetopt::csmap long_options;
 *       long_options['b'] = "append";
 *       long_options['c'] = "create:";      // expects value.
 *       long_options['h'] = "help";
 *
 *       std::map<char,std::string> help_strings;
 *       help_strings['b'] = "perform action 'append'.";
 *       help_strings['c'] = "create a new file with provided\nfilename.";
 *       help_strings['h'] = "print this message.";
 *
 *       rtt_dsxx::XGetopt program_options( argc, argv, long_options,
 *                                          help_strings );
 *
 *       int c(0);
 *       while( (c = program_options()) != -1 )
 *       {
 *          switch (c)
 *          {
 *            case 'b':
 *               //
 *               // set some flag here
 *               //
 *               break;
 *
 *            case 'c':
 *               cvalue = program_options.get_option_value();
 *               break;
 *
 *            case 'h':
 *               std::cout << program_options.display_help( "myApplication" )
 *                         << std::endl;
 *               return 0;
 *               break;
 *
 *            default:
 *               return 0; // nothing to do.
 *          }
 *       }
 * \endcode
 *
 * \todo Future: Consider replacing this class with Boost::program_options.
 * \todo Future: Consider using a constructor that takes a
 *       map<char,pair<string,string>> in place of longopts+helpstrings.
 *
 * \example ds++/test/tstXGetopt.cc
 */
//================================================================================================//
class XGetopt {
public:
  // typedefs
  using csmap = std::map<char, std::string>;

  /*!
   * \brief Xgetopt constructor that only supports single character command line arguments.
   *
   * \param[in] argc The number of command line arguments
   * \param[in] argv An array of C-style strings that represent the command line.
   * \param[in] shortopts A string of letters that represent known command line arguments.  Any
   *            letter followed by a colon represents an option that requires that an argument be
   *            provided.
   * \param[in] helpstrings_ A map of strings whose keys must match values proivded in shortopts and
   *            whose values are messages related to each option that will be printed by the member
   *            function display_help.
   */
  XGetopt(int const argc, char **&argv, std::string const &shortopts, csmap helpstrings_ = csmap())
      : optind(0), optarg(),                                     // empty string
        cmd_line_args(argv + 1, argv + argc), longopts(csmap()), // empty
        vshortopts_hasarg(std::map<char, bool>()), vshortopts(decompose_shortopts(shortopts)),
        helpstrings(std::move(helpstrings_)) {
    match_args_to_options();
  }

  // const argv version used by the unit test.
  XGetopt(int const argc, char const *const *argv, std::string const &shortopts,
          csmap helpstrings_ = csmap())
      : optind(0), optarg(),                                     // empty string
        cmd_line_args(argv + 1, argv + argc), longopts(csmap()), // empty
        vshortopts_hasarg(std::map<char, bool>()), vshortopts(decompose_shortopts(shortopts)),
        helpstrings(std::move(helpstrings_)) {
    match_args_to_options();
  }

  /*!
   * \brief Xgetopt constructor that requires a map that includes char-to-string names of options
   *        that are to be recognized.
   *
   * \param[in] argc The number of command line arguments
   * \param[in] argv An array of C-style strings that represent the command line.
   * \param[in] longopts_ A map whose keys are letters that represent known command line arguments.
   *            Associated values are strings that provide identical functionality.  These string
   *            values may be followed by a colon represents an option that requires that an
   *            additional argument be provided.
   * \param[in] helpstrings_ A map of strings whose keys must match values proivded in shortopts and
   *            whose values are messages related to each option that will be printed by the member
   *            function display_help.
   */
  XGetopt(int const argc, char **&argv, csmap const &longopts_, csmap helpstrings_ = csmap())
      : optind(0), optarg(), // empty string
        cmd_line_args(argv + 1, argv + argc), longopts(store_longopts(longopts_)),
        vshortopts_hasarg(std::map<char, bool>()), vshortopts(decompose_longopts(longopts_)),
        helpstrings(std::move(helpstrings_)) {
    match_args_to_options();
  }

  // const argv version used by the unit test
  XGetopt(int const argc, char const *const *argv, csmap const &longopts_,
          csmap helpstrings_ = csmap())
      : optind(0), optarg(), // empty string
        cmd_line_args(argv + 1, argv + argc), longopts(store_longopts(longopts_)),
        vshortopts_hasarg(std::map<char, bool>()), vshortopts(decompose_longopts(longopts_)),
        helpstrings(std::move(helpstrings_)) {
    match_args_to_options();
  }

  /*! \brief This operator should be called from within a loop that iterates through command line
   * arguments (see class documentation).  It will return the char value that matches registered
   * options. */
  int operator()() {
    if (optind >= matched_arguments.size())
      return -1;
    if (matched_arguments_values.count(matched_arguments[optind]) > 0)
      optarg = matched_arguments_values[matched_arguments[optind]];
    else
      optarg.clear();
    return matched_arguments[++optind - 1];
  }

  /*! \brief fetch the value associated with current option. E.g.: for the option '-v foobar',
   *         return string 'foobar' associated with option 'v'. */
  std::string get_option_value() const { return optarg; }

  /*! \brief Return a vector of strings that were provided on the command line but do not match any
   *         registered arguments. */
  std::vector<std::string> get_unmatched_arguments() const { return unmatched_arguments; }

  //! reset the global counters
  // void reset(void) { optind=0; optarg=std::string(""); return; }

  //! Print a help/usage message.
  std::string display_help(std::string const &appName) const;

private:
  // >>> DATA

  //! Index to the currently processed command line argument.
  size_t optind;

  /*! \brief String provided as value for option. E.g.: for the option '-v foobar', return string
             'foobar' associated with option 'v'. */
  std::string optarg;

  //! storage for command line arguments (copy of argv)
  std::vector<std::string> const cmd_line_args;

  //! A list chars that have equivalent string names ( -v | --version ).
  csmap const longopts;

  /*! \brief Vector of flags that specify if the current option requires that an value be provided
   *         on the command line. */
  std::map<char, bool> vshortopts_hasarg;

  //! A vector of characters that represent known options.
  std::vector<char> const vshortopts;

  //! A list of matched command line arguments.
  std::vector<char> matched_arguments;

  //! A list of values associated with matched command line arguments
  std::map<char, std::string> matched_arguments_values;

  //! A list of unmatched command line arguments
  std::vector<std::string> unmatched_arguments;

  //! help strings used by display_help.
  csmap const helpstrings;

  // >>> IMPLEMENTATION

  //! Convert a string into easy-to-parse vectors and strip colons.
  std::vector<char> decompose_shortopts(std::string const &shortopts);

  //! Convert a map<char,string> into easy-to-parse vectors and strip colons.
  std::vector<char> decompose_longopts(csmap const &longopts);

  //! Cleanup and save the user-provided longoptions map. Colons are stripped.
  csmap store_longopts(csmap const &longopts_);

  /*! \brief Match provided command line arguments to registered options and
   *         record the results to class data. */
  void match_args_to_options();
};

} // end namespace rtt_dsxx

#endif // rtt_dsxx_XGetopt_hh

//------------------------------------------------------------------------------------------------//
// end of XGetOpt.hh
//------------------------------------------------------------------------------------------------//
