//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/draco_info_main.cc
 * \author Kelly Thompson
 * \date   Wednesday, Nov 07, 2012, 18:49 pm
 * \brief  Small executable that prints the version and copyright strings.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "diagnostics/draco_info.hh"
#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "ds++/XGetopt.hh"
#include <iostream>

int main(int argc, char *argv[]) {
  using std::cout;
  using std::endl;
  try {
    bool version(false);
    bool brief(false);
    bool author(false);
    bool known_arg(false);
    bool use_doxygen_formatting(false);
    rtt_diagnostics::DracoInfo di;

    // Preparing to parse command line arguments.
    rtt_dsxx::XGetopt::csmap long_options;
    long_options['a'] = "author";
    long_options['b'] = "brief";
    long_options['d'] = "use_doxygen_formatting";
    long_options['h'] = "help";
    long_options['v'] = "version";
    std::map<char, std::string> help_strings;
    help_strings['a'] = "print the author list.";
    help_strings['b'] = "print a brief message.";
    help_strings['d'] = "Add extra formatting (doxygen) to output strings.";
    help_strings['v'] = "print version information and exit.";
    help_strings['h'] = "print this message.";
    rtt_dsxx::XGetopt program_options(argc, argv, long_options, help_strings);

    int c(0);
    while ((c = program_options()) != -1) {
      switch (c) {
      case 'a': // --author
        author = true;
        known_arg = true;
        break;

      case 'b': // --brief
        brief = true;
        known_arg = true;
        break;

      case 'd': // --use_doxygen_formatting
        use_doxygen_formatting = true;
        known_arg = true;
        break;

      case 'v': // --version
        version = true;
        known_arg = true;
        break;

      default: // same as option 'h':
        cout << program_options.display_help("draco_info") << endl;
        return 0;
      }
    }

    // If a bad argument is provided, print the help message and exit.
    if (argc > 1 && !known_arg) {
      cout << program_options.display_help("draco_info") << endl;
      return 0;
    }

    // Otherwise generate and print the requested report.
    if (version)
      cout << di.versionReport();
    else if (brief)
      cout << di.briefReport();
    else if (author)
      cout << rtt_dsxx::author_list(use_doxygen_formatting);
    else
      cout << di.fullReport();
    // -----------------------------------------------------------------------//
  } catch (std::exception &err) {
    cout << "ERROR: While running " << argv[0] << ", " << err.what() << endl;
    return 1;
  } catch (...) {
    cout << "ERROR: While running " << argv[0] << ", "
         << "An unknown C++ exception was thrown" << endl;
    return 1;
  }

  return 0;
}

//------------------------------------------------------------------------------------------------//
// end of draco_info_main.cc
//------------------------------------------------------------------------------------------------//
