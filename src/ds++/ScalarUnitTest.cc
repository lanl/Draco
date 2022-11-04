//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/ScalarUnitTest.cc
 * \author Kelly Thompson
 * \date   Thu May 18 17:08:54 2006
 * \brief  Provide services for scalar unit tests.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ScalarUnitTest.hh"
#include "Assert.hh"
#include "XGetopt.hh"
#include <iostream>
#include <sstream>

namespace rtt_dsxx {
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for ScalarUnitTest
 * \arg argc The number of command line arguments
 * \arg argv A list of strings containing the command line arguments
 * \arg release_ A function pointer to this package's release function.
 * \arg out_ A user specified iostream that defaults to std::cout.
 * \arg verbose_ flags whether to print messages for successful tests. Defaults to true.
 * \exception rtt_dsxx::assertion An exception with the message "Success" will be thrown if \c
 *      --version is found in the argument list.
 *
 * The constructor initializes the base class UnitTest by setting numPasses and numFails to zero.
 * It also prints a message that declares this to be a scalar unit test and provides the unit test
 * name.
 */
ScalarUnitTest::ScalarUnitTest(int &argc, char **&argv, string_fp_void release_, std::ostream &out_,
                               bool const verbose_)
    : UnitTest(argc, argv, release_, out_, verbose_) {
  using std::endl;
  using std::string;

  Require(argc > 0);
  Require(release != nullptr);

  // header
  out << "\n============================================="
      << "\n=== Scalar Unit Test: " << testName
      << "\n=============================================\n"
      << endl;

  // version tag
  out << testName << ": version " << release() << "\n" << endl;

  // Handle arguments:
  rtt_dsxx::XGetopt::csmap long_options;
  long_options['h'] = "help";
  long_options['v'] = "version";
  std::map<char, std::string> help_strings;
  help_strings['h'] = "print this message.";
  help_strings['v'] = "print version information and exit.";
  rtt_dsxx::XGetopt program_options(argc, argv, long_options, help_strings);

  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'v': // --version
      throw rtt_dsxx::assertion(string("Success"));

    case 'h': // --help
      std::cout << program_options.display_help("tstXGetopt") << std::endl;
      return;

    default:
      break; // nothing to do.
    }
  }

  Ensure(numPasses == 0);
  Ensure(numFails == 0);
  Ensure(testName.length() > 0);

  return;
}

} // end namespace rtt_dsxx

//------------------------------------------------------------------------------------------------//
// end of ScalarUnitTest.cc
//------------------------------------------------------------------------------------------------//
