//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/UnitTest.cc
 * \author Kelly Thompson
 * \date   Thu May 18 15:46:19 2006
 * \brief  Implementation file for UnitTest.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "UnitTest.hh"
#include "path.hh"
#include <fstream>
#include <sstream>

#ifdef DRACO_DIAGNOSTICS_LEVEL_3
#include "fpe_trap.hh"
#endif

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for UnitTest object.
 *
 * \param argv A list of command line arguments.
 * \param release_ A function pointer to the local package's release() function.
 * \param out_ A user selectable output stream.  By default this is std::cout.
 * \param verbose_ Print the messages for passing tests. By default, this is set to true.
 *
 * This constructor automatically parses the command line to setup the name of the unit test (used
 * when generating status reports).  The object produced by this constructor will respond to the
 * command line argument "--version."
 */
UnitTest::UnitTest(int & /* argc */, char **&argv, string_fp_void release_, std::ostream &out_,
                   bool const verbose_)
    : numPasses(0), numFails(0), fpe_trap_active(false),
      testName(getFilenameComponent(std::string(argv[0]), rtt_dsxx::FC_NAME)),
      testPath(getFilenameComponent(
          getFilenameComponent(std::string(argv[0]), rtt_dsxx::FC_REALPATH), rtt_dsxx::FC_PATH)),
      release(release_), out(out_), m_dbcRequire(false), m_dbcCheck(false), m_dbcEnsure(false),
      m_dbcNothrow(false), verbose(verbose_) {
  Require(release != nullptr);
  Ensure(numPasses == 0);
  Ensure(numFails == 0);
  Ensure(testName.length() > 0);
#if DBC & 1
  m_dbcRequire = true;
#endif
#if DBC & 2
  m_dbcCheck = true;
#endif
#if DBC & 4
  m_dbcEnsure = true;
#endif
#if DBC & 8
  m_dbcNothrow = true;
#endif

// Turn on fpe_traps at level 3.
#ifdef DRACO_DIAGNOSTICS_LEVEL_3
  // if set to false, most compilers will provide a stack trace.
  // if set to true, fpe_trap forms a simple message and calls Insist.
  bool const abortWithInsist(true);
  rtt_dsxx::fpe_trap fpeTrap(abortWithInsist);
  fpe_trap_active = fpeTrap.enable();
#endif

  return;
}

//------------------------------------------------------------------------------------------------//
//! Build the final message that will be desplayed when UnitTest is destroyed.
std::string UnitTest::resultMessage() const {
  using DT = Term::DracoTerminal;
  std::ostringstream msg;
  std::string const div("*********************************************");
  if (UnitTest::numPasses > 0 && UnitTest::numFails == 0) {
    msg << "\n"
        << Term::ccolor(DT::pass) << div << "\n"
        << "**** " << testName << " Test: PASSED.\n"
        << div << Term::ccolor(DT::reset) << "\n";
  } else {
    msg << "\n"
        << Term::ccolor(DT::fail) << div << "\n"
        << "**** " << testName << " Test: FAILED.\n"
        << div << Term::ccolor(DT::reset) << "\n";
  }
  return msg.str();
}

//------------------------------------------------------------------------------------------------//
/*!\brief Increment the failure count and print a message with the source line number.
 * \param[in] line The line number of the source code where the failure was ecnountered.
 */
bool UnitTest::failure(int line) {
  using DT = Term::DracoTerminal;
  out << Term::ccolor(DT::fail) << "Test: failed on line " << line << Term::ccolor(DT::reset)
      << std::endl;
  UnitTest::numFails++;
  return false;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Increment the failure count and print a message with the source line number.
 * \param line The line number of the source code where fail was called from.
 * \param file The name of the file where the failure occured.
 */
bool UnitTest::failure(int line, char const *file) {
  using DT = Term::DracoTerminal;
  out << Term::ccolor(DT::fail) << "Test: failed on line " << line << " in " << file
      << Term::ccolor(DT::reset) << std::endl;
  UnitTest::numFails++;
  return false;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Increment the pass count and print a message that a test passed.
 * \param passmsg The message to be printed to the iostream \c UnitTest::out.
 */
bool UnitTest::passes(const std::string &passmsg) {
  if (verbose) {
    out << "Test: passed" << std::endl;
    out << "     " << passmsg << std::endl;
  }
  UnitTest::numPasses++;
  return true;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Increment either the pass or fail count and print a test description.
 *
 * This function is intended to reduce the number of uncovered branches in the test suite when
 * performing coverage analysis. It is used as follows:
 *
 * \code
 *     ut.check(a>0.0, "a>0.0");
 * \endcode
 *
 * If a is in fact greater than zero, the pass count is incremented and the output stream receives
 *
 * \code
 *     Test: passed
 *     a>0.0
 * \endcode
 *
 * Otherwise the fail count in incremented and the output stream receives
 *
 * \code
 *     Test: failed
 *     a>0.0
 * \endcode
 *
 * No branch is visible in the calling code that will be left uncovered during coverage analysis.
 *
 * To further reduce visible branches, a failed test optionally throws an exception, so that a
 * series of tests will be terminated if it is impossible to recover. For example, if an object
 * needed for subsequent tests is not successfully created, the test for successful creation should
 * set the fatal argument so that the sequence of tests is aborted.
 *
 * \param good True if the test passed, false otherwise.
 * \param passmsg The message to be printed to the iostream \c UnitTest::out.
 * \param fatal True if the test is to throw a std::assert on failure.
 */
bool UnitTest::check(bool const good, std::string const &passmsg, bool const fatal) {
  if (good) {
    passes(passmsg);
  } else {
    failure(passmsg);
    if (fatal)
      throw assertion("assertion thrown on critical check failure");
  }
  return true;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Increment the failure count and print a message that a test failed.
 * \param failmsg The message to be printed to the iostream \c UnitTest::out.
 */
bool UnitTest::failure(const std::string &failmsg) {
  using DT = Term::DracoTerminal;
  out << Term::ccolor(DT::fail) << "Test: failed\n"
      << "     " << failmsg << Term::ccolor(DT::reset) << std::endl;
  UnitTest::numFails++;
  return false;
}

} // end namespace rtt_dsxx

//------------------------------------------------------------------------------------------------//
// end of UnitTest.cc
//------------------------------------------------------------------------------------------------//
