//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/do_exception.cc
 * \author Rob Lowrie
 * \date   Thu Oct 13 14:33:59 2005
 * \brief  Does a floating-point exception.
 * \note   Copyright (C) 2013-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Assert.hh"
#include "ds++/StackTrace.hh"
#include "ds++/fpe_trap.hh"
#include <cfenv>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

#if defined(MSVC)
#pragma warning(push)
// This test deliberately divides by zero to test trapping of IEEE errors. Silence the warning from
// MSVC about this issue.
#pragma warning(disable : 4723)
#endif

using namespace std;

//------------------------------------------------------------------------------------------------//
/*
 * Usage: 'do_exception <integer_value>'
 *
 * <integer_value> must be 0, 1, 2 or 3.
 *
 * If the required argument is 0, then simple floating point operations are done which should not
 * cause any errors.
 *
 * If the provided integral value is 1,2, or 3, then one of the test case listed below will be run.
 * Each of these test cases should cause an IEEE floating point exception.
 *
 * Test
 * case    Description (type of IEEE exception)
 * ----    --------------------------------------
 *    1    test double division by zero
 *    2    test sqrt(-1)
 *    3    test overflow
 *
 * If \c 'DRACO_DIAGNOSTICS && 0100' is non-zero (e.g.: \c DRACO_DIAGNOSTICS=7), then Draco's
 * fpe_trap should convert the IEEE FPE into a C++ exception and print a stack trace.
 *
 * When run through via \c ctest, the output from these tests is captured in the files \c
 * do_exception_[0-9]+.out and \c do_exception_[0-9]+.err.
 */
void run_test(int /*argc*/, char **argv) {

  bool const abortWithInsist(true);
  rtt_dsxx::fpe_trap fpet(abortWithInsist);
  if (fpet.enable()) {
    // Platform is supported.
    cout << "\n- fpe_trap: This platform is supported.\n";
    if (!fpet.active())
      cout << "- fpe_trap: active flag set to false was not expected.\n";
  } else {
    // Platform is not supported.
    cout << "\n- fpe_trap: This platform is not supported\n";
    if (fpet.active())
      cout << "- fpe_trap: active flag set to true was not expected.\n";
    return;
  }

  // Accept a command line argument with value 0, 1, 2 or 3.
  int test = atoi(argv[1]);
  Insist(test >= 0 && test <= 3, "Bad test value.");

  double zero(0.0); // for double division by zero
  double neg(-1.0); // for sqrt(-1.0)
  double result(-1.0);

  // Certain tests may be optimized away by the compiler, by recognizing the constants set above and
  // precomputing the results below.  So do something here to hopefully avoid this.  This tricks the
  // optimizer, at least for gnu and KCC. [2018-02-09 KT -- some web posts hint that marking zero
  // and neg as 'volatile' may also prevent removal of logic due to optimization.]

  if (test < 0) { // this should never happen
    Insist(0, "Something is very wrong.");
    zero = neg = 1.0; // trick the optimizer?
  }

  std::ostringstream msg;
  switch (test) {
  case 0: {
    // Should not throw a IEEE floating point exception.
    cout << "- Case zero: this operation should not throw a SIGFPE. The result should be 2..."
         << endl;
    result = 1.0 + zero + sqrt(-neg);
    cout << "  result = " << result << endl;
    break;
  }
  case 1: {
    cout << "- Case one: trying a div_by_zero operation..." << endl;
    result = 1.0 / zero; // should fail here
    cout << "  result = " << 1.0 * result;
    break;
  }
  case 2: {
    // http://en.cppreference.com/w/cpp/numeric/math/sqrt
    cout << "\n- Case two: trying to evaluate sqrt(-1.0)..." << endl;
    result = std::sqrt(neg); // should fail here
    cout << "  result = " << result;
    break;
  }
  case 3: {
    cout << "- Case three: trying to cause an overflow condition..." << endl;
    result = 2.0;
    std::vector<double> data;
    for (size_t i = 0; i < 100; i++) {
      // should fail at some i
      result = result * result * result * result * result;
      data.push_back(result); // force optimizer to evaluate the above line.
    }
    cout << "  result = " << result << endl;
    break;
  }
  default: {
    std::cout << "We should never get here." << std::endl;
    break;
  }
  }
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  Insist(argc == 2, "Wrong number of args.");

  try {
    run_test(argc, argv);
  } catch (exception &err) {
    if (rtt_dsxx::fpe_trap().enable()) {
      cout << "While running " << argv[0]
           << ", a SIGFPE was successfully caught.\n\twhat = " << err.what() << endl;
      return 0;
    } else {
      cout << "ERROR: While running " << argv[0] << ", An exception was caught when it was not "
           << "expected.\n\twhat = " << err.what() << endl;
    }
  } catch (...) {
    cout << "ERROR: While testing " << argv[0] << ", An unknown exception was thrown." << endl;
    return 1;
  }
  return 0;
}

//------------------------------------------------------------------------------------------------//
#if defined(MSVC)
#pragma warning(pop)
#endif

//------------------------------------------------------------------------------------------------//
// end of do_exception.cc
//------------------------------------------------------------------------------------------------//
