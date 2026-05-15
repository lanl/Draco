//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/test/tstAssert.cc
 * \author Thomas M. Evans
 * \brief  Assertion tests.
 * \note   Copyright (C) 2024-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Assert.hh"
#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/StackTrace.hh"
#include "dsxx/UnitTest.hh"
#include <iostream>
#include <new>
#include <regex>
#include <stdexcept>
#include <string>

namespace {

using namespace std;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
// The way this test article works is that each of the DBC macros are tested in a separate function.
// A false condition is asserted using each macro, and after this follows a throw.  Two catch
// clauses are available, one to catch an assertion object, and one to catch anything else.  By
// comparing the exception that is actually caught with the one that should be caught given the DBC
// setting in force, we can determine whether each test passes or fails.
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
// Make sure we can differentiate between a std::runtime_error and a rtt_dsxx::assertion.
//------------------------------------------------------------------------------------------------//

void t1(rtt_dsxx::UnitTest &ut) {
  std::cout << "t1 test: ";
  try {
    throw std::runtime_error("hello1");
  } catch (rtt_dsxx::assertion const & /*error*/) {
    FAILMSG("rtt_dsxx::assertion caught.");
  } catch (...) {
    PASSMSG("runtime_error exception caught");
  }
  return;
}

//------------------------------------------------------------------------------------------------//
// Make sure we can catch a rtt_dsxx::assertion and extract the error message.
//------------------------------------------------------------------------------------------------//

void t2(rtt_dsxx::UnitTest &ut) {
  std::cout << "t2 test: ";
  std::string error_message;
  try {
    throw rtt_dsxx::assertion("hello1", "myfile", 42);
  } catch (rtt_dsxx::assertion const &a) {
    PASSMSG("caught rtt_dsxx::assertion");
    error_message = std::string(a.what());
  } catch (...) {
    FAILMSG("failed to catch rtt_dsxx:assertion");
  }

  // Make sure we can extract the error message.
  std::string const compare_value("Assertion: hello1, failed in myfile, line 42.");
  std::regex const rgx(std::string(".*") + compare_value + ".*");
  std::smatch match;

  if (!std::regex_search(error_message, match, rgx))
    FAILMSG(std::string("compare_value = \"") + compare_value + "\"\nmatch = \"" + match[1].str() +
            "\n");
  return;
}

//------------------------------------------------------------------------------------------------//
// Test throwing and catching of a literal
//------------------------------------------------------------------------------------------------//

void t3(rtt_dsxx::UnitTest &ut) {
  std::cout << "t3 test: ";
  try {
    throw "hello"; // NOLINT [hicpp-exception-baseclass]
  } catch (rtt_dsxx::assertion const & /*error*/) {
    FAILMSG("Should not have caught an rtt_dsxx::assertion");
  } catch (const char * /*message*/) {
    PASSMSG("Caught a const char* exception.");
  } catch (...) {
    FAILMSG("Failed to catch a const char* exception.");
  }
  return;
}

//------------------------------------------------------------------------------------------------//
// Check the toss_cookies function. This function builds an error message and throws an exception.
//------------------------------------------------------------------------------------------------//
void test_toss_cookies(rtt_dsxx::UnitTest &ut) {
  {
    std::cout << "test_toss_cookies test: ";
    try {
      std::string const msg("testing toss_cookies()");
      std::string const file("DummyFile.ext");
      int const line(55);
      rtt_dsxx::toss_cookies(msg, file, line);
      // throw "Bogus!";
    } catch (rtt_dsxx::assertion const & /* error */) {
      PASSMSG("Caught rtt_dsxx::assertion thrown by toss_cookies.");
    } catch (...) {
      FAILMSG("Failed to catch rtt_dsxx::assertion thrown by toss_cookies.");
    }
  }
  {
    std::cout << "test_toss_cookies_ptr test: ";
    try {
      char const *const msg("testing toss_cookies_ptr()");
      char const *const file("DummyFile.ext");
      int const line(56);
      rtt_dsxx::toss_cookies_ptr(msg, file, line);
      //throw "Bogus!";
    } catch (rtt_dsxx::assertion const & /* error */) {
      PASSMSG("Caught rtt_dsxx::assertion thrown by toss_cookies_ptr.");
    } catch (...) {
      FAILMSG("Failed to catch rtt_dsxx::assertion thrown by toss_cookies_ptr.");
    }
  }
  return;
}

//------------------------------------------------------------------------------------------------//
// Check the check_cookies function. This function builds an error message and throws an exception.
//------------------------------------------------------------------------------------------------//
void tcheck_cookies(rtt_dsxx::UnitTest &ut) {
  {
    std::cout << "tcheck_cookies test: ";
    try {
      rtt_dsxx::check_cookies(false, "testing check_cookies()", "DummyFile.ext", 55);
      throw "Bogus!"; // NOLINT [hicpp-exception-baseclass]
    } catch (rtt_dsxx::assertion const & /* error */) {
      PASSMSG("Caught assertion thrown by check_cookies with false condition.");
    } catch (...) {
      FAILMSG("Failed to catch assertion thrown by check_cookies with false condition.");
    }
  }
  {
    std::cout << "tcheck_cookies test: ";
    try {
      rtt_dsxx::check_cookies(true, "testing check_cookies()", "DummyFile.ext", 55);
      PASSMSG("Passed check_cookies with true condition.");
    } catch (rtt_dsxx::assertion const & /* error */) {
      PASSMSG("Bogus!");
    } catch (...) {
      FAILMSG("Failed check_cookies with true condition.");
    }
  }
  return;
}

//------------------------------------------------------------------------------------------------//
// Check the show_cookies function. This function builds an error message and throws an exception.
//------------------------------------------------------------------------------------------------//
void tshow_cookies(rtt_dsxx::UnitTest &ut) {
  using namespace std;
  {
    cout << "tshow_cookies test: \n";
    try {
      string const msg("testing show_cookies()");
      string const file("DummyFile.ext");
      int const line(55);
      cout << "The following line should be an error message...\n\t";
      rtt_dsxx::show_cookies(msg, file, line);
      throw "Bogus!"; // NOLINT [hicpp-exception-baseclass]
    } catch (rtt_dsxx::assertion const & /* error */) {
      FAILMSG("Did not expect to catch a rtt_dsxx::assertion when calling show_cookies.");
    } catch (...) {
      PASSMSG("show_cookies did not throw!");
    }
  }
  return;
}

//------------------------------------------------------------------------------------------------//
// Check the operation of the Require() macro.
//------------------------------------------------------------------------------------------------//

void trequire(rtt_dsxx::UnitTest &ut) {
  std::cout << "t-Require test: \n";
  try {
    if (ut.dbcNothrow()) {
      std::cout << "(NOTHROW=ON) The next line should be the output from Require(0) w/o an "
                << "exception thrown.\n";
    }
    Require(0);
    if (!ut.dbcNothrow())
      throw "Bogus!"; // NOLINT [hicpp-exception-baseclass]
  } catch (rtt_dsxx::assertion const &a) {
    // The nothrow option should never get here.
    FAIL_IF(ut.dbcNothrow());
    if (ut.dbcRequire()) {
      PASSMSG("trequire: caught rtt_dsxx::assertion");
      std::cout << "t-Require message value test: ";
      std::string const msg(a.what());
      std::string const expected_value("Assertion: 0, failed in");
      string::size_type const idx = msg.find(expected_value);
      FAIL_IF(idx == string::npos);
    }
    // If require is off we should never get here.
    else {
      ITFAILS;
    }
  } catch (...) {
    if (ut.dbcRequire())
      ITFAILS;
  }
  return;
}

//------------------------------------------------------------------------------------------------//
// Check the operation of the Check() macro.
//------------------------------------------------------------------------------------------------//

void tcheck(rtt_dsxx::UnitTest &ut) {
  std::cout << "t-Check test: \n";
  try {
    if (ut.dbcNothrow()) {
      std::cout << "(NOTHROW=ON) The next line should be the output from Check(false) w/o an "
                << "exception thrown.\n";
    }
    Check(false);
    if (!ut.dbcNothrow())
      throw std::runtime_error(std::string("tstAssert: tcheck()"));
  } catch (rtt_dsxx::assertion const &a) {
    // The nothrow option should never get here.
    FAIL_IF(ut.dbcNothrow());
    if (ut.dbcCheck()) {
      PASSMSG("tcheck: caught rtt_dsxx::assertion");
      std::cout << "t-Check message value test: ";
      std::string const msg(a.what());
      std::string const expected_value("Assertion: false, failed in");
      string::size_type const idx = msg.find(expected_value);
      FAIL_IF(idx == string::npos);
    }
    // If check is off we should never get here.
    else {
      ITFAILS;
    }
  } catch (...) {
    FAIL_IF(ut.dbcCheck());
  }
  return;
}

//------------------------------------------------------------------------------------------------//
// Check the operation of the Ensure() macro.
//------------------------------------------------------------------------------------------------//

void tensure(rtt_dsxx::UnitTest &ut) {
  std::cout << "t-Ensure test: \n";
  try {
    if (ut.dbcNothrow()) {
      std::cout << "(NOTHROW=ON) The next line should be the output from Ensure(0) w/o an "
                << "exception thrown.\n";
    }
    Ensure(0);
    if (!ut.dbcNothrow())
      throw "Bogus!"; // NOLINT [hicpp-exception-baseclass]
  } catch (rtt_dsxx::assertion const &a) {
    // The nothrow option should never get here.
    FAIL_IF(ut.dbcNothrow());

    if (ut.dbcEnsure()) {
      PASSMSG("tensure: caught rtt_dsxx::assertion");
      std::cout << "t-Ensure message value test: ";
      std::string const msg(a.what());
      std::string const expected_value("Assertion: 0, failed in");
      string::size_type const idx = msg.find(expected_value);
      if (idx == string::npos)
        ITFAILS;
    } else {
      ITFAILS;
    }
  } catch (...) {
    FAIL_IF(ut.dbcEnsure());
  }
  return;
}

//------------------------------------------------------------------------------------------------//
// Check the operation of the Remember() macro.
//------------------------------------------------------------------------------------------------//
void test_remember(rtt_dsxx::UnitTest &ut) {
  std::cout << "t-Remember test: ";
  int x = 0;
  Remember(x = 5);
  if (ut.dbcEnsure()) {
    FAIL_IF_NOT(x == 5);
  } else {
    FAIL_IF_NOT(x == 0);
  }
  return;
}

//------------------------------------------------------------------------------------------------//
// Check the operation of the Assert() macro, which works like Check().
//------------------------------------------------------------------------------------------------//

void tassert(rtt_dsxx::UnitTest &ut) {
  std::cout << "t-Assert test: \n";
  try {
    if (ut.dbcNothrow()) {
      std::cout << "(NOTHROW=ON) The next line should be the output from Assert(0) w/o an "
                << "exception thrown.\n";
    }
    Assert(0);
    if (!ut.dbcNothrow())
      throw "Bogus!"; // NOLINT [hicpp-exception-baseclass]
  } catch (rtt_dsxx::assertion const &a) {
    // The nothrow option should never get here.
    FAIL_IF(ut.dbcNothrow());
    if (ut.dbcCheck()) {
      PASSMSG("tassert: caught rtt_dsxx::assertion");
      std::cout << "t-Assert message value test: ";
      std::string const msg(a.what());
      std::string const expected_value("Assertion: 0, failed in");
      string::size_type const idx = msg.find(expected_value);
      if (idx == string::npos)
        ITFAILS;
    } else {
      ITFAILS;
    }
  } catch (...) {
    FAIL_IF(ut.dbcCheck());
  }
  return;
}

//------------------------------------------------------------------------------------------------//
// Basic test of the Insist() macro.
//------------------------------------------------------------------------------------------------//

void test_insist(rtt_dsxx::UnitTest &ut) {
  {
    std::cout << "t-Insist test: ";
    std::string const insist_message("You must be kidding!");
    try {
      Insist(0, insist_message);
      throw "Bogus!"; // NOLINT [hicpp-exception-baseclass]
    } catch (rtt_dsxx::assertion const &a) {
      PASSMSG("test_insist: caught rtt_dsxx::assertion");
      std::cout << "t-Insist message value test: ";
      {
        bool passed(true);
        std::string const msg(a.what());
        std::string const expected_value("You must be kidding!");
        string::size_type idx(msg.find(expected_value));
        if (idx == string::npos)
          passed = false;
        idx = msg.find(insist_message);
        if (idx == string::npos)
          passed = false;
        FAIL_IF_NOT(passed);
      }
    } catch (...) {
      ITFAILS;
    }
  }

  {
    std::cout << "t-Insist ptr test: ";
    char const *const insist_message("You must be kidding!");
    try {
      Insist_ptr(0, insist_message);
      throw "Bogus!"; // NOLINT [hicpp-exception-baseclass]
    } catch (rtt_dsxx::assertion const &a) {
      PASSMSG("test_insist_ptr: caught rtt_dsxx::assertion");
      std::cout << "t-Insist ptr message value test: ";
      {
        bool passed(true);
        std::string const msg(a.what());
        std::string const expected_value("You must be kidding!");
        string::size_type idx(msg.find(expected_value));
        if (idx == string::npos)
          passed = false;
        idx = msg.find(insist_message);
        if (idx == string::npos)
          passed = false;
        FAIL_IF_NOT(passed);
      }
    } catch (...) {
      ITFAILS;
    }
  }

  // you should always get the throwing version of Insist_device when building dsxx library tests
  {
    std::cout << "Insist device test (should map to Insist in CPU builds: ";
    char const *const insist_message("You must be kidding!");
    try {
      Insist_device(0, insist_message);
      throw "Bogus!"; // NOLINT [hicpp-exception-baseclass]
    } catch (rtt_dsxx::assertion const &a) {
      PASSMSG("test_insist_device: caught rtt_dsxx::assertion");
      std::cout << "t-Insist device message value test: ";
      {
        bool passed(true);
        std::string const msg(a.what());
        std::string const expected_value("You must be kidding!");
        string::size_type idx(msg.find(expected_value));
        if (idx == string::npos)
          passed = false;
        idx = msg.find(insist_message);
        if (idx == string::npos)
          passed = false;
        FAIL_IF_NOT(passed);
      }
    } catch (...) {
      ITFAILS;
    }
  }
  return;
}

//------------------------------------------------------------------------------------------------//
// Basic test of the Insist_ptr() macro.
//------------------------------------------------------------------------------------------------//

void test_insist_ptr(rtt_dsxx::UnitTest &ut) {
  std::cout << "t-Insist test: ";
  try {
    Insist(0, "You must be kidding!");
    throw "Bogus!"; // NOLINT [hicpp-exception-baseclass]
  } catch (rtt_dsxx::assertion const &a) {
    PASSMSG("test_insist_ptr: caught  rtt_dsxx::assertion");
    std::cout << "t-Insist_ptr message value test: ";
    {
      std::string const msg(a.what());
      std::string const expected_value("You must be kidding!");
      string::size_type const idx(msg.find(expected_value));
      FAIL_IF(idx == string::npos);
    }
  } catch (...) {
    ITFAILS;
  }
  return;
}

//------------------------------------------------------------------------------------------------//
// Check the verbose_error() function.
//------------------------------------------------------------------------------------------------//

void tverbose_error(rtt_dsxx::UnitTest &ut) {
  std::string const message(rtt_dsxx::verbose_error(std::string("This is an error.")));
  std::cout << "verbose_error() test: ";
  FAIL_IF(message.find(std::string("Host")) == std::string::npos ||
          message.find(std::string("PID")) == std::string::npos);
  return;
}

//------------------------------------------------------------------------------------------------//
// test catch of std::bad_alloc
//------------------------------------------------------------------------------------------------//
void t_catch_bad_alloc(rtt_dsxx::UnitTest &ut) {

  std::cout << "tstAssert::t_catch_bad_alloc()...\n";

  try {
    // instead of 'int * big = new int(999999999999999);'
    throw std::bad_alloc();
    //FAILMSG("failed to catch std::bad_alloc exception.");
  } catch (std::bad_alloc & /*err*/) {
    PASSMSG("caught a manually thrown std::bad_alloc exception.");
    std::cout << rtt_dsxx::print_stacktrace("Caught a std::bad_alloc") << '\n';
  } catch (...) {
    FAILMSG("failed to catch std::bad_alloc exception.");
  }

  return;
}

//------------------------------------------------------------------------------------------------//
int unused(int i) {
  switch (i) { // NOLINT [hicpp-multiway-paths-covered]
  case 0:
    return 0;

  default:
    Insist(false, "bad case");
    // Should not trigger a return with no value warning, because insist is flagged as a noreturn
    // function
  }
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try { // >>> UNIT TESTS

    // Test basic throw and catch functionality.
    t1(ut);
    t2(ut);
    t3(ut);

    // Test mechanics of Assert functions.
    test_toss_cookies(ut);
    tshow_cookies(ut);
    tcheck_cookies(ut);

    // Test Design-by-Contract macros.
    trequire(ut);
    tcheck(ut);
    tensure(ut);
    test_remember(ut);
    tassert(ut);
    test_insist(ut);
    test_insist_ptr(ut);

    // fancy output
    tverbose_error(ut);

    // noreturn - called only to keep code coverage good
    unused(0);

    // catch bad_alloc
    t_catch_bad_alloc(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstAssert.cc
//------------------------------------------------------------------------------------------------//
