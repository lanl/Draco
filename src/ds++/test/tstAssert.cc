//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstAssert.cc
 * \author Thomas M. Evans
 * \date   Wed Mar 12 12:11:22 2003
 * \brief  Assertion tests.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/StackTrace.hh"
#include <regex>

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
  std::regex rgx(std::string(".*") + compare_value + ".*");
  std::smatch match;

  if (!std::regex_search(error_message, match, rgx)) {
    ITFAILS;
    std::cout << "compare_value = \"" << compare_value << "\"\nmatch = \"" << match[1] << "\n"
              << std::endl;
  }

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
void ttoss_cookies(rtt_dsxx::UnitTest &ut) {
  {
    std::cout << "ttoss_cookies test: ";
    try {
      std::string const msg("testing toss_cookies()");
      std::string const file("DummyFile.ext");
      int const line(55);
      rtt_dsxx::toss_cookies(msg, file, line);
      // throw "Bogus!";
    } catch (rtt_dsxx::assertion const & /* error */) {
      PASSMSG("Caught rtt_dsxx::assertion thrown by toss_cookies.");
    } catch (...) {
      ITFAILS;
    }
  }
  {
    std::cout << "ttoss_cookies_ptr test: ";
    try {
      char const *const msg("testing toss_cookies_ptr()");
      char const *const file("DummyFile.ext");
      int const line(56);
      rtt_dsxx::toss_cookies_ptr(msg, file, line);
      //throw "Bogus!";
    } catch (rtt_dsxx::assertion const & /* error */) {
      PASSMSG("Caught rtt_dsxx::assertion thrown by toss_cookies_ptr.");
    } catch (...) {
      ITFAILS;
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
      ITFAILS;
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
      ITFAILS;
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
      ITFAILS;
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
                << "exception thrown." << std::endl;
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
      std::string msg(a.what());
      std::string expected_value("Assertion: 0, failed in");
      string::size_type idx = msg.find(expected_value);
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
                << "exception thrown." << std::endl;
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
      std::string msg(a.what());
      std::string expected_value("Assertion: false, failed in");
      string::size_type idx = msg.find(expected_value);
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
                << "exception thrown." << std::endl;
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
      std::string msg(a.what());
      std::string expected_value("Assertion: 0, failed in");
      string::size_type idx = msg.find(expected_value);
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
void tremember(rtt_dsxx::UnitTest &ut) {
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
                << "exception thrown." << std::endl;
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
      std::string msg(a.what());
      std::string expected_value("Assertion: 0, failed in");
      string::size_type idx = msg.find(expected_value);
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

void tinsist(rtt_dsxx::UnitTest &ut) {
  {
    std::cout << "t-Insist test: ";
    std::string insist_message("You must be kidding!");
    try {
      Insist(0, insist_message);
      throw "Bogus!"; // NOLINT [hicpp-exception-baseclass]
    } catch (rtt_dsxx::assertion const &a) {
      PASSMSG("tinsist: caught rtt_dsxx::assertion");
      std::cout << "t-Insist message value test: ";
      {
        bool passed(true);
        std::string msg(a.what());
        std::string expected_value("You must be kidding!");
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
      PASSMSG("tinsist_ptr: caught rtt_dsxx::assertion");
      std::cout << "t-Insist ptr message value test: ";
      {
        bool passed(true);
        std::string msg(a.what());
        std::string expected_value("You must be kidding!");
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

  // you should always get the throwing version of Insist_device when building ds++ library tests
  {
    std::cout << "Insist device test (should map to Insist in CPU builds: ";
    char const *const insist_message("You must be kidding!");
    try {
      Insist_device(0, insist_message);
      throw "Bogus!"; // NOLINT [hicpp-exception-baseclass]
    } catch (rtt_dsxx::assertion const &a) {
      PASSMSG("tinsist_device: caught rtt_dsxx::assertion");
      std::cout << "t-Insist device message value test: ";
      {
        bool passed(true);
        std::string msg(a.what());
        std::string expected_value("You must be kidding!");
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

void tinsist_ptr(rtt_dsxx::UnitTest &ut) {
  std::cout << "t-Insist test: ";
  try {
    Insist(0, "You must be kidding!");
    throw "Bogus!"; // NOLINT [hicpp-exception-baseclass]
  } catch (rtt_dsxx::assertion const &a) {
    PASSMSG("tinsist_ptr: caught  rtt_dsxx::assertion");
    std::cout << "t-Insist_ptr message value test: ";
    {
      std::string msg(a.what());
      std::string expected_value("You must be kidding!");
      string::size_type idx(msg.find(expected_value));
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

  std::cout << "tstAssert::t_catch_bad_alloc()..." << std::endl;

  try {
    // instead of 'int * big = new int(999999999999999);'
    throw std::bad_alloc();
    //FAILMSG("failed to catch std::bad_alloc exception.");
  } catch (std::bad_alloc & /*err*/) {
    PASSMSG("caught a manually thrown std::bad_alloc exception.");
    std::cout << rtt_dsxx::print_stacktrace("Caught a std::bad_alloc") << std::endl;
  } catch (...) {
    FAILMSG("failed to catch std::bad_alloc exception.");
  }

  return;
}

//------------------------------------------------------------------------------------------------//
bool no_exception() NOEXCEPT;
bool no_exception_c() NOEXCEPT_C(true);

void tnoexcept(rtt_dsxx::UnitTest &ut) {
#if DBC
  ut.check(!noexcept(no_exception()), "with DBC on, NOEXCEPT has no effect");
  ut.check(!noexcept(no_exception_c()), "with DBC on, NOEXCEPT_C has no effect");
#else
  ut.check(noexcept(no_exception()), "with DBC off, NOEXCEPT has effect");
  ut.check(noexcept(no_exception_c()), "with DBC off, NOEXCEPT_C has effect");
#endif
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

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try { // >>> UNIT TESTS

    // Test basic throw and catch functionality.
    t1(ut);
    t2(ut);
    t3(ut);

    // Test mechanics of Assert functions.
    ttoss_cookies(ut);
    tshow_cookies(ut);
    tcheck_cookies(ut);

    // Test Design-by-Contract macros.
    trequire(ut);
    tcheck(ut);
    tensure(ut);
    tremember(ut);
    tassert(ut);
    tinsist(ut);
    tinsist_ptr(ut);

    // fancy output
    tverbose_error(ut);

    // noexcept
    tnoexcept(ut);

    // noreturn
    // called only to keep code coverage good
    unused(0);

    // catch bad_alloc
    t_catch_bad_alloc(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstAssert.cc
//------------------------------------------------------------------------------------------------//
