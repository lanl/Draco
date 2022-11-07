//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/UnitTest.hh
 * \author Kelly Thompson
 * \date   Thu May 18 15:46:19 2006
 * \brief  Provide some common functions for unit testing within Draco
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef dsxx_UnitTest_hh
#define dsxx_UnitTest_hh

#include "path.hh"
#include <list>
#include <map>
#include <memory>
#include <vector>

namespace rtt_dsxx {

//================================================================================================//
/*!
 * \class UnitTest
 * \brief Object to encapsulate unit testing of Draco classes and functions.
 *
 * This is a virtual class.  You should use one of the following UnitTest classes in your test
 * application:
 *
 * - ScalarUnitTest      - Used for testing code that does not use parallel communication (rtt_c4).
 * - ParallelUnitTest    - Used for testing code that does use parallel communications (rtt_c4).
 * - ApplicationUnitTest - Used for testing applications that run in parallel. The unit test code
 *                         is run in scalar-mode and calls mpirun to run the specified application.
 * .
 *
 * \sa UnitTest.cc for additional details.
 *
 * \par Code Sample:
 *
 * Scalar UnitTests should have the following syntax.
 * \code

 int main(int argc, char *argv[]) {
   rtt_utils::ScalarUnitTest ut( argc, argv, release );
   try { tstOne(ut); }
   UT_EPILOG(ut);
 }

 * \endcode
 *
 * \test All of the member functions of this class are tested by
 * ds++/test/tstScalarUnitTest.cc, including the early exit caused by
 * \c --version on the command line.
 */
//================================================================================================//

class UnitTest {
public:
  // NESTED CLASSES AND TYPEDEFS

  //! Typedef for function pointer to this package's release function.
  using string_fp_void = const std::string (*)();

  // CREATORS

  //! Default constructors.
  UnitTest(int & /* argc */, char **&argv, string_fp_void release_, std::ostream &out_,
           bool const verbose_ = true);

  //! The copy constructor is disabled.
  UnitTest(UnitTest const &rhs) = delete;
  //! The move constructor is disabled.
  UnitTest(UnitTest const &&rhs) = delete;
  //! The assignment copy operator is disabled.
  UnitTest &operator=(UnitTest const &rhs) = delete;
  //! The assignment move operator is disabled.
  UnitTest &operator=(UnitTest const &&rhs) = delete;

  //! Destructor is virtual because this class will be inherited from.
  virtual ~UnitTest() = default;

  // MANIPULATORS

  /*!
   * \brief Only special cases should use these (like the unit test tstScalarUnitTest.cc).
   */
  void dbcRequire(bool b) {
    m_dbcRequire = b;
    return;
  }
  void dbcCheck(bool b) {
    m_dbcCheck = b;
    return;
  }
  void dbcEnsure(bool b) {
    m_dbcEnsure = b;
    return;
  }

  // ACCESSORS
  bool failure(int line);
  bool failure(int line, char const *file);
  bool failure(std::string const &failmsg);
  bool passes(std::string const &passmsg);
  bool check(bool good, std::string const &checkmsg, bool fatal = false);
  virtual bool check_all(bool good, std::string const &checkmsg, bool fatal = false) {
    return check(good, checkmsg, fatal);
  }

  /*!
   * \brief Provide a summary of the test status
   *
   * This pure virtual function must be provided by the inherited class.  It should provide output
   * concerning the status of UnitTest.
   */
  void status() const {
    out << resultMessage() << std::endl;
    return;
  }
  //! Reset the pass and fail counts to zero.
  void reset() {
    numPasses = 0;
    numFails = 0;
    return;
  }

  bool dbcRequire() const { return m_dbcRequire; }
  bool dbcCheck() const { return m_dbcCheck; }
  bool dbcEnsure() const { return m_dbcEnsure; }
  bool dbcNothrow() const { return m_dbcNothrow; }
  bool dbcOn() const { return m_dbcRequire || m_dbcCheck || m_dbcEnsure; }
  std::string getTestPath() const { return testPath; }
  std::string getTestName() const { return testName; }
  /*!
   * \brief Returns the path of the test binary directory (useful for locating input files).
   *
   * This function depends on the cmake build system setting the COMPILE_DEFINITIONS target
   * property. This should be done in config/component_macros.cmake.
   *
   * set_target_property( unit_test_target_name COMPILE_DEFINITIONS
   *    PROJECT_BINARY_DIR="${PROJECT_BINARY_DIR}" )
   */
  static inline std::string getTestInputPath() {
#ifdef PROJECT_BINARY_DIR
    std::string sourcePath(rtt_dsxx::getFilenameComponent(PROJECT_BINARY_DIR, rtt_dsxx::FC_NATIVE));
    // if absent, append path separator.
    if (sourcePath[sourcePath.size() - 1] != rtt_dsxx::WinDirSep &&
        sourcePath[sourcePath.size() - 1] != rtt_dsxx::UnixDirSep)
      sourcePath += rtt_dsxx::dirSep;

    return sourcePath;
#else
    // We should never get here. However, when compiling ScalarUnitTest.cc, this function must be
    // valid.  ScalarUnitTest.cc is not a unit test so PROJECT_SOURCE_DIR is not defined.
    return "unknown";
#endif
  }
  /*!
   * \brief Returns the path of the test source directory (useful for locating input files).
   *
   * This function depends on the cmake build system setting the COMPILE_DEFINITIONS target
   * property. This should be done in config/component_macros.cmake.
   *
   * set_target_property( unit_test_target_name
   *    COMPILE_DEFINITIONS PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}" )
   */
  static inline std::string getTestSourcePath() {
#ifdef PROJECT_SOURCE_DIR
    std::string sourcePath(rtt_dsxx::getFilenameComponent(PROJECT_SOURCE_DIR, rtt_dsxx::FC_NATIVE));
    // if absent, append path separator.
    if (sourcePath[sourcePath.size() - 1] != rtt_dsxx::WinDirSep &&
        sourcePath[sourcePath.size() - 1] != rtt_dsxx::UnixDirSep)
      sourcePath += rtt_dsxx::dirSep;

    return sourcePath;
#else
    // We should never get here. However, when compiling ScalarUnitTest.cc, this function must be
    // valid.  ScalarUnitTest.cc is not a unit test so PROJECT_SOURCE_DIR is not defined.
    return "unknown";
#endif
  }

  // DATA
  //! The number of passes found for this test.
  unsigned numPasses;
  //! The number of failures found for this test.
  unsigned numFails;

  //! Is fpe_traping active?
  bool fpe_trap_active;

protected:
  // IMPLEMENTATION
  std::string resultMessage() const;

  // DATA

  //! The name of this unit test.
  std::string const testName;
  //! Relative path to the unit test.
  std::string const testPath;

  //! Function pointer this package's release(void) function
  string_fp_void release;

  //! Where should output be sent (default is std::cout)
  std::ostream &out;

  /*! Save the state of DBC so that it is easily accessible from within a unit
   * test.
   */
  bool m_dbcRequire;
  bool m_dbcCheck;
  bool m_dbcEnsure;
  bool m_dbcNothrow;

  /* Report successful tests? */
  bool verbose;
};

} // end namespace rtt_dsxx

#define PASSMSG(m) ut.passes(m)
#define FAILMSG(m) ut.failure(m)
#define UT_CHECK(ut, m) ut.check(m, #m);
#define UT_MSG(c, m) ut.check(c, #m);
#define ITFAILS ut.failure(__LINE__, __FILE__)
#define FAILURE ut.failure(__LINE__, __FILE__);
#define FAIL_IF_NOT(c)                                                                             \
  if (!(c))                                                                                        \
  ITFAILS
#define FAIL_IF(c)                                                                                 \
  if ((c))                                                                                         \
  ITFAILS
#define UT_EPILOG(foo)                                                                             \
  catch (rtt_dsxx::assertion & err) {                                                              \
    std::cout << "DRACO ERROR: While testing " << foo.getTestName() << ", "                        \
              << "the following error was thrown...\n"                                             \
              << err.what() << std::endl;                                                          \
    foo.numFails++;                                                                                \
  }                                                                                                \
  catch (std::exception & err) {                                                                   \
    std::cout << "ERROR: While testing " << foo.getTestName() << ", "                              \
              << "the following error was thrown...\n"                                             \
              << err.what() << std::endl;                                                          \
    foo.numFails++;                                                                                \
  }                                                                                                \
  catch (...) {                                                                                    \
    std::cout << "ERROR: While testing " << foo.getTestName() << ", "                              \
              << "An unknown exception was thrown on processor " << std::endl;                     \
    foo.numFails++;                                                                                \
  };                                                                                               \
  return foo.numFails;
#define UT_EPILOG2(foo, bar)                                                                       \
  catch (rtt_dsxx::assertion & err) {                                                              \
    std::cout << "DRACO ERROR: While testing " << foo.getTestName() << ", "                        \
              << "the following error was thrown...\n"                                             \
              << err.what() << std::endl;                                                          \
    foo.numFails++;                                                                                \
  }                                                                                                \
  catch (std::exception & err) {                                                                   \
    std::cout << "ERROR: While testing " << foo.getTestName() << ", "                              \
              << "the following error was thrown...\n"                                             \
              << err.what() << std::endl;                                                          \
    foo.numFails++;                                                                                \
  }                                                                                                \
  catch (...) {                                                                                    \
    std::cout << "ERROR: While testing " << foo.getTestName() << ", "                              \
              << "An unknown exception was thrown on processor " << std::endl;                     \
    foo.numFails++;                                                                                \
  };                                                                                               \
  return foo.numFails + bar.numFails;

#endif // dsxx_UnitTest_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/UnitTest.hh
//------------------------------------------------------------------------------------------------//
