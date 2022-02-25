//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/ParallelUnitTest.hh
 * \author Kelly Thompson
 * \date   Thu Jun  1 17:15:05 2006
 * \brief  Declaration file for encapsulation of Draco parallel unit tests.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 *
 * This file provides a definition for ParallelUnitTest.  The purpose of this class is to
 * encapsulate the keywords and behavior of DBS parallel unit tests. */
//------------------------------------------------------------------------------------------------//

#ifndef c4_ParallelUnitTest_hh
#define c4_ParallelUnitTest_hh

#include "C4_Functions.hh"
#include "ds++/ScalarUnitTest.hh"

namespace rtt_c4 {

//================================================================================================//
/*!
 * \class ParallelUnitTest
 * \brief This class encapsulates services for parallel unit tests.
 *
 * This class inherits from UnitTest.  Much of the documentation for the services of this class is
 * provided in UnitTest.hh
 *
 * \sa rtt_dsxx::UnitTest for detailed description of all the UnitTest classes.
 *
 * \par Code Sample:
 * \code
 int main(int argc, char *argv[])
 {
   rtt_c4::ParallelUnitTest ut(argc, argv, release);
   try
   {
     tstOne(ut);
     tstTwo(ut);
   }
   UT_EPILOG(ut);
 }
 * \endcode
 *
 * \test All of the member functions of this class are tested by ds++/test/tstScalarUnitTest.cc,
 *       including the early exit caused by \c --version on the command line.
 *
 * \warning The output from this class is closely tied to the DBS python script \c
 *       tools/regression_filter.py that is used during \c gmake \c check.  Changing the format or
 *       keyword in the output streams from this class should be coordinated with the regular
 *       expression matches found in \c tools/regression_filter.py.
 *
 * \warning The differences between ScalarUnitTest, ParallelUnitTest and ApplicationUnitTest are
 *       correlated to the DBS m4 macros \c AC_RUNTESTS and \c AC_TEST_APPLICATION.  Changes to
 *       these classes should be coordinated with changes to these DBS m4 macro command
 */
/*!
 * \example c4/test/tstParallelUnitTest.cc
 * This unit test demonstrates typical usage for ParallelUnitTest. *
 */
//================================================================================================//

class ParallelUnitTest : public rtt_dsxx::UnitTest {
public:
  // CREATORS

  //! Default constructor.
  ParallelUnitTest(int &argc, char **&argv, string_fp_void release_, std::ostream &out_ = std::cout,
                   bool verbose_ = true);

  //!  The copy/move constructors are disabled.
  ParallelUnitTest(ParallelUnitTest const &rhs) = delete;
  ParallelUnitTest(ParallelUnitTest &&rhs) = delete;

  //! Destructor.
  ~ParallelUnitTest() override;

  // MANIPULATORS

  //! The copy/move assignment operators are disabled.
  ParallelUnitTest &operator=(ParallelUnitTest const &rhs) = delete;
  ParallelUnitTest &operator=(ParallelUnitTest &&rhs) = delete;

  // ACCESSORS

  //! Provide a report of the number of unit test passes and fails.
  void status();

  bool check_all(bool good, std::string const &checkmsg, bool fatal = false) override;
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Run a parallel unit test.
 *
 * \param[in] argc Number of command line arguments
 * \param[in] argv Command line arguments
 * \param[in] release Release string
 * \param[in] lambda Lambda function defining the test.
 * \return EXIT_SUCCESS or EXIT_FAILURE as appropriate.
 */
template <typename... Lambda, typename Release>
int do_parallel_unit_test(int argc, char **argv, Release release, Lambda const &...lambda);

} // end namespace rtt_c4

#include "ParallelUnitTest.i.hh"

#endif // c4_ParallelUnitTest_hh

//------------------------------------------------------------------------------------------------//
// end of c4/ParallelUnitTest.hh
//------------------------------------------------------------------------------------------------//
