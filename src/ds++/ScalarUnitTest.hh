//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/ScalarUnitTest.hh
 * \author Kelly Thompson
 * \date   Thu May 18 17:08:54 2006
 * \brief  Provide services for scalar unit tests
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef dsxx_ScalarUnitTest_hh
#define dsxx_ScalarUnitTest_hh

#include "UnitTest.hh"

namespace rtt_dsxx {

//===========================================================================//
/*!
 * \class ScalarUnitTest
 * \brief This class provides services for scalar unit tests.
 *
 * This class inherits from UnitTest.  Much of the documentation for the
 * services of this class is provided in UnitTest.hh
 *
 * \sa rtt_dsxx::UnitTest for detailed description of all the UnitTest
 * classes.
 *
 * \par Code Sample:
 *
 * Scalar UnitTests should have the following syntax.
 * \code
 int main(int argc, char *argv[])
 {
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
/*!
 * \example ds++/test/tstScalarUnitTest.cc
 * The unit test for and example usage of the ScalarUnitTest class.
 */
//===========================================================================//

class ScalarUnitTest : public UnitTest {
public:
  // CREATORS

  //! Default constructors.
  DLL_PUBLIC_dsxx ScalarUnitTest(int &argc, char **&argv,
                                 string_fp_void release_,
                                 std::ostream &out_ = std::cout,
                                 bool verbose_ = true);

  //! The copy constructor is disabled.
  ScalarUnitTest(const ScalarUnitTest &rhs);

  //! Destructor.
  ~ScalarUnitTest(void) {
    out << resultMessage() << std::endl;
    return;
  };

  // MANIPULATORS

  //! The assignment operator for ScalarUnitTest is disabled.
  ScalarUnitTest &operator=(const ScalarUnitTest &rhs);
};

} // end namespace rtt_dsxx

#endif // dsxx_ScalarUnitTest_hh

//---------------------------------------------------------------------------//
// end of ds++/ScalarUnitTest.hh
//---------------------------------------------------------------------------//
