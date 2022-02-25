//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/ScalarUnitTest.i.hh
 * \author Kent Grimmett Budge
 * \date   Tue Nov  6 13:12:37 2018
 * \brief  Member definitions of class test
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef dsxx_ScalarUnitTest_i_hh
#define dsxx_ScalarUnitTest_i_hh

namespace rtt_dsxx {

template <typename Lambda> void implement_do_unit_test(UnitTest &ut, Lambda const &lambda) {
  lambda(ut);
}

template <typename First_Lambda, typename... More_Lambdas>
void implement_do_unit_test(UnitTest &ut, First_Lambda const &first_lambda,
                            More_Lambdas const &...more_lambdas) {
  first_lambda(ut);
  implement_do_unit_test(ut, more_lambdas...);
}

template <typename... Lambda, typename Release>
int do_scalar_unit_test(int argc, char **argv, Release release, Lambda const &...lambda) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    implement_do_unit_test(ut, lambda...);
  }
  UT_EPILOG(ut);
}

} // end namespace rtt_dsxx

#endif // dsxx_test_i_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/ScalarUnitTest.i.hh
//------------------------------------------------------------------------------------------------//
