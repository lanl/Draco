//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/ParallelUnitTest.i.hh
 * \author Kent Grimmett Budge
 * \date   Tue Nov  6 13:12:37 2018
 * \brief  Member definitions of class test
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef dsxx_parallel_unit_test_i_hh
#define dsxx_parallel_unit_test_i_hh

namespace rtt_c4 {

template <typename... Lambda, typename Release>
int do_parallel_unit_test(int argc, char **argv, Release release, Lambda const &...lambda) {
  rtt_c4::ParallelUnitTest ut(argc, argv, release);
  try {
    rtt_dsxx::implement_do_unit_test(ut, lambda...);
  }
  UT_EPILOG(ut);
}

} // namespace rtt_c4

#endif // dsxx_test_i_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/parallel_unit_test.i.hh
//------------------------------------------------------------------------------------------------//
