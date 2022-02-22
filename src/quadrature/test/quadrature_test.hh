//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/test/quadrature_test.hh
 * \author Kent G. Budge
 * \brief  Define class quadrature_test
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef quadrature_quadrature_test_hh
#define quadrature_quadrature_test_hh

#include "ds++/UnitTest.hh"
#include "quadrature/Quadrature.hh"
#include "quadrature/Quadrature_Interface.hh"

namespace rtt_quadrature {

using rtt_dsxx::UnitTest;
void quadrature_test(UnitTest &ut, Quadrature &quadrature, bool const cartesian_tests_only = false);
void quadrature_integration_test(UnitTest &ut, Quadrature &quadrature);

} // end namespace rtt_quadrature

#endif // quadrature_quadrature_test_hh

//------------------------------------------------------------------------------------------------//
// end of quadrature/quadrature_test.hh
//------------------------------------------------------------------------------------------------//
