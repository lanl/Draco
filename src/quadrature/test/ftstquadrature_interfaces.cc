//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/test/ftstquadrature_interfaces.cc
 * \author Allan Wollaber
 * \date   Mon May 23 15:34:18 MDT 2016
 * \brief  Test the correctness of the Fortran interface to quadrature_data
 * \note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

// forward declaration of f90 functions
extern "C" void test_quadrature_interfaces(void);

//------------------------------------------------------------------------------------------------//
int main(int /*argc*/, char * /*argv*/[]) {
  test_quadrature_interfaces();
  return 0;
}

//------------------------------------------------------------------------------------------------//
// end of ftstquadrature_interfaces.cc
//------------------------------------------------------------------------------------------------//
