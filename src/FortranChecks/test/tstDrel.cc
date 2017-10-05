//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   FortranCheck/test/tstDrel.cc
 * \author Kelly Thompson
 * \date   Tuesday, Jun 12, 2012, 16:03 pm
 * \brief  Test C++ main linking a Fortran library.
 * \note   Copyright (c) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// forward declaration of f90 functions
extern "C" void drelf90(int &nf);

//---------------------------------------------------------------------------//
int main(int /*argc*/, char * /*argv*/ []) {
  int nf(0); // number of fails returned by the Fortran subroutine.  This
             // will be used as the C++ return code to trigger ctest failure
             // when nf /= 0.
  drelf90(nf);
  return nf;
}

//---------------------------------------------------------------------------//
// end of cppmain.cc
//---------------------------------------------------------------------------//
