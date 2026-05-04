//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/phw.cc
 * \author Kelly Thompson
 * \brief  Parallel application used by the unit test for tstApplicationUnitTest.
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/C4_Functions.hh"
#include <iostream>

using namespace rtt_c4;

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  initialize(argc, argv);
  if (argc > 1) {
    std::cout << "Hello, world!\n";
    finalize();
    return 0;
  }

  std::cout << "Goodbye cruel world...\n";
  finalize();
  return 1;
}

//------------------------------------------------------------------------------------------------//
// end of phw.cc
//------------------------------------------------------------------------------------------------//
