//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ode/test/tstFunction_Traits.cc
 * \author Kent Budge
 * \date   Wed Aug 18 10:28:16 2004
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ode/Function_Traits.hh"
#include <iostream>
#include <typeinfo>

using namespace std;
using namespace rtt_ode;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// Helper class
//------------------------------------------------------------------------------------------------//

class Test_Functor {
public:
  using return_type = double;
};

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    UT_MSG(typeid(Function_Traits<double (*)(double)>::return_type) == typeid(double),
           "return_type correct");
    UT_MSG(typeid(Function_Traits<Test_Functor>::return_type) == typeid(Test_Functor::return_type),
           "return_type NOT correct");
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstFunction_Traits.cc
//------------------------------------------------------------------------------------------------//
