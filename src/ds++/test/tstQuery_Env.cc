//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstQuery_Env.cc
 * \author Kent G. Budge
 * \date   Wed Nov 10 09:35:09 2010
 * \brief  Test functions defined in ds++/draco_math.hh.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Query_Env.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstgetpath(rtt_dsxx::UnitTest &ut) {

  bool def_path{false};
  std::string path;

  std::tie(def_path, path) = rtt_dsxx::get_env_val<std::string>("PATH", path);

  if (def_path && path.size() > 0)
    PASSMSG("PATH was set in the environment.");
  else
    FAILMSG("Failed to read the PATH environment variable.");

  return;
}

//------------------------------------------------------------------------------------------------//
void tstgetfoobar(rtt_dsxx::UnitTest &ut) {

  bool def_foobar{false};
  std::string foobar;

  std::tie(def_foobar, foobar) = rtt_dsxx::get_env_val<std::string>("FOOBAR", foobar);

  FAIL_IF(def_foobar);
  FAIL_IF(foobar.size() > 1);

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tstgetpath(ut);
    tstgetfoobar(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstQuery_Env.cc
//------------------------------------------------------------------------------------------------//
