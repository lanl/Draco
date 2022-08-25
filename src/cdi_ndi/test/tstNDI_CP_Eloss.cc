//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/test/tstNDI_CP_Eloss.cc
 * \author Ben R. Ryan
 * \date   2020 Jun 3
 * \brief  NDI_CP_Eloss test
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi/CDI.hh"
#include "cdi_ndi/NDI_CP_Eloss.hh"
#include "ds++/Query_Env.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/SystemCall.hh"
#include "ds++/dbc.hh"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using rtt_cdi_ndi::NDI_CP_Eloss;
using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void ndi_dedx_test(rtt_dsxx::UnitTest &ut) {
  int proton_zaid = 1001;
  double proton_mass = 1.6726219e-24;
  rtt_cdi::CParticle target(proton_zaid, proton_mass);
  int alpha_zaid = 2004;
  double alpha_mass = 6.64465723e-24;
  rtt_cdi::CParticle projectile(alpha_zaid, alpha_mass);

  // Write a custom gendir file to deal with NDI-required absolute path to data
  std::string gendir_in = "gendir_tmp.dedx";
  std::string gendir_path = ut.getTestInputPath() + gendir_in;
  std::string data_path = ut.getTestSourcePath() + "dedx_he4_example";
  std::ofstream gendir_tmp_file;
  gendir_tmp_file.open(gendir_path);
  gendir_tmp_file << "dedx\n";
  gendir_tmp_file << "  z=2004.000dx  d=2020-06-03  l=rpa_cut\n";
  gendir_tmp_file << "    f=" << data_path << "  ft=asc  ln=2\n";
  gendir_tmp_file << "    o=45  ng=91  aw=4.001510E+00  awr=3.96713510838 end\n";
  gendir_tmp_file << "end\n";
  gendir_tmp_file.close();

  std::string library_in = "rpa_cut";

  NDI_CP_Eloss eloss(gendir_path, library_in, target, projectile);

  // Get eloss value for almost first (1,1,1) grid point
  {
    double energy = 1.384272;
    double density = 3.344490e-01;
    double temperature = 3.981051e-04;
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss.getEloss(temperature, density, energy),
                                     2.311354474121679232e+04, 1.e-8));
  }

  // Get eloss value for almost last (2,3,4) grid point
  {
    double energy = 1.384273e+01;
    double density = 3.344495e+03;
    double temperature = 2.511868e+01;
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss.getEloss(temperature, density, energy),
                                     1.272473147179571250e+15, 1.e-8));
  }

  // Get eloss value for a point between grid points (1.5,2.5,3.5, i.e. requiring linear
  // interpolation)
  {
    double energy = 4.377453e+00;
    double density = 3.344494e+02;
    double temperature = 3.981044e+00;
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss.getEloss(temperature, density, energy),
                                     8.794247704820802689e+09, 1.e-8));
  }

  if (ut.numFails == 0) {
    PASSMSG("NDI_CP_Eloss test passes.");
  } else {
    FAILMSG("NDI_CP_Eloss test fails.");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    std::string gendir_default;
    bool def_gendir{false};
    std::tie(def_gendir, gendir_default) = rtt_dsxx::get_env_val<std::string>("NDI_GENDIR_PATH");

    if (def_gendir && rtt_dsxx::fileExists(gendir_default)) {
      ndi_dedx_test(ut);
    } else {
      PASSMSG("==> ENV{NDI_GENDIR_PATH} not set. Some tests were not run.");
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// End tstNDI_CP_Eloss.cc
//------------------------------------------------------------------------------------------------//
