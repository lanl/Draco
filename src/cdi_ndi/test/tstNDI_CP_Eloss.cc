//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/test/tstNDI_CP_Eloss.cc
 * \author Ben R. Ryan
 * \date   2020 Jun 3
 * \brief  NDI_CP_Eloss test
 * \note   Copyright (C) 2020-2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Query_Env.hh"
#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "dsxx/path.hh"
#include "cdi/CPCommon.hh"
#include "cdi_ndi/NDI_CP_Eloss.hh"
#include <cmath>
#include <fstream>
#include <string>
#include <tuple>

// cSpell:ignore Eloss, dedx

using rtt_cdi_ndi::NDI_CP_Eloss;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void ndi_dedx_test(rtt_dsxx::UnitTest &ut) {
  int const proton_zaid = 1001;
  double const proton_mass = 1.6726219e-24;
  rtt_cdi::CParticle const target(proton_zaid, proton_mass);
  int const alpha_zaid = 2004;
  double const alpha_mass = 6.64465723e-24;
  rtt_cdi::CParticle const projectile(alpha_zaid, alpha_mass);

  // Write a custom gendir file to deal with NDI-required absolute path to data
  std::string const gendir_in = "gendir_tmp.dedx";
  std::string const gendir_path = ut.getTestInputPath() + gendir_in;
  std::string const data_path = ut.getTestSourcePath() + "dedx_he4_example";
  std::ofstream gendir_tmp_file;
  gendir_tmp_file.open(gendir_path);
  gendir_tmp_file << "dedx\n";
  gendir_tmp_file << "  z=2004.000dx  d=2020-06-03  l=rpa_cut\n";
  gendir_tmp_file << "    f=" << data_path << "  ft=asc  ln=2\n";
  gendir_tmp_file << "    o=45  ng=91  aw=4.001510E+00  awr=3.96713510838 end\nend\n";
  gendir_tmp_file.close();

  std::string const library_in = "rpa_cut";

  NDI_CP_Eloss const eloss(gendir_path, library_in, target, projectile);

  // Get eloss value for slightly above (E=1,rho=1,T=1) grid point.
  {
    double const speed = 0.0982085;
    double const density = 0.00167264;
    double const temperature = 0.0000100;
    // Exponentiate the table value (in log space) at (E=1,rho=1,T=1)
    double const base_val = exp(-4.857735E+01);
    // Calculate the unit conversion (see getEloss function):
    const double multiplier{1000. * density / proton_mass * speed};
    double const reference = multiplier * base_val;
    // KPL: using loose tolerance here so we can compare directly to the data point.
    FAIL_IF_NOT(
        rtt_dsxx::soft_equiv(eloss.getEloss(temperature, density, speed), reference, 1.e-4));
  }

  // Get eloss value for slightly above (E=2,rho=3,T=4) grid point
  {
    double const speed = 0.9820851;
    double const density = 16.726301;
    double const temperature = 0.631000;
    // Exponentiate the table value (in log space) at (E=2,rho=3,T=4)
    double const base_val = exp(-5.147755E+01);
    // Calculate the unit conversion (see getEloss function):
    const double multiplier{1000. * density / proton_mass * speed};
    double const reference = multiplier * base_val;
    // KPL: using loose tolerance here so we can compare directly to the data point.
    FAIL_IF_NOT(
        rtt_dsxx::soft_equiv(eloss.getEloss(temperature, density, speed), reference, 1.e-4));
  }

  // Get eloss value for a point between grid points (E=1.5,rho=2.5,T=3.5, i.e. requiring linear
  // interpolation)
  {
    // Halfway between table points 1 and 2:
    double const speed = 0.5 * (0.0982085 + 0.982085);
    // Halfway between table points 2 and 3:
    double const density = 0.5 * (0.167263 + 16.7263);
    // Halfway between table points 3 and 4:
    double const temperature = 0.5 * (0.0158489 + 0.630957);
    // Exponentiate the interpolated table value
    double const base_val = exp(-3.787682963E+01);
    // Calculate the unit conversion (see getEloss function):
    const double multiplier{1000. * density / proton_mass * speed};
    double const reference = multiplier * base_val;
    FAIL_IF_NOT(
        rtt_dsxx::soft_equiv(eloss.getEloss(temperature, density, speed), reference, 1.e-8));
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
