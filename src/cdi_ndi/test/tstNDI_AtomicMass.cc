//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/test/tstNDI_AtomicMass.cc
 * \author Ben R. Ryan
 * \date   2020 Mar 6
 * \brief  NDI_AtomicMass test
 * \note   Copyright (C) 2020-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Query_Env.hh"
#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "dsxx/path.hh"
#include "cdi_ndi/NDI_AtomicMass.hh"
#include <fstream>
#include <string>
#include <tuple>

// cSpell:ignore mendf

using rtt_cdi_ndi::NDI_AtomicMass;
using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void amw_test(rtt_dsxx::UnitTest &ut) {

  // Write a custom gendir file to deal with NDI-required absolute path to data
  std::string const gendir_in = "gendir_tmp.all";
  std::string const gendir_tmp_path = ut.getTestInputPath() + gendir_in;
  std::string const data_path = ut.getTestSourcePath() + "ndi_data";
  std::ofstream gendir_tmp_file;
  gendir_tmp_file.open(gendir_tmp_path);
  gendir_tmp_file << "multigroup_neutron\n";
  gendir_tmp_file << "  z=1001.700nm  d=11/19/2010  l=mendf71x\n";
  gendir_tmp_file << "    f=fake/data/path\n";
  gendir_tmp_file << "    ft=asc  ln=2  o=28\n";
  gendir_tmp_file << "    ng=618  t=2.5300642359999999e-08  s0=10000000000\n";
  gendir_tmp_file << "    aw=1.0078249887344399  awr=0.99916729999999998  end\n";
  gendir_tmp_file.close();

  NDI_AtomicMass const ndi_amw(gendir_tmp_path);

  int const proton_zaid = 1001;
  double const proton_amw = ndi_amw.get_amw(proton_zaid);

  FAIL_IF_NOT(soft_equiv(proton_amw, 1.673532619349924133e-24, 1.e-8));

  int const electron_zaid = -1;
  double const electron_amw = ndi_amw.get_amw(electron_zaid);

  FAIL_IF_NOT(soft_equiv(electron_amw, 9.109382909999999302e-28, 1.e-8));

  if (ut.numFails == 0) {
    PASSMSG("NDI_AtomicMass test passes.");
  } else {
    FAILMSG("NDI_AtomicMass test fails.");
  }
}

//------------------------------------------------------------------------------------------------//
void amw_default_test(rtt_dsxx::UnitTest &ut) {

  NDI_AtomicMass const ndi_amw;

  int const proton_zaid = 1001;
  double const proton_amw = ndi_amw.get_amw(proton_zaid);

  // This might change between NDI releases hence the low tolerance
  FAIL_IF_NOT(soft_equiv(proton_amw, 1.673532619349924133e-24, 1.e-4));

  int const electron_zaid = -1;
  double const electron_amw = ndi_amw.get_amw(electron_zaid);

  FAIL_IF_NOT(soft_equiv(electron_amw, 9.109382909999999302e-28, 1.e-8));

  if (ut.numFails == 0) {
    PASSMSG("NDI_AtomicMass (default gendir path) test passes.");
  } else {
    FAILMSG("NDI_AtomicMass (default_gendir path) test fails.");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    amw_test(ut);
    std::string gendir_default;
    bool def_gendir{false};
    std::tie(def_gendir, gendir_default) = rtt_dsxx::get_env_val<std::string>("NDI_GENDIR_PATH");

    if (def_gendir && rtt_dsxx::fileExists(gendir_default)) {
      amw_default_test(ut);
    } else {
      PASSMSG("==> ENV{NDI_GENDIR_PATH} not set. Some tests were not run.");
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstNDI_AtomicMass.cc
//------------------------------------------------------------------------------------------------//
