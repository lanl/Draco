//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/test/tstNDI_TNReaction.cc
 * \author Ben R. Ryan
 * \date   2019 Nov 18
 * \brief  NDI_TNReaction test
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi/CDI.hh"
#include "cdi_ndi/NDI_TNReaction.hh"
#include "ds++/Query_Env.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/SystemCall.hh"
#include "ds++/dbc.hh"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using rtt_cdi_ndi::NDI_TNReaction;
using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void gendir_test(rtt_dsxx::UnitTest &ut) {

  // Write a custom gendir file to deal with NDI-required absolute path to data
  std::string gendir_in = "gendir_tmp.all";
  std::string gendir_path = ut.getTestInputPath() + gendir_in;
  std::string data_path = ut.getTestSourcePath() + "ndi_data";
  std::ofstream gendir_tmp_file;
  gendir_tmp_file.open(gendir_path);
  gendir_tmp_file << "tndata\n";
  gendir_tmp_file << "  z=n+be7->p+li7.040ztn  d=12/20/2004  l=lanl04\n";
  gendir_tmp_file << "    f=" << data_path << "\n";
  gendir_tmp_file << "    ft=asc  ln=73  o=3372  end\n";
  gendir_tmp_file.close();

  std::string library_in = "lanl04";
  std::string reaction_in = "n+be7->p+li7";

  std::vector<double> mg_e_bounds = {17.e3, 7.79e3, 2.232e3, 0.184e3, 1.67e-1};
  NDI_TNReaction tn(gendir_path, library_in, reaction_in, mg_e_bounds);

  // Check return values of getters
  FAIL_IF(tn.get_gendir().find(gendir_in) == std::string::npos);
  FAIL_IF_NOT(tn.get_dataset() == "tn");
  FAIL_IF_NOT(tn.get_library() == "lanl04");
  FAIL_IF_NOT(tn.get_reaction() == "n+be7->p+li7");
  FAIL_IF_NOT(tn.get_reaction_name() == "n+be7->p+li7.040ztn");
  FAIL_IF_NOT(tn.get_num_products() == 2);
  auto products = tn.get_products();
  FAIL_IF_NOT(products.size() == 2 && products[0] == 1001 && products[1] == 3007);
  auto multiplicities = tn.get_product_multiplicities();
  FAIL_IF_NOT(multiplicities.size() == 2 && multiplicities[0] == 1 && multiplicities[1] == 1);
  auto reaction_temperature = tn.get_reaction_temperature();
  FAIL_IF_NOT(reaction_temperature.size() == 3);
  FAIL_IF_NOT(soft_equiv(reaction_temperature.front(), 1.000000000000000056e-01, 1.e-8));
  FAIL_IF_NOT(soft_equiv(reaction_temperature.back(), 1.165914400000000045e-01, 1.e-8));
  auto einbar = tn.get_einbar();
  FAIL_IF_NOT(einbar.size() == 3);
  FAIL_IF_NOT(soft_equiv(einbar.front(), 2.968071330000000008e-01, 1.e-8));
  FAIL_IF_NOT(soft_equiv(einbar.back(), 3.458878690000000145e-01, 1.e-8));
  auto sigvbar = tn.get_sigvbar();
  FAIL_IF_NOT(sigvbar.size() == 3);
  FAIL_IF_NOT(soft_equiv(sigvbar.front() / 1.e-23, 7.504850620000000827, 1.e-8));
  FAIL_IF_NOT(soft_equiv(sigvbar.back() / 1.e-23, 7.467488500000000044, 1.e-8));
  FAIL_IF_NOT(soft_equiv(tn.get_reaction_q(), 1.644289999999999964e+03, 1.e-8));
  FAIL_IF_NOT(tn.get_num_groups() == 4);
  auto group_bounds = tn.get_group_bounds();
  FAIL_IF_NOT(group_bounds.size() == 5);
  FAIL_IF_NOT(soft_equiv(group_bounds.front(), 1.700000000000000000e+04, 1.e-8));
  FAIL_IF_NOT(soft_equiv(group_bounds.back(), 1.669999999999999818e-01, 1.e-8));
  auto group_energies = tn.get_group_energies();
  FAIL_IF_NOT(group_energies.size() == 4);
  FAIL_IF_NOT(soft_equiv(group_energies.front(), 1.239500000000000000e+04, 1.e-8));
  FAIL_IF_NOT(soft_equiv(group_energies.back(), 9.208350000000000080e+01, 1.e-8));

  int product_zaid_1 = 1001;            // proton
  int product_zaid_2 = 3007;            // Lithium-7
  double material_temperature = 1.1e-1; // keV

  auto pdf_1 = tn.get_PDF(product_zaid_1, material_temperature);
  for (int n = 0; n < tn.get_num_groups(); n++) {
    if (n == 2) {
      FAIL_IF_NOT(soft_equiv(pdf_1[n], 1., 1.e-8));
    } else {
      FAIL_IF_NOT(soft_equiv(pdf_1[n], 0., 1.e-8));
    }
  }

  auto pdf_2 = tn.get_PDF(product_zaid_2, material_temperature);
  for (int n = 0; n < tn.get_num_groups(); n++) {
    if (n == 3) {
      FAIL_IF_NOT(soft_equiv(pdf_2[n], 1., 1.e-8));
    } else {
      FAIL_IF_NOT(soft_equiv(pdf_2[n], 0., 1.e-8));
    }
  }

  // Check that non-monotonically-decreasing multigroup bounds fail
  try {
    std::vector<double> increasing_mg_e_bounds = {1.e0, 1.e1, 1.e2};
    NDI_TNReaction bad_tn(gendir_path, library_in, reaction_in, increasing_mg_e_bounds);
    FAILMSG("Did not catch expected assertion for non-monotonically-decreasing "
            "multigroup bounds");
  } catch (...) {
    PASSMSG("Expected assertion caught for non-monotonically-decreasing "
            "multigroup bounds");
  }

  if (ut.numFails == 0) {
    PASSMSG("NDI_TNReaction test passes.");
  } else {
    FAILMSG("NDI_TNReaction test fails.");
  }
}

//------------------------------------------------------------------------------------------------//
void gendir_default_test(rtt_dsxx::UnitTest &ut) {
  std::string library_in = "lanl04";
  std::string reaction_in = "d+t->n+a";

  std::vector<double> mg_e_bounds = {17.e3, 7.79e3, 2.232e3, 0.184e3, 1.67e-1};
  NDI_TNReaction tn(library_in, reaction_in, mg_e_bounds);

  if (ut.numFails == 0) {
    PASSMSG("NDI_TNReaction (default gendir path) test passes.");
  } else {
    FAILMSG("NDI_TNReaction (default gendir path) test fails.");
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
      gendir_test(ut);
      gendir_default_test(ut);
    } else {
      PASSMSG("==> ENV{NDI_GENDIR_PATH} not set. Some tests were not run.");
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstNDI_TNReaction.cc
//------------------------------------------------------------------------------------------------//
