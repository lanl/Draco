//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/test/tstNDI_TNReaction.cc
 * \author Ben R. Ryan
 * \brief  NDI_TNReaction test
 * \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Query_Env.hh"
#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "cdi_ndi/NDI_TNReaction.hh"
#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

using rtt_cdi_ndi::NDI_TNReaction;
using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

// cSpell:ignore tndata, einbar

namespace {

void gendir_test(rtt_dsxx::UnitTest &ut) {

  // Write a custom gendir file to deal with NDI-required absolute path to data
  std::string const gendir_in = "gendir_tmp.all";
  std::string const gendir_path = ut.getTestInputPath() + gendir_in;
  std::string const data_path = ut.getTestSourcePath() + "ndi_data";
  std::ofstream gendir_tmp_file;
  gendir_tmp_file.open(gendir_path);
  gendir_tmp_file << "tndata\n";
  gendir_tmp_file << "  z=n+be7->p+li7.040ztn  d=12/20/2004  l=lanl04\n";
  gendir_tmp_file << "    f=" << data_path << "\n";
  gendir_tmp_file << "    ft=asc  ln=73  o=3372  end\n";
  gendir_tmp_file.close();

  std::string const library_in = "lanl04";
  std::string const reaction_in = "n+be7->p+li7";

  std::vector<double> const mg_e_bounds = {17.e3, 7.79e3, 2.232e3, 0.184e3, 1.67e-1};
  NDI_TNReaction const tn(gendir_path, library_in, reaction_in, mg_e_bounds);

  int const product_zaid_1 = 1001; // proton
  int const product_zaid_2 = 3007; // Lithium-7

  // Check return values of getters
  FAIL_IF(tn.get_gendir().find(gendir_in) == std::string::npos);
  FAIL_IF_NOT(tn.get_dataset() == "tn");
  FAIL_IF_NOT(tn.get_library() == "lanl04");
  FAIL_IF_NOT(tn.get_reaction() == "n+be7->p+li7");
  FAIL_IF_NOT(tn.get_reaction_name() == "n+be7->p+li7.040ztn");
  FAIL_IF_NOT(tn.get_num_products() == 2);
  auto products = tn.get_products();
  FAIL_IF_NOT(products.size() == 2 && products[0] == product_zaid_1 &&
              products[1] == product_zaid_2);
  auto multiplicities = tn.get_product_multiplicities();
  FAIL_IF_NOT(multiplicities.size() == 2 && multiplicities[0] == 1 && multiplicities[1] == 1);
  auto reaction_temperature = tn.get_reaction_temperature();
  FAIL_IF_NOT(reaction_temperature.size() == 3);
  FAIL_IF_NOT(soft_equiv(reaction_temperature.front(), 1.000000000000000056e-01, 1.e-8));
  FAIL_IF_NOT(soft_equiv(reaction_temperature.back(), 1.165914400000000045e-01, 1.e-8));
  auto einbar = tn.get_einbar();
  FAIL_IF_NOT(einbar.size() == 3);
  FAIL_IF_NOT(soft_equiv(einbar.front(), 2.968071330000000008e-04, 1.e-8));
  FAIL_IF_NOT(soft_equiv(einbar.back(), 3.458878690000000145e-04, 1.e-8));
  auto sigvbar = tn.get_sigvbar();
  FAIL_IF_NOT(sigvbar.size() == 3);
  FAIL_IF_NOT(soft_equiv(sigvbar.front() / 1.e-23, 7.504850620000000827, 1.e-8));
  FAIL_IF_NOT(soft_equiv(sigvbar.back() / 1.e-23, 7.467488500000000044, 1.e-8));
  FAIL_IF_NOT(soft_equiv(tn.get_reaction_q(), 1.644289999999999964e+03, 1.e-8));
  FAIL_IF_NOT(tn.get_num_groups(product_zaid_1) == 4);
  FAIL_IF_NOT(tn.get_num_groups(product_zaid_2) == 4);
  auto group_bounds = tn.get_mg_e_bounds(product_zaid_1);
  FAIL_IF_NOT(group_bounds.size() == 5);
  FAIL_IF_NOT(soft_equiv(group_bounds.front(), 1.700000000000000000e+01, 1.e-8));
  FAIL_IF_NOT(soft_equiv(group_bounds.back(), 1.669999999999999818e-04, 1.e-8));
  auto group_energies = tn.get_group_energies(product_zaid_2);
  FAIL_IF_NOT(group_energies.size() == 4);
  FAIL_IF_NOT(soft_equiv(group_energies.front(), 1.239500000000000000e+01, 1.e-8));
  FAIL_IF_NOT(soft_equiv(group_energies.back(), 9.208350000000000080e-02, 1.e-8));

  double const material_temperature = 1.1e-1; // keV

  const auto [ebins_1, pdf_1] = tn.get_ebins_and_PDF(product_zaid_1, material_temperature);
  for (int n = 0; n < tn.get_num_groups(product_zaid_1); n++) {
    if (n == 2) {
      FAIL_IF_NOT(soft_equiv(pdf_1[n], 1., 1.e-8));
    } else {
      FAIL_IF_NOT(soft_equiv(pdf_1[n], 0., 1.e-8));
    }
  }

  const auto [ebins_2, pdf_2] = tn.get_ebins_and_PDF(product_zaid_2, material_temperature);
  for (int n = 0; n < tn.get_num_groups(product_zaid_2); n++) {
    if (n == 3) {
      FAIL_IF_NOT(soft_equiv(pdf_2[n], 1., 1.e-8));
    } else {
      FAIL_IF_NOT(soft_equiv(pdf_2[n], 0., 1.e-8));
    }
  }

  // Check that non-monotonically-decreasing multigroup bounds fail
  try {
    std::vector<double> const increasing_mg_e_bounds = {1.e0, 1.e1, 1.e2};
    NDI_TNReaction const bad_tn(gendir_path, library_in, reaction_in, increasing_mg_e_bounds);
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
  std::string const library_in = "lanl04";
  std::string reaction_in = "d+t->n+a";

  std::vector<double> const mg_e_bounds = {17.e3, 7.79e3, 2.232e3, 0.184e3, 1.67e-1};
  NDI_TNReaction const tn(library_in, reaction_in, mg_e_bounds);

  if (ut.numFails == 0) {
    PASSMSG("NDI_TNReaction (default gendir path) test passes.");
  } else {
    FAILMSG("NDI_TNReaction (default gendir path) test fails.");
  }
}

//------------------------------------------------------------------------------------------------//
void gendir_default_bounds_test(rtt_dsxx::UnitTest &ut) {
  std::string const library_in = "lanl04";
  std::string reaction_in = "d+d->p+t";

  // Don't set bounds - let the class set them
  NDI_TNReaction const tn(library_in, reaction_in);

  // Retrieve the PDF for alphas at 10 keV and compute an approximate mean:
  constexpr int prot{1001};
  constexpr int trit{1003};
  const auto p_group_bounds = tn.get_mg_e_bounds(prot);
  const auto t_group_bounds = tn.get_mg_e_bounds(trit);

  constexpr std::array<double, 6> temperatures{1.0, 5.0, 10.0, 50.0, 100.0, 500.0};
  double last_p_avg_e{0.0};
  double last_t_avg_e{0.0};
  for (auto &t : temperatures) {
    const auto [p_e_bins, p_pdf] = tn.get_ebins_and_PDF(prot, t);
    const auto [t_e_bins, t_pdf] = tn.get_ebins_and_PDF(trit, t);
    // Check bins vs. PDF size:
    FAIL_IF_NOT(p_e_bins.size() == p_pdf.size() + 1);
    FAIL_IF_NOT(t_e_bins.size() == t_pdf.size() + 1);
    // Check PDF normalization:
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(std::accumulate(p_pdf.begin(), p_pdf.end(), 0.0), 1.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(std::accumulate(t_pdf.begin(), t_pdf.end(), 0.0), 1.0));
    // Soft check to be sure average product E increases with higher plasma T:
    double p_avg_e{0.0};
    double t_avg_e{0.0};
    for (size_t bin = 0; bin < p_pdf.size(); bin++) {
      p_avg_e += p_pdf[bin] * 0.5 * (p_e_bins[bin] + p_e_bins[bin + 1]);
    }
    for (size_t bin = 0; bin < t_pdf.size(); bin++) {
      t_avg_e += t_pdf[bin] * 0.5 * (t_e_bins[bin] + t_e_bins[bin + 1]);
    }
    // Triton should be less energetic than proton:
    FAIL_IF(t_avg_e > p_avg_e);
    // Product avg energies should increase:
    FAIL_IF_NOT(p_avg_e > last_p_avg_e);
    FAIL_IF_NOT(t_avg_e > last_t_avg_e);

    last_p_avg_e = p_avg_e;
    last_t_avg_e = t_avg_e;
  }

  if (ut.numFails == 0) {
    PASSMSG("NDI_TNReaction (default gendir path and MG bounds) test passes.");
  } else {
    FAILMSG("NDI_TNReaction (default gendir path and MG bounds) test fails.");
  }
}
} // unnamed namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    std::string gendir_default;
    bool def_gendir{false};
    std::tie(def_gendir, gendir_default) = rtt_dsxx::get_env_val<std::string>("NDI_GENDIR_PATH");

    if (def_gendir && std::filesystem::exists(gendir_default)) {
      gendir_test(ut);
      gendir_default_test(ut);
      gendir_default_bounds_test(ut);
    } else {
      PASSMSG("==> ENV{NDI_GENDIR_PATH} not set. Some tests were not run.");
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstNDI_TNReaction.cc
//------------------------------------------------------------------------------------------------//
