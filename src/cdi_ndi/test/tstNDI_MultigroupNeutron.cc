//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/test/tstNDI_MultigroupNeutron.cc
 * \author Kyle S. Beling
 * \date   2025 March 18
 * \brief  NDI_MultigroupNeutron test
 * \note   Copyright (C) 2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Query_Env.hh"
#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/UnitTest.hh"
#include "dsxx/path.hh"
#include "cdi_ndi/NDI_MultigroupNeutron.hh"
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

using rtt_cdi_ndi::NDI_MultigroupNeutron;
using rtt_dsxx::soft_equiv;

void gendir_test(rtt_dsxx::UnitTest &ut) {

  // Write a custom gendir file to deal with NDI-required absolute path to data
  std::string const gendir_in = "gendir_tmp.all";
  std::string const gendir_path = ut.getTestInputPath() + gendir_in;
  std::string const data_path = ut.getTestSourcePath() + "001-H-002_example";
  std::ofstream gendir_tmp_file;
  gendir_tmp_file.open(gendir_path);
  gendir_tmp_file << "multigroup_neutron\n";
  gendir_tmp_file << "  z=1002.011nm d=10/20/2011 l=cp2011\n";
  gendir_tmp_file << "  f=" << data_path << "\n";
  gendir_tmp_file << "  ft=asc ln=2 o=28 ng=618 aw=2.01410208030720E+00\n";
  gendir_tmp_file << "  t=2.530064236e-08 s0=10000000000.0 awr=1.9968 end\n";
  gendir_tmp_file.close();

  std::string const library_in = "cp2011";
  const int target_in = 1002;

  const std::vector<double> energy_grid{1.0e+1, 1.0e+0, 1.0e-1};
  NDI_MultigroupNeutron const neutron(gendir_path, library_in, target_in, energy_grid);

  // GENERAL DATA INFO
  FAIL_IF(neutron.get_gendir().find(gendir_in) == std::string::npos);
  FAIL_IF_NOT(neutron.get_dataset() == "multigroup_neutron");
  FAIL_IF_NOT(neutron.get_library() == "cp2011");
  FAIL_IF_NOT(neutron.targetZAID() == 1002);

  // GROUP STRUCTURE
  auto group_bounds = neutron.groupBounds();
  FAIL_IF_NOT(group_bounds.size() == 3);
  FAIL_IF_NOT(soft_equiv(group_bounds[0], 1.0000e+01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(group_bounds[1], 9.9013e-01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(group_bounds[2], 1.0115e-01, 1.0e-4));

  // GROUP VELOCITIES
  auto group_velocities = neutron.groupVelocities();
  FAIL_IF_NOT(group_velocities.size() == 2);
  FAIL_IF_NOT(soft_equiv(group_velocities[0], 1.8699e+01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(group_velocities[1], 8.0359e+00, 1.0e-4));

  //--------------------------------------------------
  // TOTAL CROSS SECTION DATA
  //--------------------------------------------------
  // MULTIGROUP TOTAL CROSS SECTION
  auto sigma_total = neutron.totalCrossSection();
  FAIL_IF_NOT(sigma_total.size() == 2);
  FAIL_IF_NOT(soft_equiv(sigma_total[0], 2.4898e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(sigma_total[1], 3.0278e+00, 1.0e-4));

  // ENERGY DEPENDENT TOTAL CROSS SECTION
  FAIL_IF_NOT(soft_equiv(neutron.totalCrossSection(5.0e+0), 2.4898e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(neutron.totalCrossSection(5.0e-1), 3.0278e+00, 1.0e-4));

  // TOTAL GROUP KERMAS
  auto group_kermas = neutron.kermas();
  FAIL_IF_NOT(group_kermas.size() == 2);
  FAIL_IF_NOT(soft_equiv(group_kermas[0], 2.3030e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(group_kermas[1], 7.2514e-01, 1.0e-4));

  // TOTAL SCATTERING MATRIX
  auto scattering_matrix = neutron.scatteringMatrix();
  FAIL_IF_NOT(scattering_matrix.size() == 5);
  FAIL_IF_NOT(scattering_matrix[0].size() == 2);
  FAIL_IF_NOT(scattering_matrix[0][0].size() == 2);

  // l = 0
  FAIL_IF_NOT(soft_equiv(scattering_matrix[0][0][0], 9.8462e-01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix[0][0][1], 1.5099e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix[0][1][0], 0.0000e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix[0][1][1], 3.0278e+00, 1.0e-4));

  // l = 4
  FAIL_IF_NOT(soft_equiv(scattering_matrix[4][0][0], 6.3679e-02, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix[4][0][1], 1.3396e-03, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix[4][1][0], 0.0000e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix[4][1][1], 1.1539e-02, 1.0e-4));

  // TOTAL PDF
  auto pdf = neutron.pdf(5.0e+0);
  FAIL_IF_NOT(pdf.size() == 2);
  FAIL_IF_NOT(soft_equiv(pdf[0], 3.9471e-01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(pdf[1], 6.0529e-01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(pdf[0] + pdf[1], 1.0000e+00, 1.0e-4));

  //--------------------------------------------------
  // REACTION SPECIFIC DATA
  //--------------------------------------------------
  auto reactions = neutron.listOfReactions();
  FAIL_IF_NOT(reactions.size() == 3);
  FAIL_IF_NOT(reactions[0] == 2 && reactions[1] == 16 && reactions[2] == 102);

  // MT = 2
  auto products_mt2 = neutron.listOfProducts(2);
  FAIL_IF_NOT(products_mt2.size() == 2);
  FAIL_IF_NOT(products_mt2[0] == 1 && products_mt2[1] == 1002);
  auto multiplicites_mt2 = neutron.productMultiplicities(2);
  FAIL_IF_NOT(multiplicites_mt2.size() == 2);
  FAIL_IF_NOT(multiplicites_mt2[0] == 1 && multiplicites_mt2[1] == 1);
  FAIL_IF_NOT(soft_equiv(neutron.qValue(2), 0.0000e+00, 1.0e-4));
  auto sigma_total_mt2 = neutron.totalCrossSection(2);
  FAIL_IF_NOT(sigma_total_mt2.size() == 2);
  FAIL_IF_NOT(soft_equiv(sigma_total_mt2[0], 2.4851e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(sigma_total_mt2[1], 3.0278e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(neutron.totalCrossSection(2, 5.0e+0), 2.4851e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(neutron.totalCrossSection(2, 5.0e-1), 3.0278e+00, 1.0e-4));

  // MT = 16
  auto products_mt16 = neutron.listOfProducts(16);
  FAIL_IF_NOT(products_mt16.size() == 2);
  FAIL_IF_NOT(products_mt16[0] == 1 && products_mt16[1] == 1001);
  auto multiplicites_mt16 = neutron.productMultiplicities(16);
  FAIL_IF_NOT(multiplicites_mt16.size() == 2);
  FAIL_IF_NOT(multiplicites_mt16[0] == 2 && multiplicites_mt16[1] == 1);
  FAIL_IF_NOT(soft_equiv(neutron.qValue(16), -2.2250e+00, 1.0e-4));
  auto sigma_total_mt16 = neutron.totalCrossSection(16);
  FAIL_IF_NOT(sigma_total_mt16.size() == 2);
  FAIL_IF_NOT(soft_equiv(sigma_total_mt16[0], 4.6933e-03, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(sigma_total_mt16[1], 0.0000e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(neutron.totalCrossSection(16, 5.0e+0), 4.6933e-03, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(neutron.totalCrossSection(16, 5.0e-1), 0.0000e+00, 1.0e-4));

  // MT = 102
  auto products_mt102 = neutron.listOfProducts(102);
  FAIL_IF_NOT(products_mt102.size() == 1 && products_mt102[0] == 1003);
  auto multiplicites_mt102 = neutron.productMultiplicities(102);
  FAIL_IF_NOT(multiplicites_mt16.size() == 2 && multiplicites_mt102[0] == 1);
  FAIL_IF_NOT(soft_equiv(neutron.qValue(102), 6.2574e+00, 1.0e-4));
  auto sigma_total_mt102 = neutron.totalCrossSection(102);
  FAIL_IF_NOT(sigma_total_mt102.size() == 2);
  FAIL_IF_NOT(soft_equiv(sigma_total_mt102[0], 7.7065e-06, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(sigma_total_mt102[1], 3.8919e-06, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(neutron.totalCrossSection(102, 5.0e+0), 7.7065e-06, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(neutron.totalCrossSection(102, 5.0e-1), 3.8919e-06, 1.0e-4));

  //--------------------------------------------------
  // PRODUCT SPECIFIC DATA
  //--------------------------------------------------
  auto products = neutron.listOfAllProducts();
  FAIL_IF_NOT(products.size() == 4);
  FAIL_IF_NOT(products[0] == 0 && products[1] == 1001 && products[2] == 1002 &&
              products[3] == 1003);

  // ZAID = 0
  auto kermas_0 = neutron.kermas(0);
  FAIL_IF_NOT(kermas_0.size() == 2);
  FAIL_IF_NOT(soft_equiv(kermas_0[0], 6.0704e-05, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(kermas_0[1], 2.5736e-05, 1.0e-4));

  auto scattering_matrix_0 = neutron.scatteringMatrix(0);
  FAIL_IF_NOT(scattering_matrix_0.size() == 1);
  FAIL_IF_NOT(scattering_matrix_0[0].size() == 2);
  FAIL_IF_NOT(scattering_matrix_0[0][0].size() == 2);

  // l = 0
  FAIL_IF_NOT(soft_equiv(scattering_matrix_0[0][0][0], 7.7065e-06, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_0[0][0][1], 0.0000e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_0[0][1][0], 3.8919e-06, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_0[0][1][1], 0.0000e+00, 1.0e-4));

  // pdf
  auto pdf_0 = neutron.pdf(0, 5.0e+0);
  FAIL_IF_NOT(pdf_0.size() == 2);
  FAIL_IF_NOT(soft_equiv(pdf_0[0], 1.0000e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(pdf_0[1], 0.0000e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(pdf_0[0] + pdf_0[1], 1.0000e+00, 1.0e-4));

  // ZAID = 1001
  auto kermas_1001 = neutron.kermas(1001);
  FAIL_IF_NOT(kermas_1001.size() == 2);
  FAIL_IF_NOT(soft_equiv(kermas_1001[0], 6.3591e-03, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(kermas_1001[1], 0.0000e+00, 1.0e-4));

  auto scattering_matrix_1001 = neutron.scatteringMatrix(1001);
  FAIL_IF_NOT(scattering_matrix_1001.size() == 5);
  FAIL_IF_NOT(scattering_matrix_1001[0].size() == 2);
  FAIL_IF_NOT(scattering_matrix_1001[0][0].size() == 2);

  // l = 0
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1001[0][0][0], 2.4919e-03, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1001[0][0][1], 2.2014e-03, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1001[0][1][0], 0.0000e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1001[0][1][1], 0.0000e+00, 1.0e-4));

  // l = 4
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1001[4][0][0], 2.7970e-04, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1001[4][0][1], 8.7211e-05, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1001[4][1][0], 0.0000e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1001[4][1][1], 0.0000e+00, 1.0e-4));

  // pdf
  auto pdf_1001 = neutron.pdf(1001, 5.0e+0);
  FAIL_IF_NOT(pdf_1001.size() == 2);
  FAIL_IF_NOT(soft_equiv(pdf_1001[0], 5.3096e-01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(pdf_1001[1], 4.6904e-01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(pdf_1001[0] + pdf_1001[1], 1.0000e+00, 1.0e-4));

  // ZAID = 1002
  auto kermas_1002 = neutron.kermas(1002);
  FAIL_IF_NOT(kermas_1002.size() == 2);
  FAIL_IF_NOT(soft_equiv(kermas_1002[0], 2.2966e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(kermas_1002[1], 7.2513e-01, 1.0e-4));

  auto scattering_matrix_1002 = neutron.scatteringMatrix(1002);
  FAIL_IF_NOT(scattering_matrix_1002.size() == 5);
  FAIL_IF_NOT(scattering_matrix_1002[0].size() == 2);
  FAIL_IF_NOT(scattering_matrix_1002[0][0].size() == 2);

  // l = 0
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1002[0][0][0], 9.1313e-01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1002[0][0][1], 1.5720e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1002[0][1][0], 0.0000e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1002[0][1][1], 3.0278e+00, 1.0e-4));

  // l = 4
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1002[4][0][0], -3.2589e-02, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1002[4][0][1], 1.1382e-01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1002[4][1][0], 0.0000e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1002[4][1][1], -2.0629e-01, 1.0e-4));

  // pdf
  auto pdf_1002 = neutron.pdf(1002, 5.0e+0);
  FAIL_IF_NOT(pdf_1002.size() == 2);
  FAIL_IF_NOT(soft_equiv(pdf_1002[0], 3.6744e-01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(pdf_1002[1], 6.3256e-01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(pdf_1002[0] + pdf_1002[1], 1.0000e+00, 1.0e-4));

  // ZAID = 1003
  auto kermas_1003 = neutron.kermas(1003);
  FAIL_IF_NOT(kermas_1003.size() == 2);
  FAIL_IF_NOT(soft_equiv(kermas_1003[0], 6.3651e-06, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(kermas_1003[1], 7.5438e-07, 1.0e-4));

  auto scattering_matrix_1003 = neutron.scatteringMatrix(1003);
  FAIL_IF_NOT(scattering_matrix_1003.size() == 5);
  FAIL_IF_NOT(scattering_matrix_1003[0].size() == 2);
  FAIL_IF_NOT(scattering_matrix_1003[0][0].size() == 2);

  // l = 0
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1003[0][0][0], 1.9821e-06, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1003[0][0][1], 5.7243e-06, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1003[0][1][0], 0.0000e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1003[0][1][1], 3.8919e-06, 1.0e-4));

  // l = 4
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1003[4][0][0], 1.8826e-06, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1003[4][0][1], 5.1549e-06, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1003[4][1][0], 0.0000e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(scattering_matrix_1003[4][1][1], 2.6085e-06, 1.0e-4));

  // pdf
  auto pdf_1003 = neutron.pdf(1003, 5.0e+0);
  FAIL_IF_NOT(pdf_1003.size() == 2);
  FAIL_IF_NOT(soft_equiv(pdf_1003[0], 2.5720e-01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(pdf_1003[1], 7.4280e-01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(pdf_1003[0] + pdf_1003[1], 1.0000e+00, 1.0e-4));

  //--------------------------------------------------
  // CHECK IF TEST PASSED
  //--------------------------------------------------
  if (ut.numFails == 0) {
    PASSMSG("NDI_MultigroupNeutron test passes.");
  } else {
    FAILMSG("NDI_MultigroupNeutron test fails.");
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
    } else {
      PASSMSG("==> ENV{NDI_GENDIR_PATH} not set. Some tests were not run.");
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstNDI_MultigroupNeutron.cc
//------------------------------------------------------------------------------------------------//
