//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   compton_interface/test/tCSK_Interface_dummy_data.cc
 * \author Andrew Till
 * \date   2020 Oct 14
 * \brief  Implementation file for tCSK_Interface_dummy_data
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "compton_interface/CSK_Interface.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace rtt_compton_interface_dummy_data_test {

using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

//!  Tests the CSK_Interface constructor and a couple of access routines.
void compton_file_test(rtt_dsxx::UnitTest &ut) {

  // Make true if golds need updating
  const bool do_print = false;
  auto print_lambda = [](const std::vector<double> &variable, const std::string &variable_name) {
    auto print = [](double a) { std::cout << a << ", "; };
    std::cout << variable_name << ":\n";
    std::cout << std::setprecision(14);
    std::cout << '\n';
    std::for_each(variable.begin(), variable.end(), print);
    std::cout << std::endl;
  };

  // Tolerance used for checks
  const double tol = 1e-11;

  // Start the test.

  std::cout << "\n---------------------------------------------------------\n"
            << "   Test Draco code calling CSK_generator routines\n"
            << "---------------------------------------------------------\n";

  // open a small mg opacity file:
  const std::string filename = ut.getTestSourcePath() + "../../compton_tools/test/dummy_data";
  std::cout << "Attempting to construct a CSK_Interface object...\n" << std::endl;
  std::unique_ptr<rtt_compton_interface::CSK_Interface> compton_test;

  try {
    compton_test = std::make_unique<rtt_compton_interface::CSK_Interface>(filename);
  } catch (int /*asrt*/) {
    FAILMSG("Failed to construct a CSK_Interface object!");
    // if construction fails, there is no reason to continue testing...
    return;
  }
  std::cout << "\n(...Success!)" << std::endl;

  // Check some of the data in the CSK_generator-opened file:
  const std::vector<double> grp_bds = compton_test->get_group_bounds();
  const std::vector<double> T_evals = compton_test->get_etemp_pts();

  // Unitless (divided by mec2)
  // NB: These values can be read directly from the 3rd line of the ASCII csk data files
  std::vector<double> grp_bds_gold = {1.57311251e-06, 3.14622503e-04, 7.86556258e-04,
                                      1.57311251e-03, 3.14622503e-02};
  // NB: These values can be read directly from the ASCII csk data files
  std::vector<double> T_evals_gold = {1.57311251e-05, 1.57311251e-04, 3.30353629e-04,
                                      6.60707256e-04};
  // NB: First and last temperature from line 2 of the ASCII csk data files
  std::vector<double> line2_Ts_gold = {1.41580126e-05, 7.26777982e-04};

  // Sizes
  const size_t num_groups_gold = grp_bds_gold.size() - 1U;
  const size_t num_T_evals_gold = T_evals_gold.size();
  const size_t num_evals_gold = 4; // out_lin, in_lin, out_nonlin, in_nonlin
  const size_t num_leg_moments_gold = 2;

  // Multiply by electron rest-mass energy (keV; using CSK value)
  const double mec2 = 510.998;
  for (double &grp_bd : grp_bds_gold) {
    grp_bd *= mec2;
  }
  // Interface does not scale temperatures, so no multiplication needed:
  /*
  for (size_t i = 0; i < T_evals_gold.size(); ++i)
  {
      T_evals_gold[i] *= mec2;
  }
  line2_Ts_gold[0] *= mec2;
  line2_Ts_gold[1] *= mec2;
  */

  FAIL_IF(grp_bds.size() != (num_groups_gold + 1U));
  ut.check(soft_equiv(grp_bds, grp_bds_gold, tol), "checked group boundaries");

  FAIL_IF(T_evals.size() != num_T_evals_gold);
  ut.check(soft_equiv(T_evals, T_evals_gold, tol), "checked temperature grid");

  FAIL_IF_NOT(soft_equiv(compton_test->get_min_etemp(), line2_Ts_gold[0]));
  FAIL_IF_NOT(soft_equiv(compton_test->get_max_etemp(), line2_Ts_gold[1]));

  if (ut.numFails == 0) {
    std::cout << "\nCorrectly read group bounds and electron temps!" << std::endl;
  }

  // Test data retrieval: interpolate to a grid point in temperature
  {
    const double interp_T = T_evals[num_T_evals_gold - 1U];
    std::cout << "Testing interpolation at T = " << (interp_T * mec2) << " keV\n";

    std::vector<std::vector<std::vector<std::vector<double>>>> interp_data =
        compton_test->interpolate_csk(interp_T);

    const size_t sz0 = interp_data.size();
    const size_t sz1 = interp_data[0].size();
    const size_t sz2 = interp_data[0][0].size();
    const size_t sz3 = interp_data[0][0][0].size();
    const size_t tot_size = sz0 * sz1 * sz2 * sz3;
    std::vector<double> flat_interp_data(tot_size, 0.0);

    for (size_t i0 = 0; i0 < sz0; ++i0) {
      for (size_t i1 = 0; i1 < sz1; ++i1) {
        for (size_t i2 = 0; i2 < sz2; ++i2) {
          for (size_t i3 = 0; i3 < sz3; ++i3) {
            const size_t i = i3 + sz3 * (i2 + sz2 * (i1 + sz1 * (i0)));
            flat_interp_data[i] = interp_data[i0][i1][i2][i3];
          }
        }
      }
    }

    // Test sizes
    ut.check(sz0 == num_evals_gold, "tested evals size");
    ut.check((sz1 == num_groups_gold) && (sz2 == num_groups_gold), "tested groups size");
    ut.check(sz3 == num_leg_moments_gold, "tested Legendre moments size");

    // Print result (useful if golds need updating)
    if (do_print) {
      print_lambda(flat_interp_data, "flat_interp_data");
    }

    // NB: These values come directly from the ASCII data files
    std::vector<double> flat_interp_gold = {
        // out_lin
        1.023686968316, 0.003573655955675, 0.1408569058113, -0.003385368827333, 0, 0, 0, 0,
        0.01368733137026, -0.000337289996441, 0.8580132416506, 0.00346237220459, 0.1150544176069,
        -0.002695776219916, 0, 0, 0, 0, 0.03016927525143, -0.000745487878607, 0.7395544180756,
        0.003537535339196, 0.08790804121747, -0.001974567477495, 0, 0, 0, 0, 0.02760683687763,
        -0.0006747774819197, 0.7154380590368, 0.002500085526454,
        // in_lin
        1.029549163582, 0.003360819879047, 0.1682617986774, -0.004498320348364, 0, 0, 0, 0,
        0.01155720790042, -0.000254680100335, 0.8615393494377, 0.003337305544177, 0.1370933267956,
        -0.003577889661548, 0, 0, 0, 0, 0.02546236651279, -0.0005628307049739, 0.74122185377,
        0.003480918190072, 0.104392775669, -0.002618892286625, 0, 0, 0, 0, 0.02330836203736,
        -0.0005092109709329, 0.7158405069922, 0.002485525111488,
        // out_nonlin
        2.417350052736e-20, 5.333857003773e-23, 1.569215319107e-22, -3.474039416897e-24, 0, 0, 0, 0,
        2.025140617993e-22, -5.27716088378e-24, 4.991731616824e-22, 1.576887628939e-24,
        7.937214293507e-24, -1.661082680445e-25, 0, 0, 0, 0, 1.15611911544e-23, -3.08379217909e-25,
        2.906510852958e-23, 1.161410570085e-25, 3.652307822369e-25, -6.9213492314e-27, 0, 0, 0, 0,
        6.471151078044e-25, -1.788886184269e-26, 9.067645193159e-25, 3.288019868804e-27,
        // in_nonlin
        2.412582001649e-20, 5.506360234252e-23, 1.858984933376e-22, -4.598103855448e-24, 0, 0, 0, 0,
        1.70028840385e-22, -3.978450052949e-24, 4.964897275093e-22, 1.671777660171e-24,
        9.356135881693e-24, -2.194753402066e-25, 0, 0, 0, 0, 9.680369923644e-24, -2.32267253788e-25,
        2.875946793361e-23, 1.265594484598e-25, 4.274264612141e-25, -9.142968009692e-27, 0, 0, 0, 0,
        5.392751048597e-25, -1.344066881751e-26, 8.87128629256e-25, 3.986007245261e-27};

    ut.check(soft_equiv(flat_interp_data, flat_interp_gold, tol), "checked data retrieval");
  }

  // Test interpolation
  {
    const double alpha = 0.4;
    const double interp_T = alpha * T_evals[0] + (1.0 - alpha) * T_evals[1];
    std::cout << "Testing interpolation at T = " << (interp_T * mec2) << " keV\n";

    std::vector<std::vector<std::vector<std::vector<double>>>> interp_data =
        compton_test->interpolate_csk(interp_T);

    const size_t sz0 = interp_data.size();
    const size_t sz1 = interp_data[0].size();
    const size_t sz2 = interp_data[0][0].size();
    const size_t sz3 = interp_data[0][0][0].size();
    const size_t tot_size = sz0 * sz1 * sz2 * sz3;
    std::vector<double> flat_interp_data(tot_size, 0.0);

    for (size_t i0 = 0; i0 < sz0; ++i0) {
      for (size_t i1 = 0; i1 < sz1; ++i1) {
        for (size_t i2 = 0; i2 < sz2; ++i2) {
          for (size_t i3 = 0; i3 < sz3; ++i3) {
            const size_t i = i3 + sz3 * (i2 + sz2 * (i1 + sz1 * (i0)));
            flat_interp_data[i] = interp_data[i0][i1][i2][i3];
          }
        }
      }
    }

    // Test sizes
    ut.check(sz0 == num_evals_gold, "tested evals size");
    ut.check((sz1 == num_groups_gold) && (sz2 == num_groups_gold), "tested groups size");
    ut.check(sz3 == num_leg_moments_gold, "tested Legendre moments size");

    // Print result (useful if golds need updating)
    if (do_print) {
      print_lambda(flat_interp_data, "flat_interp_data");
    }

    std::vector<double> flat_interp_gold = {
        //out_lin
        1.4895048336024, 0.0013157633414251, 0.04145112467826, -0.00099234375153514, 0, 0, 0, 0,
        0.048935718979052, -0.001201759773222, 1.1427518955345, 0.0020642528309125,
        0.015831053793403, -0.00036191022639507, 0, 0, 0, 0, 0.060076056533251, -0.0014978191061645,
        0.95368017957811, 0.0024195173860433, 0.0042167573598475, -3.7346606985974e-05, 0, 0, 0, 0,
        0.12606997604209, -0.0029756533421552, 0.70308466958259, 0.0057984236887092,
        // in_lin
        1.4902075219537, 0.0012890343406445, 0.044522269416302, -0.0011149333854142, 0, 0, 0, 0,
        0.046790157320168, -0.0011165953164094, 1.1422078133786, 0.002071967455139,
        0.016959526258933, -0.00040597846461434, 0, 0, 0, 0, 0.055452234941061, -0.0013152449665463,
        0.95025760553927, 0.002535497503532, 0.0044495495432538, -4.5183679054704e-05, 0, 0, 0, 0,
        0.11605523138559, -0.0025883973572953, 0.67868811344021, 0.0063593874505928,
        // out_nonlin
        4.6630056595372e-18, 8.5759248352589e-23, 1.1673356563462e-22, -2.6167370751723e-24, 0, 0,
        0, 0, 1.5017646436551e-22, -3.9459638253177e-24, 9.1490085641636e-21, 5.3409650219456e-23,
        5.484626607511e-24, -1.0596457235234e-25, 0, 0, 0, 0, 1.0314177969686e-23,
        -2.9818990024258e-25, 4.2638762788237e-22, 8.5042714161135e-25, 1.8073497542315e-25,
        -2.4297847629037e-27, 0, 0, 0, 0, 7.6278018408129e-25, -2.5980363364632e-26,
        1.5794222622867e-23, 1.5887774510949e-25,
        // in_nonlin
        4.6615665467427e-18, 1.7873552446889e-22, 1.2663879157434e-22, -2.9982049815629e-24, 0, 0,
        0, 0, 1.3770981199909e-22, -3.4369299797817e-24, 9.1220263142097e-21, 5.4325270040856e-23,
        5.890752806507e-24, -1.2076068345897e-25, 0, 0, 0, 0, 9.5432783203125e-24,
        -2.6517064022462e-25, 4.2104110096026e-22, 1.0546962404339e-24, 1.9308590430364e-25,
        -2.8485633209286e-27, 0, 0, 0, 0, 7.0910968932275e-25, -2.3423156364316e-26,
        1.4873126610905e-23, 1.8374675563093e-25};

    ut.check(soft_equiv(flat_interp_data, flat_interp_gold, tol), "checked data interpolation");
  }

  if (ut.numFails == 0)
    std::cout << "\nCorrectly read multigroup data points!" << std::endl;

  if (ut.numFails == 0) {
    PASSMSG("Successfully linked Draco against CSK_generator.");
  } else {
    FAILMSG("Did not successfully link Draco against CSK_generator.");
  }
}

//------------------------------------------------------------------------------------------------//
//!  Tests CSK_Interface's error-handling on a non-existent file.
void compton_fail_test(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n---------------------------------------------------------\n"
            << "    Test CSK_Interface bad file handling    \n"
            << "---------------------------------------------------------\n";
  // open a small mg opacity file:
  std::string filename = ut.getTestSourcePath() + "non_existent.compton";
  std::cout << "Testing with a non-existent file...\n" << std::endl;
  std::unique_ptr<rtt_compton_interface::CSK_Interface> compton_test;

  bool caught = false;
  try {
    compton_test = std::make_unique<rtt_compton_interface::CSK_Interface>(filename);
  } catch (rtt_dsxx::assertion &asrt) {
    std::cout << "Draco exception thrown: " << asrt.what() << std::endl;
    // We successfully caught the bad file!
    caught = true;
  } catch (const int & /*asrt*/) {
    std::cout << "CSK exception thrown. " << std::endl;
    // We successfully caught the bad file!
    caught = true;
  }

  if (!caught)
    ITFAILS;

  if (ut.numFails == 0) {
    PASSMSG("Successfully caught a CSK_generator exception.");
  } else {
    FAILMSG("Did not successfully catch a CSK_generator exception.");
  }
}

} // namespace rtt_compton_interface_dummy_data_test

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    rtt_compton_interface_dummy_data_test::compton_file_test(ut);
    rtt_compton_interface_dummy_data_test::compton_fail_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// End of test/tCSK_Interface_dummy_data.cc
//------------------------------------------------------------------------------------------------//
