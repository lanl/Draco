//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   compton_tools/test/tCompton_Native.cc
 * \author Andrew Till
 * \date   14 Oct 2020
 * \brief  Implementation file for tCompton_Native
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "compton_tools/Compton_Native.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace rtt_compton_tools_test {

using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

//!  Simple test
void test(rtt_dsxx::UnitTest &ut) {
  // Make true if golds need updating
  const bool do_print = false;
  auto print_lambda = [](const std::vector<double> &variable, const std::string &variable_name) {
    auto print = [](double a) { std::cout << a << ", "; };
    std::cout << std::setprecision(14);
    std::cout << "\n" << variable_name << "\n";
    std::for_each(variable.begin(), variable.end(), print);
    std::cout << std::endl;
  };
  auto print_diff_lambda = [](const std::vector<double> &val, const std::vector<double> &ref,
                              const std::string &val_name, const std::string &ref_name,
                              double tol) {
    std::cout << std::setprecision(14);
    std::cout << "\n" << val_name << "/" << ref_name << " - 1\n";
    for (size_t i = 0; i < val.size(); ++i)
      std::cout << (val[i] - ref[i]) / (std::fabs(ref[i]) + tol) << ", ";
    std::cout << std::endl;
  };

  // Tolerance used for checks
  const double tol = 1e-11;

  // Start the test.

  std::cout << "\n---------------------------------------------------------\n"
            << "             Test Draco Compton_Native data routines\n"
            << "---------------------------------------------------------\n";

  // open a small mg opacity file:
  const std::string filename = ut.getTestSourcePath() + "dummy_data_gold_b";
  std::cout << "Attempting to construct a Compton_Native object...\n" << std::endl;
  std::unique_ptr<rtt_compton_tools::Compton_Native> compton_test;

  try {
    compton_test = std::make_unique<rtt_compton_tools::Compton_Native>(filename);
  } catch (int /*asrt*/) {
    FAILMSG("Failed to construct a Compton_Native object!");
    // if construction fails, there is no reason to continue testing...
    return;
  }
  std::cout << "(...Success!)\n" << std::endl;

  // Check some of the data in the CSK_generator-opened file:
  const std::vector<double> grp_bds = compton_test->get_Egs();
  const std::vector<double> T_evals = compton_test->get_Ts();

  // Unitless (divided by mec2)
  // NB: These values can be read directly from the 3rd line in the ASCII data files
  std::vector<double> grp_bds_gold = {1.57311251e-06, 3.14622503e-04, 7.86556258e-04,
                                      1.57311251e-03, 3.14622503e-02};
  // NB: These values are scattered throughout the ASCII data files as headers
  std::vector<double> T_evals_gold = {1.57311251e-05, 1.57311251e-04, 3.30353629e-04,
                                      6.60707256e-04};

  // Sizes
  const size_t num_groups_gold = grp_bds_gold.size() - 1U;
  const size_t num_T_evals_gold = T_evals_gold.size();
  const size_t num_evals_gold = 3; // in_lin, out_lin, diff_nonlin
  const size_t num_leg_moments_gold = 2;
  // A point is a (Legendre moment, evaluation) pair
  // first eval (in_lin) has all Leg moments and all others have only the 0th moment
  const size_t num_points_gold = num_leg_moments_gold + (num_evals_gold - 1U);

  // Multiply by electron rest-mass energy (keV; using CSK value)
  const double mec2 = 510.998;
  for (double &grp_bd : grp_bds_gold) {
    grp_bd *= mec2;
  }
  for (double &T_eval : T_evals_gold) {
    T_eval *= mec2;
  }

  ut.check(grp_bds.size() == (num_groups_gold + 1U), "checked size of group bounds vector");
  ut.check(soft_equiv(grp_bds, grp_bds_gold, tol), "checked group boundaries");
  ut.check(soft_equiv(T_evals, T_evals_gold, tol), "checked temperature grid");

  // Test size accessor functions
  ut.check(compton_test->get_num_temperatures() == num_T_evals_gold,
           "checked number of temperatures");
  ut.check(compton_test->get_num_groups() == num_groups_gold, "checked number of groups");
  ut.check(compton_test->get_num_leg_moments() == num_leg_moments_gold,
           "checked number of Legendre moments");
  ut.check(compton_test->get_num_evals() == num_evals_gold, "checked number of evaluations");
  ut.check(compton_test->get_num_points() == num_points_gold, "checked number of points");
  ut.check(compton_test->get_highest_leg_moment() == (num_leg_moments_gold - 1U),
           "checked highest Legendre moment");

  if (ut.numFails == 0) {
    std::cout << "Correctly read sizes, group bounds, and electron temps!\n" << std::endl;
  }

  // Test data retrieval: interpolate to a grid point in temperature
  {
    const double interp_T_keV = T_evals[num_T_evals_gold - 1U];
    std::cout << "Testing interpolation at T = " << interp_T_keV << " keV\n";

    const size_t G = compton_test->get_num_groups();
    const size_t L = compton_test->get_num_leg_moments();
    std::vector<double> inscat(G * G * L, -1.0);
    std::vector<double> outscat(G, -1.0);
    std::vector<double> nl_diff(G, 0.0);

    // Returns flattened 3D inscat array with order [moment, group_to, group_from]
    compton_test->interp_dense_inscat(inscat, interp_T_keV, G);

    // Returns 1D outscat array [group_from]
    compton_test->interp_linear_outscat(outscat, interp_T_keV);

    // Use dummy flux for nonlinear component
    const double phival = 1.732984; // random
    const double phiscale = phival * double(G);
    std::vector<double> phi(compton_test->get_num_groups(), phival);
    // Returns 1D nonlinear difference array [group_from] (performs matrix multiplication with phi)
    compton_test->interp_nonlin_diff_and_add(nl_diff, interp_T_keV, phi, phiscale);

    // Print result (useful if golds need updating)
    if (do_print) {
      print_lambda(inscat, "inscat");
      print_lambda(outscat, "outscat");
      print_lambda(nl_diff, "nl_diff");
    }

    // in_lin last temperature
    std::vector<double> inscat_gold = {0.077335961983675,
                                       0.00086813512403298,
                                       3.7589542861865e-12,
                                       0,
                                       0.0126392099825,
                                       0.064715680156296,
                                       0.0019126397051273,
                                       4.2280402573895e-09,
                                       7.5743500030441e-10,
                                       0.010297948543218,
                                       0.055677870598416,
                                       0.0017508387789385,
                                       0,
                                       1.7871075876201e-07,
                                       0.0078416029229916,
                                       0.05377131679901,
                                       0.00025245247880701,
                                       -1.9130636257311e-05,
                                       -3.3342618791366e-13,
                                       0,
                                       -0.00033789734745753,
                                       0.00025068616810335,
                                       -4.2277780938278e-05,
                                       -3.556692022306e-10,
                                       -6.7240010226582e-11,
                                       -0.00026875796575325,
                                       0.00026147382401739,
                                       -3.8250062923392e-05,
                                       0,
                                       -1.5099779118058e-08,
                                       -0.00019672159570614,
                                       0.00018670354202682};

    std::vector<double> outscat_gold = {0.087476270446994, 0.074121503186901, 0.064422154448942,
                                        0.055814819950424};

    std::vector<double> nl_diff_gold = {0.017286371393823, 0.0087938222181754, 0.0011691552185996,
                                        0.00012086891440242};

    // Print diff (useful if golds need updating)
    if (do_print) {
      print_diff_lambda(inscat, inscat_gold, "inscat", "inscat_gold", tol);
      print_diff_lambda(outscat, outscat_gold, "outscat", "outscat_gold", tol);
      print_diff_lambda(nl_diff, nl_diff_gold, "nl_diff", "nl_diff_gold", tol);
    }

    ut.check(soft_equiv(inscat, inscat_gold, tol), "checked data retrieval for inscat");
    ut.check(soft_equiv(outscat, outscat_gold, tol), "checked data retrieval for outscat");
    ut.check(soft_equiv(nl_diff, nl_diff_gold, tol), "checked data retrieval for nl_diff");
  }
  // Test interpolation
  {
    const double alpha = 0.4;
    const double interp_T_keV = alpha * T_evals[0] + (1.0 - alpha) * T_evals[1];
    std::cout << "Testing interpolation at T = " << interp_T_keV << " keV\n";

    const size_t G = compton_test->get_num_groups();
    const size_t L = compton_test->get_num_leg_moments();
    std::vector<double> inscat(G * G * L, -1.0);
    std::vector<double> outscat(G, -1.0);
    std::vector<double> nl_diff(G, 0.0);

    // Returns flattened 3D inscat array with order [moment, group_to, group_from]
    compton_test->interp_dense_inscat(inscat, interp_T_keV, G);

    // Returns 1D outscat array [group_from]
    compton_test->interp_linear_outscat(outscat, interp_T_keV);

    // Use dummy flux for nonlinear component
    const double phival = 2.394559; // random
    const double phiscale = phival * double(G);
    std::vector<double> phi(compton_test->get_num_groups(), phival);
    // Returns 1D nonlinear difference array [group_from] (performs matrix multiplication with phi)
    compton_test->interp_nonlin_diff_and_add(nl_diff, interp_T_keV, phi, phiscale);

    // Print result (useful if golds need updating)
    if (do_print) {
      print_lambda(inscat, "inscat");
      print_lambda(outscat, "outscat");
      print_lambda(nl_diff, "nl_diff");
    }

    std::vector<double> inscat_gold = {0.11237635504081,
                                       0.003892290232165,
                                       0,
                                       0,
                                       0.0032474360053123,
                                       0.085567151564651,
                                       0.003918277640946,
                                       0,
                                       0,
                                       0.0013177739461996,
                                       0.071555894009793,
                                       0.0080465483923832,
                                       0,
                                       0,
                                       0.00035657507884865,
                                       0.050765688411202,
                                       9.3287262136785e-05,
                                       -9.2936704368115e-05,
                                       0,
                                       0,
                                       -8.1502684127148e-05,
                                       0.00016450195164588,
                                       -9.3405055656181e-05,
                                       0,
                                       0,
                                       -3.1620447251114e-05,
                                       0.00018889454989741,
                                       -0.00017952275276575,
                                       0,
                                       0,
                                       -4.6914663077267e-06,
                                       0.0004902506924202};

    std::vector<double> outscat_gold = {0.11533649452761, 0.090894893695562, 0.076437945828947,
                                        0.061604582646344};

    std::vector<double> nl_diff_gold = {3.7130057760837e-05, 3.0320968944752e-05,
                                        4.7840481582021e-06, 4.108948064536e-07};

    // Print diff (useful if golds need updating)
    if (do_print) {
      print_diff_lambda(inscat, inscat_gold, "inscat", "inscat_gold", tol);
      print_diff_lambda(outscat, outscat_gold, "outscat", "outscat_gold", tol);
      print_diff_lambda(nl_diff, nl_diff_gold, "nl_diff", "nl_diff_gold", tol);
    }

    ut.check(soft_equiv(inscat, inscat_gold, tol), "checked data retrieval for inscat");
    ut.check(soft_equiv(outscat, outscat_gold, tol), "checked data retrieval for outscat");
    ut.check(soft_equiv(nl_diff, nl_diff_gold, tol), "checked data retrieval for nl_diff");
  }
}

//------------------------------------------------------------------------------------------------//
//!  Tests Compton's error-handling on a non-existent file.
void bad_file_test(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n---------------------------------------------------------\n"
            << "    Test Compton_Native bad file handling    \n"
            << "---------------------------------------------------------\n";
  // open a small mg opacity file:
  std::string filename = ut.getTestSourcePath() + "non_existent_b";
  std::cout << "Testing with a non-existent file...\n" << std::endl;
  std::unique_ptr<rtt_compton_tools::Compton_Native> compton_test;

  bool caught = false;
  try {
    compton_test = std::make_unique<rtt_compton_tools::Compton_Native>(filename);
  } catch (rtt_dsxx::assertion &asrt) {
    std::cout << "Draco exception thrown: " << asrt.what() << std::endl;
    // We successfully caught the bad file!
    caught = true;
  } catch (const std::ifstream::failure & /* error */) {
    std::cout << "ifstream failure caught!" << std::endl;
    caught = true;
  } catch (...) {
    std::cout << "Unidentified exception caught!" << std::endl;
    caught = true;
  }

  if (!caught)
    ITFAILS;

  if (ut.numFails == 0) {
    PASSMSG("Successfully caught a bad file exception.");
  } else {
    FAILMSG("Did not successfully catch a bad file exception.");
  }
}
} // namespace rtt_compton_tools_test

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    rtt_compton_tools_test::test(ut);
    rtt_compton_tools_test::bad_file_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// End of test/tCompton_Native.cc
//------------------------------------------------------------------------------------------------//
