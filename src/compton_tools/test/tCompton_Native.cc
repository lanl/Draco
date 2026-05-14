//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   compton_tools/test/tCompton_Native.cc
 * \author Andrew Till
 * \brief  Implementation file for tCompton_Native
 * \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Assert.hh"
#include "dsxx/Release.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "c4/ParallelUnitTest.hh"
#include "compton_tools/Compton_Native.hh"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// cSpell:ignore bds, evals, mec

namespace {

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
    std::ranges::for_each(variable, print);
    std::cout << '\n';
  };
  auto print_diff_lambda = [](const std::vector<double> &val, const std::vector<double> &ref,
                              const std::string &val_name, const std::string &ref_name,
                              double tol) {
    std::cout << std::setprecision(14);
    std::cout << "\n" << val_name << "/" << ref_name << " - 1\n";
    for (size_t i = 0; i < val.size(); ++i)
      std::cout << (val[i] - ref[i]) / (std::fabs(ref[i]) + tol) << ", ";
    std::cout << '\n';
  };

  // Tolerance used for checks
  const double tol = 1e-11;

  // Start the test.

  std::cout << "\n---------------------------------------------------------\n"
            << "             Test Draco Compton_Native data routines\n"
            << "---------------------------------------------------------\n";

  // open a small mg opacity file:
  const std::string filename = ut.getTestSourcePath() + "dummy_data_gold_b";
  std::cout << "Attempting to construct a Compton_Native object...\n\n";
  std::unique_ptr<rtt_compton_tools::Compton_Native> compton_test;

  try {
    compton_test = std::make_unique<rtt_compton_tools::Compton_Native>(filename);
  } catch (int /*assertion*/) {
    FAILMSG("Failed to construct a Compton_Native object!");
    // if construction fails, there is no reason to continue testing...
    return;
  }
  std::cout << "(...Success!)\n\n";

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
    std::cout << "Correctly read sizes, group bounds, and electron temps!\n\n";
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
    std::vector<double> outscat_matrix(G * G, -1.0);
    std::vector<double> outscat_matrix_diff(G * G, 0.0);
    std::vector<double> amp_matrix(G * G, -1.0);
    std::vector<double> nl_amp_matrix(G * G, -1.0);

    // Returns flattened 3D inscat array with order [moment, group_to, group_from]
    compton_test->interp_dense_inscat(inscat, interp_T_keV, G);

    // Returns 1D outscat array [group_from]
    compton_test->interp_linear_outscat(outscat, interp_T_keV);

    // Returns the flattened 2D outscattering matrix [incident_group, outbound_group]
    compton_test->interp_linear_outscat_matrix(outscat_matrix, interp_T_keV);

    // Returns the flattened 2D Compton amplification (nu_prime/nu) [incident_group, outbound_group]
    compton_test->interp_linear_compton_amplification(amp_matrix, interp_T_keV);

    // Use dummy flux for nonlinear component
    const double phival = 1.732984; // random // NOLINT(modernize-use-std-numbers)
    const double phiscale = phival * double(G);
    std::vector<double> const phi(compton_test->get_num_groups(), phival);
    // Returns 1D nonlinear difference array [group_from] (performs matrix multiplication with phi)
    compton_test->interp_nonlin_diff_and_add(nl_diff, interp_T_keV, phi, phiscale);

    // Update the linear outscattering matrix to include induced effects
    // (e.g. if a matrix of zeroes is passed in, on return it will be equal to the induced contribution)
    compton_test->interp_nonlin_diff_for_outscat_matrix(outscat_matrix_diff, interp_T_keV, phi,
                                                        phiscale);

    // Calculate the non-linear form of the Compton amplification factor
    compton_test->interp_linear_compton_amplification_nonlinear(nl_amp_matrix, interp_T_keV, phi,
                                                                phiscale);

    // Print result (useful if golds need updating)
    if (do_print) {
      print_lambda(inscat, "inscat");
      print_lambda(outscat, "outscat");
      print_lambda(nl_diff, "nl_diff");
      print_lambda(outscat_matrix, "outscat_matrix");
      print_lambda(outscat_matrix_diff, "outscat_matrix_diff");
      print_lambda(amp_matrix, "amp_matrix");
      print_lambda(nl_amp_matrix, "nl_amp_matrix");
    }

    // in_lin last temperature
    std::vector<double> const inscat_gold = {0.077335961983675,
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

    std::vector<double> const outscat_gold = {0.087476270446994, 0.074121503186901,
                                              0.064422154448942, 0.055814819950424};

    std::vector<double> const nl_diff_gold = {0.017286371393823, 0.0087938222181754,
                                              0.0011691552185996, 0.00012086891440242};
    std::vector<double> const outscat_matrix_gold = {0.076895615348207,
                                                     0.010580654813084,
                                                     2.8570204118996e-10,
                                                     0,
                                                     0.0010281421965568,
                                                     0.064450811855277,
                                                     0.0086424664123312,
                                                     8.2722735959847e-08,
                                                     1.0034858631934e-11,
                                                     0.002266205448414,
                                                     0.055552618936782,
                                                     0.0066033300537106,
                                                     0,
                                                     9.1038258953769e-09,
                                                     0.0020737244638516,
                                                     0.053741086382746};

    std::vector<double> const outscat_matrix_diff_gold = {0.023851451092196,
                                                          -0.0065650795529254,
                                                          -1.4544764008884e-10,
                                                          0,
                                                          0.0083242188042929,
                                                          0.0013414766097549,
                                                          -0.00087186671535155,
                                                          -6.4805208357506e-09,
                                                          3.8184464112549e-10,
                                                          0.001103724452621,
                                                          0.00015302427784575,
                                                          -8.7593893711727e-05,
                                                          0,
                                                          1.4153487374215e-08,
                                                          0.00011066101344709,
                                                          1.0193747467961e-05};

    std::vector<double> const amp_matrix_gold = {1.0057265506424,
                                                 1.1945583903626,
                                                 2.6511361177178,
                                                 0,
                                                 0.84437262368994,
                                                 1.0041096193111,
                                                 1.1915520468236,
                                                 2.1603584152337,
                                                 0.37458966030913,
                                                 0.84398336720346,
                                                 1.0022546490882,
                                                 1.1875224862621,
                                                 0,
                                                 0.46442455138961,
                                                 0.84429672767933,
                                                 1.0005625196344};

    std::vector<double> const nl_amp_matrix_gold = {0.76762495143641,
                                                    3.147546531602,
                                                    5.4004365949151,
                                                    0,
                                                    0.092825236745471,
                                                    0.98363625382467,
                                                    1.3252450190196,
                                                    2.3439869693074,
                                                    0.0095921177007557,
                                                    0.56756067968651,
                                                    0.99950143980737,
                                                    1.2034868709283,
                                                    0,
                                                    0.18179401070018,
                                                    0.80152463799734,
                                                    1.00037276635};

    // Print diff (useful if golds need updating)
    if (do_print) {
      print_diff_lambda(inscat, inscat_gold, "inscat", "inscat_gold", tol);
      print_diff_lambda(outscat, outscat_gold, "outscat", "outscat_gold", tol);
      print_diff_lambda(nl_diff, nl_diff_gold, "nl_diff", "nl_diff_gold", tol);
      print_diff_lambda(outscat_matrix, outscat_matrix_gold, "outscat_matrix",
                        "outscat_matrix_gold", tol);
      print_diff_lambda(outscat_matrix_diff, outscat_matrix_diff_gold, "outscat_matrix_diff",
                        "outscat_matrix_diff_gold", tol);
      print_diff_lambda(amp_matrix, amp_matrix_gold, "amp_matrix", "amp_matrix_gold", tol);
      print_diff_lambda(nl_amp_matrix, nl_amp_matrix_gold, "nl_amp_matrix", "nl_amp_matrix_gold",
                        tol);
    }

    ut.check(soft_equiv(inscat, inscat_gold, tol), "checked data retrieval for inscat");
    ut.check(soft_equiv(outscat, outscat_gold, tol), "checked data retrieval for outscat");
    ut.check(soft_equiv(nl_diff, nl_diff_gold, tol), "checked data retrieval for nl_diff");
    ut.check(soft_equiv(outscat_matrix, outscat_matrix_gold, tol), "checked data retrieval for "
                                                                   "outscat_matrix");
    ut.check(soft_equiv(outscat_matrix_diff, outscat_matrix_diff_gold, tol),
             "checked data retrieval for outscat_matrix_diff");
    ut.check(soft_equiv(amp_matrix, amp_matrix_gold, tol), "checked data retrieval for amp_matrix");
    ut.check(soft_equiv(nl_amp_matrix, nl_amp_matrix_gold, tol),
             "checked data retrieval for nl_amp_matrix");
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
    std::vector<double> outscat_matrix(G * G, -1.0);
    std::vector<double> outscat_matrix_diff(G * G, 0.0);
    std::vector<double> amp_matrix(G * G, -1.0);
    std::vector<double> nl_amp_matrix(G * G, -1.0);

    // Returns flattened 3D inscat array with order [moment, group_to, group_from]
    compton_test->interp_dense_inscat(inscat, interp_T_keV, G);

    // Returns 1D outscat array [group_from]
    compton_test->interp_linear_outscat(outscat, interp_T_keV);

    // Returns the flattend 2D outscattering matrix [incident_group, outbound_group]
    compton_test->interp_linear_outscat_matrix(outscat_matrix, interp_T_keV);

    // Returns the flattend 2D Compton amplification (nu_prim/nu) [incident_group, outbound_group]
    compton_test->interp_linear_compton_amplification(amp_matrix, interp_T_keV);

    // Use dummy flux for nonlinear component
    const double phival = 2.394559; // random
    const double phiscale = phival * double(G);
    std::vector<double> const phi(compton_test->get_num_groups(), phival);
    // Returns 1D nonlinear difference array [group_from] (performs matrix multiplication with phi)
    compton_test->interp_nonlin_diff_and_add(nl_diff, interp_T_keV, phi, phiscale);

    // Updates the outscattering matrix to include induced effects (accumulate to a zero array to
    // see just the induced contribution, typically we would pass in the linear outscat_matrix and
    // this would append the difference to that).
    compton_test->interp_nonlin_diff_for_outscat_matrix(outscat_matrix_diff, interp_T_keV, phi,
                                                        phiscale);

    // Calculate the non-linear form of the Compton amplification factor
    compton_test->interp_linear_compton_amplification_nonlinear(nl_amp_matrix, interp_T_keV, phi,
                                                                phiscale);

    // Print result (useful if golds need updating)
    if (do_print) {
      print_lambda(inscat, "inscat");
      print_lambda(outscat, "outscat");
      print_lambda(nl_diff, "nl_diff");
      print_lambda(outscat_matrix, "outscat_matrix");
      print_lambda(outscat_matrix_diff, "outscat_matrix_diff");
      print_lambda(amp_matrix, "amp_matrix");
      print_lambda(nl_amp_matrix, "nl_amp_matrix");
    }

    std::vector<double> const inscat_gold = {0.11237635504081,
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

    std::vector<double> const outscat_gold = {0.11533649452761, 0.090894893695562,
                                              0.076437945828947, 0.061604582646344};

    std::vector<double> const nl_diff_gold = {3.7130057760837e-05, 3.0320968944752e-05,
                                              4.7840481582021e-06, 4.108948064536e-07};

    std::vector<double> const outscat_matrix_gold = {
        0.11232257652334,  0.0030139180042694,  0, 0, 0.0040602444028484,
        0.085612046787231, 0.0012226025054832,  0, 0, 0.0042451435441658,
        0.071857006651994, 0.00033579563278701, 0, 0, 0.0087471769906374,
        0.05285740565570};

    std::vector<double> const outscat_matrix_diff_gold = {
        6.3095694580059e-05, -2.5965636819222e-05, 0, 0, 2.9066470259542e-05,
        4.6825132824432e-06, -3.428014597233e-06,  0, 0, 3.8895862938046e-06,
        1.1360439980447e-06, -2.415821336472e-07,  0, 0, 2.7283569671083e-07,
        1.3805910974278e-07};

    std::vector<double> const amp_matrix_gold = {
        1.0004787863592,  1.0774798785873, 0, 0, 0.95863446777599,
        0.99947559690179, 1.0778433221669, 0, 0, 0.92300239089228,
        0.99580955767251, 1.0618812278444, 0, 0, 0.91990231831319,
        0.96042716780069};

    std::vector<double> const nl_amp_matrix_gold = {
        0.99991709639734, 1.0868432980019, 0, 0, 0.951820576362,
        0.99942093401293, 1.0808739487668, 0, 0, 0.92215746986236,
        0.99579381438402, 1.0626457292182, 0, 0, 0.9198736262721,
        0.96042465925205};

    // Print diff (useful if golds need updating)
    if (do_print) {
      print_diff_lambda(inscat, inscat_gold, "inscat", "inscat_gold", tol);
      print_diff_lambda(outscat, outscat_gold, "outscat", "outscat_gold", tol);
      print_diff_lambda(nl_diff, nl_diff_gold, "nl_diff", "nl_diff_gold", tol);
      print_diff_lambda(outscat_matrix, outscat_matrix_gold, "outscat_matrix",
                        "outscat_matrix_gold", tol);
      print_diff_lambda(outscat_matrix_diff, outscat_matrix_diff_gold, "outscat_matrix_diff",
                        "outscat_matrix_diff_gold", tol);
      print_diff_lambda(amp_matrix, amp_matrix_gold, "amp_matrix", "amp_matrix_gold", tol);
      print_diff_lambda(nl_amp_matrix, nl_amp_matrix_gold, "nl_amp_matrix", "nl_amp_matrix_gold",
                        tol);
    }

    ut.check(soft_equiv(inscat, inscat_gold, tol), "checked data retrieval for inscat");
    ut.check(soft_equiv(outscat, outscat_gold, tol), "checked data retrieval for outscat");
    ut.check(soft_equiv(nl_diff, nl_diff_gold, tol), "checked data retrieval for nl_diff");
    ut.check(soft_equiv(outscat_matrix, outscat_matrix_gold, tol),
             "checked data retrieval for outscat_matrix");
    ut.check(soft_equiv(outscat_matrix_diff, outscat_matrix_diff_gold, tol),
             "checked data retrieval for outscat_matrix_diff");
    ut.check(soft_equiv(amp_matrix, amp_matrix_gold, tol), "checked data retrieval for amp_matrix");
    ut.check(soft_equiv(nl_amp_matrix, nl_amp_matrix_gold, tol),
             "checked data retrieval for nl_amp_matrix");
  }
}

//------------------------------------------------------------------------------------------------//
//!  Tests Compton's error-handling on a non-existent file.
void bad_file_test(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n---------------------------------------------------------\n"
            << "    Test Compton_Native bad file handling    \n"
            << "---------------------------------------------------------\n";
  // open a small mg opacity file:
  std::string const filename = ut.getTestSourcePath() + "non_existent_b";
  std::cout << "Testing with a non-existent file...\n\n";
  std::unique_ptr<rtt_compton_tools::Compton_Native> compton_test;

  bool caught = false;
  try {
    compton_test = std::make_unique<rtt_compton_tools::Compton_Native>(filename);
  } catch (rtt_dsxx::assertion &assertion) {
    std::cout << "Draco exception thrown: " << assertion.what() << '\n';
    // We successfully caught the bad file!
    caught = true;
  } catch (const std::ifstream::failure & /* error */) {
    std::cout << "ifstream failure caught!\n";
    caught = true;
  } catch (...) {
    std::cout << "Unidentified exception caught!\n";
    caught = true;
  }

  FAIL_IF_NOT(caught);

  if (ut.numFails == 0) {
    PASSMSG("Successfully caught a bad file exception.");
  } else {
    FAILMSG("Did not successfully catch a bad file exception.");
  }
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    test(ut);
    bad_file_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// End of test/tCompton_Native.cc
//------------------------------------------------------------------------------------------------//
