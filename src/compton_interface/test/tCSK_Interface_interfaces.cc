//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   compton_interface/test/tCSK_Interface_interfaces.cc
 * \author Kendra Keady
 * \date   2017 Feb 10
 * \brief  Implementation file for tCSK_Interface_interfaces
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "compton_interface/CSK_Interface.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <array>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace rtt_compton_interface_interfaces_test {

using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

//!  Tests the CSK_Interface constructor and a couple of access routines.
void compton_file_test(rtt_dsxx::UnitTest &ut) {
  // Start the test.

  std::cout << "\n---------------------------------------------------------\n"
            << "   Test Draco code calling CSK_generator routines\n"
            << "---------------------------------------------------------\n";

  // open a small mg opacity file:
  const std::string filename = ut.getTestSourcePath() + "mg_ascii.compton";
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
  const std::vector<double> etemp_evals = compton_test->get_etemp_pts();

  FAIL_IF(grp_bds.size() != 2);
  FAIL_IF(etemp_evals.size() != 7);

  if (!soft_equiv(grp_bds[0], 19.9998996222))
    FAILMSG("Lower group bound read incorrectly!");
  if (!soft_equiv(grp_bds[1], 29.9998749832))
    FAILMSG("Upper group bound read incorrectly!");

  if (!soft_equiv(etemp_evals[0], 1.76377944e-05))
    FAILMSG("Etemp eval read incorrectly!");
  if (!soft_equiv(etemp_evals[6], 6.75507064e-04))
    FAILMSG("Etemp eval read incorrectly!");

  if (!soft_equiv(compton_test->get_min_etemp(), 1e-6))
    FAILMSG("Min etemp read incorrectly!");
  if (!soft_equiv(compton_test->get_max_etemp(), 7e-4))
    FAILMSG("Max etemp read incorrectly!");

  if (ut.numFails == 0) {
    std::cout << "\nCorrectly read group bounds and electron temps!" << std::endl;
  }

  // try "interpolating" at one of the exact eval points in the test library,
  // and check the result:
  const double test_etemp = 4.87227167e-04;
  std::vector<std::vector<std::vector<double>>> interp_data =
      compton_test->interpolate_csk(test_etemp)[0];
  // get interpolated nu_ratios
  std::vector<std::vector<double>> interp_nu_data = compton_test->interpolate_nu_ratio(test_etemp);

  // Check the size of the returned data:
  FAIL_IF(interp_nu_data.size() != 1);
  FAIL_IF(interp_nu_data[0].size() != 1);
  FAIL_IF(interp_data.size() != 1);
  FAIL_IF(interp_data[0].size() != 1);
  FAIL_IF(interp_data[0][0].size() != 4);

  // Check that the data is actually correct:
  if (!soft_equiv(interp_data[0][0][0], 4.45668383e+00))
    ITFAILS;
  if (!soft_equiv(interp_data[0][0][1], 3.17337784e-01))
    ITFAILS;
  if (!soft_equiv(interp_data[0][0][2], 4.50133379e-01))
    ITFAILS;
  if (!soft_equiv(interp_data[0][0][3], 3.59663442e-02))
    ITFAILS;

  if (!soft_equiv(interp_nu_data[0][0], 1.5000000e+00))
    ITFAILS;

  const double test_etemp2 = 6.75507064e-04;
  // get interpolated nu_ratio for a different
  std::vector<std::vector<double>> interp_nu_data2 =
      compton_test->interpolate_nu_ratio(test_etemp2);

  if (!soft_equiv(interp_nu_data2[0][0], 2.0000000e+00))
    ITFAILS;

  // get the number of xi evals in the library (we know it should be 4)
  if (compton_test->get_num_xi() != 4)
    ITFAILS;

  if (compton_test->get_num_groups() != 1)
    ITFAILS;

  if (ut.numFails == 0)
    std::cout << "\nCorrectly read multigroup data points!" << std::endl;

  if (ut.numFails == 0) {
    PASSMSG("Successfully linked Draco against CSK_generator.");
  } else {
    FAILMSG("Did not successfully link Draco against CSK_generator.");
  }
}

//------------------------------------------------------------------------------------------------//
void const_compton_file_test(rtt_dsxx::UnitTest &ut) {
  // Start the test.

  std::cout << "\n---------------------------------------------------------\n"
            << " Test Draco code calling CSK_generator routines -- const \n"
            << "---------------------------------------------------------\n";

  // open a small mg opacity file:
  const std::string filename = ut.getTestSourcePath() + "mg_ascii.compton";
  std::cout << "Attempting to construct a const CSK_Interface object...\n" << std::endl;
  std::unique_ptr<const rtt_compton_interface::CSK_Interface> compton_test;

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
  const std::vector<double> etemp_evals = compton_test->get_etemp_pts();

  FAIL_IF(grp_bds.size() != 2);
  FAIL_IF(etemp_evals.size() != 7);

  if (!soft_equiv(grp_bds[0], 19.9998996222))
    FAILMSG("Lower group bound read incorrectly!");
  if (!soft_equiv(grp_bds[1], 29.9998749832))
    FAILMSG("Upper group bound read incorrectly!");

  if (!soft_equiv(etemp_evals[0], 1.76377944e-05))
    FAILMSG("Etemp eval read incorrectly!");
  if (!soft_equiv(etemp_evals[6], 6.75507064e-04))
    FAILMSG("Etemp eval read incorrectly!");

  if (!soft_equiv(compton_test->get_min_etemp(), 1e-6))
    FAILMSG("Min etemp read incorrectly!");
  if (!soft_equiv(compton_test->get_max_etemp(), 7e-4))
    FAILMSG("Max etemp read incorrectly!");

  if (ut.numFails == 0) {
    std::cout << "\nCorrectly read group bounds and electron temps!" << std::endl;
  }

  // try "interpolating" at one of the exact eval points in the test library,
  // and check the result:
  const double test_etemp = 4.87227167e-04;
  // get interpolated csks
  std::vector<std::vector<std::vector<double>>> interp_data =
      compton_test->interpolate_csk(test_etemp)[0];
  // get interpolated nu_ratios
  std::vector<std::vector<double>> interp_nu_data = compton_test->interpolate_nu_ratio(test_etemp);

  // Check the size of the returned data:
  FAIL_IF(interp_nu_data.size() != 1);
  FAIL_IF(interp_nu_data[0].size() != 1);
  FAIL_IF(interp_data.size() != 1);
  FAIL_IF(interp_data[0].size() != 1);
  FAIL_IF(interp_data[0][0].size() != 4);

  // Check that the data is actually correct:
  if (!soft_equiv(interp_data[0][0][0], 4.45668383e+00))
    ITFAILS;
  if (!soft_equiv(interp_data[0][0][1], 3.17337784e-01))
    ITFAILS;
  if (!soft_equiv(interp_data[0][0][2], 4.50133379e-01))
    ITFAILS;
  if (!soft_equiv(interp_data[0][0][3], 3.59663442e-02))
    ITFAILS;

  if (!soft_equiv(interp_nu_data[0][0], 1.5000000e+00))
    ITFAILS;

  const double test_etemp2 = 6.75507064e-04;
  // get interpolated nu_ratio for a different
  std::vector<std::vector<double>> interp_nu_data2 =
      compton_test->interpolate_nu_ratio(test_etemp2);

  if (!soft_equiv(interp_nu_data2[0][0], 2.0000000e+00))
    ITFAILS;

  // get the number of xi evals in the library (we know it should be 4)
  if (compton_test->get_num_xi() != 4)
    ITFAILS;

  if (compton_test->get_num_groups() != 1)
    ITFAILS;

  if (ut.numFails == 0)
    std::cout << "\nCorrectly read multigroup data points!" << std::endl;

  if (ut.numFails == 0) {
    PASSMSG("Successfully linked Draco against CSK_generator (const version).");
  } else {
    FAILMSG("Did not successfully link Draco against CSK_generator (const "
            "version).");
  }
}

//------------------------------------------------------------------------------------------------//
//!  Tests the CSK_Interface mg build capability
void compton_build_test(rtt_dsxx::UnitTest &ut) {
  // Start the test.

  std::cout << "\n---------------------------------------------------------\n"
            << "Test Draco call to CSK_generator mg opacity builder\n"
            << "---------------------------------------------------------\n";

  // open a small pointwise opacity file:
  const std::string filename = ut.getTestSourcePath() + "lagrange_csk_ascii.compton";
  std::cout << "Attempting to construct a CSK_Interface object..." << std::endl;

  // make an uninitialized pointer...
  std::shared_ptr<rtt_compton_interface::CSK_Interface> compton_test;

  // make a small fake group structure to pass in:
  const std::vector<double> test_groups = {20.0, 30.0, 40.0, 50.0, 60.0};
  const std::string opac_type = "jayenne";
  const std::string wt_func = "planck";
  const bool induced = false;
  const bool det_bal = false;

  // set the number of angular points to retrieve (Legendre or otherwise)
  const size_t nxi = 5;

  try {
    // (This call has some output of its own, so we print some newlines around
    // it)
    std::cout << "\n\n";
    compton_test = std::make_unique<rtt_compton_interface::CSK_Interface>(
        filename, test_groups, opac_type, wt_func, induced, det_bal, nxi);
    std::cout << "\n\n";
  } catch (rtt_dsxx::assertion & /*asrt*/) {
    FAILMSG("Failed to construct a CSK_Interface object!");
    // if construction fails, there is no reason to continue testing...
    return;
  }
  std::cout << "(...Success!)" << std::endl;

  // Check some of the data in the CSK_generator-formed data file:
  const std::vector<double> grp_bds = compton_test->get_group_bounds();

  if (!soft_equiv(grp_bds.begin(), grp_bds.end(), test_groups.begin(), test_groups.end()))
    ITFAILS;

  if (!soft_equiv(compton_test->get_min_etemp(), 0.0))
    FAILMSG("Min etemp read incorrectly!");
  if (!soft_equiv(compton_test->get_max_etemp(), 0.0223607))
    FAILMSG("Max etemp read incorrectly!");

  if (ut.numFails == 0) {
    std::cout << "\nCorrectly stored group bounds and electron temps! " << std::endl;
  }

  // check that the files were constructed with the correct name, and
  // that they're open-able:
  std::string libfile = ut.getTestSourcePath() + "lagrange_csk_ascii.compton_4group";
  std::string nufile = ut.getTestSourcePath() + "lagrange_csk_ascii.compton_4group_nu_ratios";

  // "interpolate" the built MG data to check values:
  std::vector<std::vector<std::vector<std::vector<double>>>> csk_gen =
      compton_test->interpolate_csk(0.01);
  std::vector<std::vector<double>> nu_gen = compton_test->interpolate_nu_ratio(0.01);

  // open generated output files and perform same operation for comparison:
  compton_test = std::make_unique<rtt_compton_interface::CSK_Interface>(libfile);
  std::vector<std::vector<std::vector<std::vector<double>>>> csk_exist =
      compton_test->interpolate_csk(0.01);
  std::vector<std::vector<double>> nu_exist = compton_test->interpolate_nu_ratio(0.01);

  for (size_t a = 0; a < csk_gen.size(); a++) {
    for (size_t b = 0; b < csk_gen[a].size(); b++) {
      for (size_t c = 0; c < csk_gen[a][b].size(); c++) {
        if (!soft_equiv(csk_exist[a][b][c].begin(), csk_exist[a][b][c].end(),
                        csk_gen[a][b][c].begin(), csk_gen[a][b][c].end()))
          ITFAILS;
      }
    }
  }

  for (size_t a = 0; a < nu_gen.size(); a++) {
    if (!soft_equiv(nu_exist[a].begin(), nu_exist[a].end(), nu_gen[a].begin(), nu_gen[a].end()))
      ITFAILS;
  }

  if (!std::ifstream(libfile).good()) {
    ITFAILS;
  } else { // remove the file-- we don't actually want it sitting around...
    std::remove(libfile.c_str());
    std::cout << "Successfully wrote MG file -- deleting! " << std::endl;
  }
  if (!std::ifstream(nufile).good()) {
    ITFAILS;
  } else { // remove the file-- we don't actually want it sitting around...
    std::remove(nufile.c_str());
    std::cout << "Successfully wrote MG nu_ratio file -- deleting! " << std::endl;
  }

  if (ut.numFails == 0) {
    PASSMSG("Successfully built a CSK_generator mg library.");
  } else {
    FAILMSG("Did not successfully build a CSK_generator mg library.");
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

//------------------------------------------------------------------------------------------------//
//!  Tests Comptohn interface to LLNL-style CSK_Interface data.
void llnl_compton_test(rtt_dsxx::UnitTest &ut) {
  // Start the test.s
  std::cout << "\n---------------------------------------------------------\n"
            << " Test Draco calling CSK_generator LLNL-style routines \n"
            << "---------------------------------------------------------\n";

  const std::string filename = ut.getTestSourcePath() + "llnl_ascii.compton";
  const bool llnl_style = true;
  std::cout << "Attempting to construct a CSK_Interface object...\n" << std::endl;
  std::unique_ptr<rtt_compton_interface::CSK_Interface> compton_test;
  try {
    compton_test = std::make_unique<rtt_compton_interface::CSK_Interface>(filename, llnl_style);
  } catch (int /*asrt*/) {
    FAILMSG("Failed to construct an LLNL-style CSK_Interface object!");
    // if construction fails, there is no reason to continue testing...
    return;
  }
  std::cout << "\n(...Success!)" << std::endl;

  // Test the two sets of interpolators, and make sure they match.

  // Use temperatures that are very close to the library eval points:
  const std::vector<double> cell_temps = {1.50001, 2.49999};
  // Ensure densities are unit; these are only folded into the evaluations when
  // the data is pre-interpolated (not in the on-the-fly case)
  const std::vector<double> cell_dens = {1.0, 1.0};

  compton_test->interpolate_precycle(cell_temps, cell_dens);

  const std::array<double, 2> test_freq = {12.4233, 183.43};

  // reference sol'n for erec:
  const std::vector<double> ref_sigc = {1.91083e-01, 1.25551e-01, 1.91042e-01, 1.25395e-01};
  const std::vector<double> ref_erec = {-1.22843e-02, -2.00713e-01, -5.00821e-03, -1.974413e-01};

  double otf_sigcval;
  double otf_erecval;
  double pre_sigcval;
  double pre_erecval;

  size_t i = 0;
  for (size_t k = 0; k < cell_temps.size(); k++) {
    for (const double freq : test_freq) {
      // sigma_compton value:
      otf_sigcval = compton_test->interpolate_sigc(cell_temps[k], freq);
      pre_sigcval = compton_test->interpolate_cell_sigc(static_cast<int64_t>(k), freq);

      // expected relative energy change value:
      otf_erecval = compton_test->interpolate_erec(cell_temps[k], freq);
      pre_erecval = compton_test->interpolate_cell_erec(static_cast<int64_t>(k), freq);

      // compare the values to each other for consistency:
      FAIL_IF_NOT(soft_equiv(otf_sigcval, pre_sigcval));
      FAIL_IF_NOT(soft_equiv(otf_erecval, pre_erecval));

      // compare the values to the expected answer for accuracy:
      // (use a loose tolerance, because our points are close -- but not
      // equal to -- the evaluation points in the library)
      FAIL_IF_NOT(soft_equiv(otf_erecval, ref_erec[i], 1.0e-4));
      FAIL_IF_NOT(soft_equiv(otf_sigcval, ref_sigc[i], 1.0e-4));

      i++;
    }
  }

  if (ut.numFails == 0) {
    PASSMSG("Successfully read an LLNL-style CSK library.");
  } else {
    FAILMSG("Did not successfully read an LLNL-style CSK library.");
  }
}
} // namespace rtt_compton_interface_interfaces_test

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    rtt_compton_interface_interfaces_test::compton_file_test(ut);
    rtt_compton_interface_interfaces_test::const_compton_file_test(ut);
    rtt_compton_interface_interfaces_test::compton_build_test(ut);
    rtt_compton_interface_interfaces_test::compton_fail_test(ut);
    rtt_compton_interface_interfaces_test::llnl_compton_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// End of test/tCSK_Interface_interfaces.cc
//------------------------------------------------------------------------------------------------//
