//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   compton_tools/test/tCompton_Edep.cc
 * \author Kendra P. Long
 * \date   2022 Oct 7
 * \brief  Implementation file for tCompton_Edep
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "compton_tools/Compton_Edep.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace rtt_compton_edep_test {

using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

//!  Tests the CSK_Interface constructor and a couple of access routines.
void compton_edep_data_test(rtt_dsxx::UnitTest &ut) {

  // Start the test.
  std::cout << "\n---------------------------------------------------------\n"
            << "   Test Draco code reading / interping Edep data \n"
            << "---------------------------------------------------------\n";

  // open a small Edep library:
  const std::string filename = ut.getTestSourcePath() + "dummy_ascii";
  std::cout << "Attempting to construct an Edep_Container object...\n" << std::endl;
  std::shared_ptr<compton_edep::Edep_Container> edep_sp;

  try {
    edep_sp =
        std::make_shared<compton_edep::Edep_Container>(compton_edep::read_edep_file(filename));
  } catch (int /*asrt*/) {
    FAILMSG("Failed to construct an Edep_Container object!");
    // if construction fails, there is no reason to continue testing...
    return;
  }
  std::cout << "\n(...Success!)" << std::endl;

  // Use temperatures/frequencies that are close to the library eval points:
  const std::vector<double> cell_temps = {1.50001, 2.49999};
  const std::array<double, 3> test_freq = {12.4233, 183.43, 25.55};

  // Form the interpolated-in-temperature data from the pointwise data:
  compton_edep::Preinterp_Container preinterped_data(
      cell_temps, std::vector<double>(cell_temps.size(), 1.0), *edep_sp);

  // reference sol'n for sigma_c and EREC:
  const std::vector<double> ref_sigc = {1.91083e-01, 1.25551e-01, 1.82477e-01,
                                        1.91042e-01, 1.25395e-01, 1.82404e-01};
  const std::vector<double> ref_erec = {-1.22843e-02, -2.00713e-01, -3.51937e-02,
                                        -5.00821e-03, -1.97441e-01, -2.85429e-02};

  std::pair<double, double> pre_vals;

  for (size_t k = 0; k < cell_temps.size(); k++) {
    for (size_t i = 0; i < test_freq.size(); i++) {
      // get pair containing sigma_compton value and expected relative energy change value:
      pre_vals = compton_edep::interpolate_sigc_erec((k + 1), test_freq[i], preinterped_data);
      // compare the values to the expected answer for accuracy:
      // (use a loose tolerance, because our points are close -- but not
      // equal to -- the evaluation points in the library)
      FAIL_IF_NOT(soft_equiv(pre_vals.first, ref_sigc[3 * k + i], 1.0e-4));
      FAIL_IF_NOT(soft_equiv(pre_vals.second, ref_erec[3 * k + i], 1.0e-4));
    }
  }

  // Test a corner case where the temp is right on top of a library point:
  compton_edep::Preinterp_Container preinterped_onepoint_data({2.5}, {1.0}, *edep_sp);
  auto onepoint_vals = compton_edep::interpolate_sigc_erec(1, 400.0, preinterped_onepoint_data);
  // Include additional normalization factor on sigma.
  FAIL_IF_NOT(soft_equiv(onepoint_vals.first, (2.529432e+00 * 3.755728e-02), 1.0e-04));
  FAIL_IF_NOT(soft_equiv(onepoint_vals.second, -3.04825e-01, 1.0e-04));

  // Test a corner case where the frequency and temp are right on top of a library point:
  compton_edep::Preinterp_Container preinterped_twopoint_data({1.5}, {1.0}, *edep_sp);
  auto twopoint_vals = compton_edep::interpolate_sigc_erec(1, 459.68526, preinterped_twopoint_data);
  // Include additional normalization factor on sigma.
  FAIL_IF_NOT(soft_equiv(twopoint_vals.first, (2.39488e+00 * 3.755728e-02), 1.0e-4));
  FAIL_IF_NOT(soft_equiv(twopoint_vals.second, -3.26781061e-01, 1.0e-4));

  if (ut.numFails == 0)
    std::cout << "\nCorrectly interpolated all data points!!" << std::endl;

  if (ut.numFails == 0) {
    PASSMSG("Successfully read and interpolated CSK Edep data.");
  } else {
    FAILMSG("Did not successfully read and interpolate CSK Edep data.");
  }
}

//------------------------------------------------------------------------------------------------//
//!  Tests CSK_Interface's error-handling on a non-existent file.
void compton_fail_test(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n---------------------------------------------------------\n"
            << "    Test Compton_Edep bad file handling    \n"
            << "---------------------------------------------------------\n";
  // open a small mg opacity file:
  std::string filename = ut.getTestSourcePath() + "non_existent.compton";
  std::cout << "Testing with a non-existent file...\n" << std::endl;

  bool caught = false;
  try {
    compton_edep::Edep_Container mycontainer = compton_edep::read_edep_file(filename);
  } catch (rtt_dsxx::assertion &asrt) {
    std::cout << "Draco exception thrown: " << asrt.what() << std::endl;
    // We successfully caught the bad file!
    caught = true;
  }

  if (!caught)
    ITFAILS;

  if (ut.numFails == 0) {
    PASSMSG("Successfully caught a bad library file.");
  } else {
    FAILMSG("Did not successfully catch a bad library file.");
  }
}

} // namespace rtt_compton_edep_test

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    rtt_compton_edep_test::compton_edep_data_test(ut);
    rtt_compton_edep_test::compton_fail_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// End of test/tCompton_Edep.cc
//------------------------------------------------------------------------------------------------//
