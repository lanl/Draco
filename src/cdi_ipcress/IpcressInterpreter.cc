//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ipcress/IpcressInterpreter.cc
 * \author Allan Wollaber
 * \date   Fri Oct 12 15:39:39 2001
 * \brief  Basic reader to print info in IPCRESS files.
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "IpcressFile.hh"
#include "IpcressGrayOpacity.hh"
#include "IpcressMultigroupOpacity.hh"
#include "cdi/OpacityCommon.hh"
#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "ds++/XGetopt.hh"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

using rtt_cdi_ipcress::IpcressFile;
using rtt_cdi_ipcress::IpcressGrayOpacity;
using rtt_cdi_ipcress::IpcressMultigroupOpacity;
using std::cout;
using std::endl;
using std::ios;
using std::string;

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Basic reader to print info in IPCRESS files.
 *
 * Information about the IPCRESS IPCRESS Interface Library may be found on X-Division's Physical
 * Data Team web site:
 *
 * http://velvet.lanl.gov/PROJECTS/DATA/atomic/ipcress/intro.lasso
 *
 * We have a slightly modified copy of the ipcress libraries (explicitly added compiler libraries)
 * located at:
 *
 * /radtran/vendors/ipcress
 *
 * To use this package Draco must be compiled with Ipcress
 *
 * Modification of this executable to make it more useful is encouraged.
 */
//------------------------------------------------------------------------------------------------//
void ipcress_file_read(std::string const &op_data_file) {
  // Ipcress data filename (IPCRESS format required)

  std::shared_ptr<IpcressFile> spGF;
  try {
    spGF = std::make_shared<rtt_cdi_ipcress::IpcressFile>(op_data_file);
  } catch (rtt_dsxx::assertion const & /*excpt*/) {
    std::cerr << "Error: Can't open file " << op_data_file << ". Aborting" << endl;
    throw;
  }

  // Test the new object to verify the constructor and accessors.
  size_t numMaterials = spGF->getNumMaterials();
  cout << "This opacity file has " << numMaterials << " materials:" << endl;

  // Print the Mat IDs.
  std::vector<size_t> matIDs = spGF->getMatIDs();
  for (size_t i = 0; i < numMaterials; ++i)
    cout << "Material " << i + 1 << " has ID number " << matIDs[i] << endl;

  size_t i = 0;
  cout << "=========================================\n"
       << "For Material " << i + 1 << " with ID number " << matIDs[i]
       << "\n=========================================" << endl;
  auto spMgOp = std::make_shared<IpcressMultigroupOpacity>(spGF, matIDs[i], rtt_cdi::ROSSELAND,
                                                           rtt_cdi::TOTAL);

  // set precision
  cout.precision(7);
  cout.setf(ios::scientific);

#if defined(MSVC) && MSVC_VERSION < 1900
  // [2015-02-06 KT]: By default, MSVC uses a 3-digit exponent (presumably because
  // numeric_limits<double>::max() has a 3-digit exponent.)  Enable two-digit exponent format to
  // stay consistent with GNU and Intel on Linux.(requires <stdio.h>).
  unsigned old_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

  // Print the density grid
  std::vector<double> dens = spMgOp->getDensityGrid();
  cout << "Density grid\n";
  for (size_t tIndex = 0; tIndex < dens.size(); ++tIndex)
    cout << tIndex + 1 << "\t" << dens[tIndex] << endl;

  // Print the temperature grid
  std::vector<double> temps = spMgOp->getTemperatureGrid();
  cout << "Temperature grid\n";
  for (size_t tIndex = 0; tIndex < temps.size(); ++tIndex)
    cout << tIndex + 1 << "\t" << temps[tIndex] << endl;

  // Print the frequency grid
  std::vector<double> groups = spMgOp->getGroupBoundaries();
  cout << "Frequency grid\n";
  for (size_t gIndex = 0; gIndex < groups.size(); ++gIndex)
    cout << gIndex + 1 << "\t" << groups[gIndex] << endl;
  cout << endl;

  int keepGoing = 1;
  size_t matID(1);
  size_t selID(1);
  while (keepGoing) {
    cout << "Enter 0 or q to quit." << endl;
    cout << "Please select a material from 1 to " << numMaterials << endl;
    std::cin >> keepGoing;
    if (keepGoing == 0)
      break;
    matID = keepGoing - 1;
    Insist((matID < numMaterials), "Invalid material index");

    double density(0.0);
    cout << "Please enter a density from " << dens[0] << " to " << dens[dens.size() - 1] << endl;
    std::cin >> density;
    if (density <= 0.0)
      break;

    double temperature(0.0);
    cout << "Please select a temperature from " << temps[0] << " to " << temps[temps.size() - 1]
         << endl;
    std::cin >> temperature;
    if (temperature <= 0.0)
      break;

    cout << "Choose your opacity type:\n"
         << "1: Rosseland Absorption (Gray), 2: Planck Absorption (Gray)\n"
         << "3: Rosseland Absorption (MG),   4: Planck Absorption (MG)\n"
         << "5: Rosseland Total (Gray)" << endl;
    std::cin >> keepGoing;
    if (keepGoing == 0)
      break;
    selID = keepGoing;
    Insist((selID < 6), "Invalid choice.");

    if (selID == 1) {
      auto spGOp = std::make_shared<IpcressGrayOpacity>(spGF, matIDs[matID], rtt_cdi::ROSSELAND,
                                                        rtt_cdi::ABSORPTION);
      cout << "The Gray Rosseland Absorption Opacity for\n"
           << "material " << matID << " Id(" << matIDs[matID] << ") at density " << density
           << ", temperature " << temperature << " is " << spGOp->getOpacity(temperature, density)
           << endl;
    } else if (selID == 2) {
      auto spGOp = std::make_shared<IpcressGrayOpacity>(spGF, matIDs[matID], rtt_cdi::PLANCK,
                                                        rtt_cdi::ABSORPTION);
      cout << "The Gray Planck Absorption Opacity for\n"
           << "material " << matID << " Id(" << matIDs[matID] << ") at density " << density
           << ", temperature " << temperature << " is " << spGOp->getOpacity(temperature, density)
           << endl;

    } else if (selID == 3) {
      auto spMGOp = std::make_shared<IpcressMultigroupOpacity>(
          spGF, matIDs[matID], rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION);
      cout << "The Multigroup Rosseland Absorption Opacity for\n"
           << "material " << matID << " Id(" << matIDs[matID] << ") at density " << density
           << ", temperature " << temperature << " is: " << endl;
      std::vector<double> opData = spMGOp->getOpacity(temperature, density);
      cout << "Index \t Group Center \t\t Opacity" << endl;
      for (size_t g = 0; g < opData.size(); ++g)
        cout << g + 1 << "\t " << 0.5 * (groups[g] + groups[g + 1]) << "   \t " << opData[g]
             << endl;
    } else if (selID == 4) {
      auto spMGOp = std::make_shared<IpcressMultigroupOpacity>(spGF, matIDs[matID], rtt_cdi::PLANCK,
                                                               rtt_cdi::ABSORPTION);
      cout << "The Multigroup Planck Absorption Opacity for\n"
           << "material " << matID << " Id(" << matIDs[matID] << ") at density " << density
           << ", temperature " << temperature << " is: " << endl;
      std::vector<double> opData = spMGOp->getOpacity(temperature, density);
      cout << "Index \t Group Center  \t\t Opacity" << endl;
      for (size_t g = 0; g < opData.size(); ++g)
        cout << g + 1 << "\t " << 0.5 * (groups[g] + groups[g + 1]) << "   \t " << opData[g]
             << endl;
    }
    if (selID == 5) {
      auto spGOp = std::make_shared<IpcressGrayOpacity>(spGF, matIDs[matID], rtt_cdi::ROSSELAND,
                                                        rtt_cdi::TOTAL);
      cout << "The Gray Rosseland Total Opacity for\n"
           << "material " << matID << " Id(" << matIDs[matID] << ") at density " << density
           << ", temperature " << temperature << " is " << spGOp->getOpacity(temperature, density)
           << endl;
    }
  }

#if defined(MSVC) && MSVC_VERSION < 1900
  // Disable two-digit exponent format
  _set_output_format(old_exponent_format);
#endif

  cout << "Ending session." << endl;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  // Process known command line arguments:
  rtt_dsxx::XGetopt::csmap long_options;
  long_options['h'] = "help";
  long_options['v'] = "version";
  std::map<char, std::string> help_strings;
  help_strings['h'] = "print this message.";
  help_strings['v'] = "print version information and exit.";
  rtt_dsxx::XGetopt program_options(argc, argv, long_options, help_strings);

  std::string const helpstring("\nUsage: IpcressInterpreter [-hv] <ipcress file>\nFollow the "
                               "prompts to print opacity data to the screen.");

  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'v': { // --version
      cout << argv[0] << ": version " << rtt_dsxx::release() << endl;
      return 0;
    }
    case 'h': { // --help
      cout << argv[0] << ": version " << rtt_dsxx::release() << helpstring << endl;
      return 0;
    }
    default: {
      cout << argv[0] << ": version " << rtt_dsxx::release() << helpstring << endl;
      return 0;
    }
    }
  }

  // Assume last command line argument is the name of the ipcress file.
  std::string const filename = string(argv[argc - 1]);

  try {
    // >>> UNIT TESTS
    ipcress_file_read(filename);
  } catch (rtt_dsxx::assertion &excpt) {
    cout << "While attempting to read an opacity file, " << excpt.what() << endl;
    return 1;
  }

  return 0;
}

//------------------------------------------------------------------------------------------------//
// end of IpcressInterpreter.cc
//------------------------------------------------------------------------------------------------//
