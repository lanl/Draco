//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_eospac/QueryEospac.cc
 * \author Kelly Thompson
 * \date   Friday, Nov 09, 2012, 13:02 pm
 * \brief  An interactive program for querying data from EOSPAC.
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Eospac.hh"
#include "EospacException.hh"
#include "SesameTables.hh"
#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "ds++/XGetopt.hh"
#include <cmath>
#include <cstdlib> // defines atoi
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

//------------------------------------------------------------------------------------------------//
void query_eospac() {
  using std::cout;
  using std::endl;

  cout << "Starting QueryEospac...\n\nMaterial tables are published at "
       << "http://xweb.lanl.gov/projects/data.\n"
       << endl;

  bool keepGoing(true);
  while (keepGoing) {

    // Table ID
    std::string entry;
    cout << "What material/table id (0 or q to quit)? ";
    std::cin >> entry;
    if (entry == "q" || entry == "0") {
      std::cout << std::endl;
      break;
    }
    // Convert 'string entry' to 'int tableID'
    int tableID = atoi(entry.c_str());

    // Create a SesameTable
    rtt_cdi_eospac::SesameTables SesameTab;

    // Query user for table(s) to query.
    std::string eosprop("");
    cout << "What property (h for help)? ";
    std::cin >> eosprop;
    if (eosprop == std::string("h")) {
      //printeosproplist();
      SesameTab.printEosTableList();
      cout << "What property (h for help)? ";
      std::cin >> eosprop;
    }

    // Register some EOS tables...
    if (eosprop == std::string("Uic_DT") || eosprop == std::string("EOS_Uic_DT"))
      SesameTab.Uic_DT(tableID);
    else if (eosprop == std::string("Ktc_DT") || eosprop == std::string("EOS_Ktc_DT"))
      SesameTab.Ktc_DT(tableID);
    else {
      cout << "Requested EOS property unknown or currently unsupported.  "
           << "Please select one of (you can ommit the prefix EOS_):" << endl;
      SesameTab.printEosTableList();
      continue;
    }

    // Generate EOS Table
    std::shared_ptr<rtt_cdi_eospac::Eospac const> spEospac(new rtt_cdi_eospac::Eospac(SesameTab));

    // Parameters
    double temp(0.0);
    double dens(0.0);
    cout << "Evaluate at\n  Temperature (keV): ";
    std::cin >> temp;
    cout << "  Density (g/cm^3): ";
    std::cin >> dens;

    // Result
    cout << "For table "
         << tableID
         // << " (" << SesameTab.tableName[ SesameTab.returnTypes(tableID) ] << ")"
         << endl;
    if (eosprop == std::string("Uic_DT")) {
      cout << "  Specific Ion Internal Energy = "
           << spEospac->getSpecificIonInternalEnergy(temp, dens) << " kJ/g\n"
           << "  Ion Heat Capacity            = " << spEospac->getIonHeatCapacity(temp, dens)
           << "kJ/g/keV" << endl;
    } else if (eosprop == std::string("Ktc_DT")) {
      cout << "  Electron thermal conductivity = "
           << spEospac->getElectronThermalConductivity(temp, dens) << " /s/cm." << endl;
    }
  }

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  // Process command line arguments:
  rtt_dsxx::XGetopt::csmap long_options;
  long_options['v'] = "version";
  long_options['h'] = "help";
  std::map<char, std::string> help_strings;
  help_strings['h'] = "print this message.";
  help_strings['v'] = "print version information and exit.";
  rtt_dsxx::XGetopt program_options(argc, argv, long_options, help_strings);

  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'v': // --version
      std::cout << argv[0] << ": version " << rtt_dsxx::release() << std::endl;
      return 0;

    case 'h': // --help
      std::cout << argv[0] << ": version " << rtt_dsxx::release() << "\nUsage: QueryEospac\n"
                << "Follow the prompts to print equation-of-state data to the screen." << std::endl;
      return 0;

    default:
      break; // nothing to do.
    }
  }

  try {
    // >>> Run the application
    query_eospac();
  } catch (rtt_cdi_eospac::EospacException &err) {
    std::cout << "EospacException ERROR: While running " << argv[0] << ", " << err.what()
              << std::endl;
    return 1;
  } catch (rtt_dsxx::assertion &err) {
    std::cout << "ERROR: While running " << argv[0] << ", " << err.what() << std::endl;
    return 1;
  } catch (...) {
    std::cout << "ERROR: While running " << argv[0] << ", "
              << "An unknown exception was thrown on processor " << std::endl;
    return 1;
  }
  return 0;
}

//------------------------------------------------------------------------------------------------//
// end of QueryEospac.cc
//------------------------------------------------------------------------------------------------//
