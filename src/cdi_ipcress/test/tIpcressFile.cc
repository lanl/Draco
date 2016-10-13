//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_ipcress/test/tIpcressFile.cc
 * \author Kelly Thompson
 * \date   Fri Oct 12 15:39:39 2001
 * \brief  Ipcress file test
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../IpcressFile.hh"
#include "cdi_ipcress_test.hh"
#include "ds++/Release.hh"
#include "ds++/SP.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"

#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

using rtt_cdi_ipcress::IpcressFile;
using rtt_dsxx::SP;

#define PASSMSG(m) ut.passes(m)
#define FAILMSG(m) ut.failure(m)

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
/*!
 * \brief Tests the IpcressFile constructor and access routines.
 */
void ipcress_file_test(rtt_dsxx::ScalarUnitTest &ut) {
  // Ipcress data filename (IPCRESS format required)
  const std::string op_data_file = "Al_BeCu.ipcress";

  // Start the test.

  std::cout << std::endl
            << "Testing the IpcressFile component of the "
            << "cdi_ipcress package." << std::endl;

  // Create a IpcressFile Object

  std::cout << "Creating a Ipcress File object\n" << std::endl;

  SP<IpcressFile> spGF(new rtt_cdi_ipcress::IpcressFile(op_data_file));

  // Test the new object to verify the constructor and accessors.

  std::vector<size_t> matIDs = spGF->getMatIDs();
  if (matIDs[0] == 10001 && matIDs[1] == 10002)
    PASSMSG("Found two materials in IPCRESS file with expected IDs.");
  else
    FAILMSG("Did not find materials with expected IDs in IPCRESS file.");

  if (spGF->materialFound(10001))
    PASSMSG("Looks like material 10001 is in the data file.");
  else
    FAILMSG("Can't find material 10001 in the data file.");

  if (spGF->materialFound(5500)) // should fail
  {
    ostringstream message;
    message << "Material 5500 shouldn't exist in the data file."
            << "\n\tLooks like we have a problem.";
    FAILMSG(message.str());
  } else {
    ostringstream message;
    message << "Access function correctly identified material 5500"
            << "\n\tas being absent from IPCRESS file.";
    PASSMSG(message.str());
  }

  if (spGF->getDataFilename() == op_data_file) {
    PASSMSG("Data filename set and retrieved correctly.");
  } else {
    ostringstream message;
    message << "Data filename either not set correctly or not "
            << "retrieved correctly.";
    FAILMSG(message.str());
  }

  if (spGF->getNumMaterials() == 2) {
    PASSMSG("Found the correct number of materials in the data file.");
  } else {
    ostringstream message;
    message << "Did not find the correct number of materials in "
            << "the data file.";
    FAILMSG(message.str());
  }

  std::cout << std::endl << "Materials found in the data file:" << std::endl;

  for (size_t i = 0; i < spGF->getNumMaterials(); ++i)
    std::cout << "  Material " << i << " has the identification number "
              << spGF->getMatIDs()[i] << std::endl;

  // Retrieve a list of fields available for material 10001
  {
    size_t const matid(10001);
    std::vector<std::string> fieldNames = spGF->listDataFieldNames(matid);
    std::cout << "\nMaterial 0 (10001) provides the following fields:\n";
    for (size_t i = 0; i < fieldNames.size(); ++i)
      std::cout << "   " << fieldNames[i] << "\n";
    std::cout << std::endl;

    if (fieldNames[0] == std::string("tgrid"))
      PASSMSG("Found tgrid in the list of fields for mat 10001.");
    else
      FAILMSG("Did not find tgrid in the list of fields for mat 10001.");

    std::vector<double> tgrid = spGF->getData(matid, fieldNames[0]);
    std::cout << "\nMaterial 0 (10001)'s tgrid field has " << tgrid.size()
              << " entries: \n{ ";
    for (size_t i = 0; i < tgrid.size() - 1; ++i)
      std::cout << tgrid[i] << ", ";
    std::cout << tgrid[tgrid.size() - 1] << " }\n" << std::endl;

    std::vector<double> tgrid_expect(10);
    tgrid_expect[0] = 0.0005;
    tgrid_expect[1] = 0.0015;
    tgrid_expect[2] = 0.004;
    tgrid_expect[3] = 0.0125;
    tgrid_expect[4] = 0.04;
    tgrid_expect[5] = 0.125;
    tgrid_expect[6] = 0.4;
    tgrid_expect[7] = 1.25;
    tgrid_expect[8] = 4.0;
    tgrid_expect[9] = 15.0;
    if (rtt_dsxx::soft_equiv(tgrid_expect.begin(), tgrid_expect.end(),
                             tgrid.begin(), tgrid.end()))
      PASSMSG("tgrid for mat 10001 has the expected values.");
    else
      FAILMSG("tgrid for mat 10001 does not have the expected values.");
  }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    ipcress_file_test(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tIpcressFile.cc
//---------------------------------------------------------------------------//
