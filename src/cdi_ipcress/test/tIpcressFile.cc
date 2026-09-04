//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ipcress/test/tIpcressFile.cc
 * \author Kelly Thompson
 * \brief  Ipcress file test
 * \note   Copyright (C) 2011-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "cdi_ipcress/IpcressFile.hh"
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// cSpell:ignore IPCRESS, matid, tgrid

namespace {

using namespace std;
using rtt_cdi_ipcress::IpcressFile;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

//! Tests the IpcressFile constructor and access routines.
void ipcress_file_test(rtt_dsxx::ScalarUnitTest &ut) {

  const string op_data_file = ut.getTestSourcePath() + "two-mats.ipcress";

  // Start the test.
  cout << "\nTesting the IpcressFile component of the cdi_ipcress package.\n";

  // Create a IpcressFile Object
  cout << "Creating a Ipcress File object\n\n";
  shared_ptr<IpcressFile> const spGF = std::make_shared<rtt_cdi_ipcress::IpcressFile>(op_data_file);

  // Test the new object to verify the constructor and accessors.

  vector<size_t> matIDs = spGF->getMatIDs();
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
    FAILMSG("Material 5500 shouldn't exist in the data file. Looks like we have a problem.");
  } else {
    PASSMSG("Function correctly identified material 5500 as being absent from IPCRESS file.");
  }

  if (spGF->getDataFilename() == op_data_file) {
    PASSMSG("Data filename set and retrieved correctly.");
  } else {
    FAILMSG("Data filename either not set correctly or not retrieved correctly.");
  }

  if (spGF->getNumMaterials() == 2) {
    PASSMSG("Found the correct number of materials in the data file.");
  } else {
    FAILMSG("Did not find the correct number of materials in the data file.");
  }

  cout << "\nMaterials found in the data file:\n";

  for (size_t i = 0; i < spGF->getNumMaterials(); ++i)
    cout << "  Material " << i << " has the identification number " << spGF->getMatIDs()[i] << '\n';

  // Retrieve a list of fields available for material 10001
  {
    size_t const matid(10001);
    vector<string> fieldNames = spGF->listDataFieldNames(matid);
    cout << "\nMaterial 0 (10001) provides the following fields:\n";
    for (auto &fieldName : fieldNames)
      cout << "   " << fieldName << "\n";
    cout << '\n';

    if (fieldNames[0] == string("tgrid"))
      PASSMSG("Found tgrid in the list of fields for mat 10001.");
    else
      FAILMSG("Did not find tgrid in the list of fields for mat 10001.");

    vector<double> tgrid = spGF->getData(matid, fieldNames[0]);
    cout << "\nMaterial 0 (10001)'s tgrid field has " << tgrid.size() << " entries: \n{ ";
    for (size_t i = 0; i < tgrid.size() - 1; ++i)
      cout << tgrid[i] << ", ";
    cout << tgrid[tgrid.size() - 1] << " }\n\n";

    vector<double> tgrid_expect = {0.01, 0.2575, 0.505, 0.7525, 1.0};
    if (rtt_dsxx::soft_equiv(tgrid_expect.begin(), tgrid_expect.end(), tgrid.begin(), tgrid.end()))
      PASSMSG("tgrid for mat 10001 has the expected values.");
    else
      FAILMSG("tgrid for mat 10001 does not have the expected values.");
  }
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    ipcress_file_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tIpcressFile.cc
//------------------------------------------------------------------------------------------------//
