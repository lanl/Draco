//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ipcress/test/tIpcressFile.cc
 * \author Kelly Thompson
 * \date   Fri Oct 12 15:39:39 2001
 * \brief  Ipcress file test
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "cdi_ipcress_test.hh"
#include "cdi_ipcress/IpcressFile.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"

using namespace std;

using rtt_cdi_ipcress::IpcressFile;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

//! Tests the IpcressFile constructor and access routines.
void ipcress_file_test(rtt_dsxx::ScalarUnitTest &ut) {

  const string op_data_file = ut.getTestSourcePath() + "two-mats.ipcress";

  // Start the test.

  cout << "\nTesting the IpcressFile component of the "
       << "cdi_ipcress package." << endl;

  // Create a IpcressFile Object

  cout << "Creating a Ipcress File object\n" << endl;

  shared_ptr<IpcressFile> spGF(new rtt_cdi_ipcress::IpcressFile(op_data_file));

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

  cout << "\nMaterials found in the data file:" << endl;

  for (size_t i = 0; i < spGF->getNumMaterials(); ++i)
    cout << "  Material " << i << " has the identification number " << spGF->getMatIDs()[i] << endl;

  // Retrieve a list of fields available for material 10001
  {
    size_t const matid(10001);
    vector<string> fieldNames = spGF->listDataFieldNames(matid);
    cout << "\nMaterial 0 (10001) provides the following fields:\n";
    for (auto &fieldName : fieldNames)
      cout << "   " << fieldName << "\n";
    cout << endl;

    if (fieldNames[0] == string("tgrid"))
      PASSMSG("Found tgrid in the list of fields for mat 10001.");
    else
      FAILMSG("Did not find tgrid in the list of fields for mat 10001.");

    vector<double> tgrid = spGF->getData(matid, fieldNames[0]);
    cout << "\nMaterial 0 (10001)'s tgrid field has " << tgrid.size() << " entries: \n{ ";
    for (size_t i = 0; i < tgrid.size() - 1; ++i)
      cout << tgrid[i] << ", ";
    cout << tgrid[tgrid.size() - 1] << " }\n" << endl;

    vector<double> tgrid_expect = {0.01, 0.2575, 0.505, 0.7525, 1.0};
    if (rtt_dsxx::soft_equiv(tgrid_expect.begin(), tgrid_expect.end(), tgrid.begin(), tgrid.end()))
      PASSMSG("tgrid for mat 10001 has the expected values.");
    else
      FAILMSG("tgrid for mat 10001 does not have the expected values.");
  }
}

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
