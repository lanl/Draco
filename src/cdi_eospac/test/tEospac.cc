//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_eospac/test/tEospac.cc
 * \author Kelly Thompson
 * \date   Mon Apr 2 14:20:14 2001
 * \brief  Implementation file for tEospac
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi_eospac/Eospac.hh"
#include "cdi_eospac/EospacException.hh"
#include "ds++/DracoStrings.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <iomanip>
#include <sstream>

namespace rtt_cdi_eospac_test {

using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

//!  Tests the Eospac constructor and access routines.
void cdi_eospac_test(rtt_dsxx::UnitTest &ut) {
  // Start the test.

  std::cout << "\nTest of C++ code calling EOSPAC routines\n" << std::endl;

  // ---------------------------- //
  // Create a SesameTables object //
  // ---------------------------- //

  // The user must create a SesameTables object that links material ID numbers to EOSPAC data types
  // (each SesameTables object only contains lookups for one material).  If the user needs heat
  // capacity values for Al then he/she must create a SesameTables object for Aluminum and then
  // assign an aluminum material ID (e.g. 3717) to the enelc EOSPAC data type.  See the tests below
  // for more details.

  // Set the material identifier
  // This one is for Aluminum (03717)
  // Category-1 data (0) + Mat# 371 (Al) + Version # 7

  // See http://xweb.lanl.gov/projects/data/ for material ID information.

  // This matID for Al has lookup tables for prtot, entot, tptot, tntot, pntot, eptot, prelc, enelc,
  // tpelc, tnelc pnelc, epelc, prcld, and encld (see SesameTables.hh for an explanantion of these
  // keywords).  I need the table that contains enion lookups so that I can query for Cve() values.

  // Sesame Number 3717 provides data tables: 101 102 201 301 303 304 305 306 401
  int const Al3717 = 3717;

  // I also want to lookup the mean ion charge (EOS_Zfc_DT) which is found in sesame table 601.  Add
  // sesame number 23714 which provides data tables: 101 102 201 601 602 603 604
  int const Al23714 = 23714;

  // Create a SesameTables object for Aluminum.
  rtt_cdi_eospac::SesameTables AlSt;

  // Print a list of EosTables
  AlSt.printEosTableList();

  // Assign matID Al3717 to enion lookups (used for Cvi) for AlSt.  We can also assign these tables
  // when the Eospac object is created (see example below).

  // Also assign matID Al23714 for temperature-based electron thermal conductivity (tconde).
  AlSt.Uic_DT(Al3717).Ktc_DT(Al23714);
  AlSt.T_DUe(Al3717);  // getElectronTemperature(rho,Ue)
  AlSt.T_DUic(Al3717); // getIonTemperature(rho,Ue)

  // Verify that the assignments were made correctly.

  // Cvi (returnType=EOS_Uic_DT (=ES4enion)) should point to matID 3717. The user should never need
  // to access this function.  However Eospac.cc does and we need to test this functionality.

  if (AlSt.matID(EOS_Uic_DT) != 3717)
    FAILMSG("AlSt.matID(EOS_Uic_DT) points to the wrong matID.");

  // The temperature-based electron thermal conductivity (returnType=27=EOS_Ktc_DT) should point to
  // matID 23714.  The user should never need to access this function.  However Eospac.cc does and
  // we need to test this functionality.

  if (AlSt.matID(EOS_Ktc_DT) != 23714)
    FAILMSG("AlSt.matID(27) points to the wrong matID.");

  // ----------------------- //
  // Create an Eospac object //
  // ----------------------- //

  // An Eospac object allows the user to access EoS information about a material that has been
  // constructed in a SesameTable object.  The constructor for Eospac takes one argument: a
  // SesameTables object.

  std::shared_ptr<rtt_cdi_eospac::Eospac const> spEospac;

  // Try to instantiate the new Eospac object.  Simultaneously, we are assigned material IDs to more
  // SesameTable values.

  spEospac = std::make_shared<rtt_cdi_eospac::Eospac>(AlSt.Ue_DT(Al3717).Zfc_DT(Al23714));

  if (spEospac) {
    PASSMSG("shared_ptr to new Eospac object created.");
  } else {
    FAILMSG("Unable to create shared_ptr to new Eospac object.");

    // if construction fails, there is no reason to continue testing...
    return;
  }

  // --------------------------- //
  // Test scalar access routines //
  // --------------------------- //

  double const K2keV = 1.0 / 1.1604412E+7; // keV/Kelvin

  // All of these tests request an EoS value given a single temperature and a single density.

  // Retrieve an Electron internal energy value;

  double density = 1.0;      // g/cm^3
  double temperature = 5800; // K
  temperature *= K2keV;      // convert temps to keV

  double refValue = 1.0507392783; // kJ/g

  double specificElectronInternalEnergy =
      spEospac->getSpecificElectronInternalEnergy(temperature, density);
  double const tol(1.0e-10);

  if (soft_equiv(specificElectronInternalEnergy, refValue, tol))
    PASSMSG("getSpecificElectronInternalEnergy() test passed.");
  else {
    std::cout.precision(12);
    std::cout << "refValue = " << refValue << "\ntabValue = " << specificElectronInternalEnergy
              << std::endl;
    FAILMSG("getSpecificElectronInternalEnergy() test failed.");
  }

  // Retrieve an electron heat capacity (= dE/dT)

  // old refValue = 3146.719924188898; // kJ/g/keV
  refValue = 4101.9991645; // kJ/g/keV

  double heatCapacity = spEospac->getElectronHeatCapacity(temperature, density);

  if (soft_equiv(heatCapacity, refValue, tol))
    PASSMSG("getElectronHeatCapacity() test passed.");
  else
    FAILMSG("getElectronHeatCapacity() test failed.");

  // Retrieve an Ion Internal Energy

  refValue = 5.23391652028; // kJ/g

  double specificIonInternalEnergy = spEospac->getSpecificIonInternalEnergy(temperature, density);

  if (soft_equiv(specificIonInternalEnergy, refValue, tol))
    PASSMSG("getSpecificIonInternalEnergy() test passed.");
  else
    FAILMSG("getSpecificIonInternalEnergy() test failed.");

  // Retrieve an ion based heat capacity

  refValue = 6748.7474603; // kJ/g/keV

  heatCapacity = spEospac->getIonHeatCapacity(temperature, density);

  if (soft_equiv(heatCapacity, refValue, tol))
    PASSMSG("getIonHeatCapacity() test passed.");
  else
    FAILMSG("getIonHeatCapacity() test failed.");

  // Retrieve the number of free electrons per ion

  refValue = 12.8992087458; // electrons per ion

  double nfree = spEospac->getNumFreeElectronsPerIon(temperature, density);

  if (soft_equiv(nfree, refValue, tol))
    PASSMSG("getNumFreeElectronsPerIon() test passed.");
  else
    FAILMSG("getNumFreeElectronsPerIon() test failed.");

  // Retrieve the electron based thermal conductivity

  refValue = 1.38901721467e+29; // 1/s/cm

  double chie = spEospac->getElectronThermalConductivity(temperature, density);

  if (soft_equiv(chie, refValue, tol))
    PASSMSG("getElectronThermalConductivity() test passed.");
  else
    FAILMSG("getElectronThermalConductivity() test failed.");

  // Test the getElectronTemperature function

  double SpecificElectronInternalEnergy(1.0);                                              // kJ/g
  double Tout = spEospac->getElectronTemperature(density, SpecificElectronInternalEnergy); // keV

  double const Tegold(0.000487303450297301); // keV
  if (soft_equiv(Tout, Tegold))
    PASSMSG("getElectronTemperature() test passed for scalar.");
  else {
    std::ostringstream msg;
    msg << "getElectronTemperature() test failed for scalar.\n\tTout  = " << std::setprecision(16)
        << Tout << " keV\n\tTegold = " << Tegold << " keV";
    FAILMSG(msg.str());
  }

  // Test the getElectronTemperature function

  double SpecificIonInternalEnergy(10.0);                                 // kJ/g
  Tout = spEospac->getIonTemperature(density, SpecificIonInternalEnergy); // keV

  double const Tigold(0.001205608722470064); // keV
  if (soft_equiv(Tout, Tigold))
    PASSMSG("getIonTemperature() test passed for scalar.");
  else {
    std::ostringstream msg;
    msg << "getIonTemperature() test failed for scalar.\n\tTout  = " << std::setprecision(16)
        << Tout << " keV\n\tTigold = " << Tigold << " keV";
    FAILMSG(msg.str());
  }

  // --------------------------- //
  // Test vector access routines //
  // --------------------------- //

  // Set up simple temp and density vectors.  vtemp(i) will always be associated with vdensities(i).
  // In this case both tuples have identical data so that the returned results will also be
  // identical.

  std::vector<double> vtemps(2);
  std::vector<double> vdensities(2);

  vtemps[0] = temperature;
  vtemps[1] = temperature;
  vdensities[0] = density;
  vdensities[1] = density;

  // Retrieve electron based heat capacities for each set of (density, temperature) values.

  std::vector<double> vCve(2);
  vCve = spEospac->getElectronHeatCapacity(vtemps, vdensities);

  // Since the i=0 and i=1 tuples of density and temperature are identical the two returned heat
  // capacities should also soft_equiv.

  if (soft_equiv(vCve[0], vCve[1], tol))
    PASSMSG("getElectronHeatCapacity() test passed for vector state values.");
  else {
    std::cout.precision(12);
    std::cout << "refValue = " << refValue << "\ntabValue = " << vCve[0] << std::endl;
    FAILMSG("getElectronHeatCapacity() test failed for vector state values.");
  }

  // Retrieve the electron based thermal conductivity. This result should also match the scalar
  // value calculated above.

  std::vector<double> vchie(2);
  vchie = spEospac->getElectronThermalConductivity(vtemps, vdensities);

  if (soft_equiv(vchie[0], vchie[1], tol))
    PASSMSG("getElectronThermalConductivity() test passed for vector.");
  else
    FAILMSG("getElectronThermalConductivity() test failed for vector.");

  // Retrieve the ion based heat capacity. This result should also match the scalar value calculated
  // above.

  std::vector<double> vCvi(2);
  vCvi = spEospac->getIonHeatCapacity(vtemps, vdensities);

  if (soft_equiv(vCvi[0], vCvi[1], tol))
    PASSMSG("getIonHeatCapacity() test passed for vector.");
  else
    FAILMSG("getIonHeatCapacity() test failed for vector.");

  // Retrieve the number of free electrons per ion. This result should also match the scalar value
  // calculated above.

  std::vector<double> vnfepi(2);
  vnfepi = spEospac->getNumFreeElectronsPerIon(vtemps, vdensities);

  if (soft_equiv(vnfepi[0], vnfepi[1], tol))
    PASSMSG("getNumFreeElectronsPerIon() test passed for vector.");
  else
    FAILMSG("getNumFreeElectronsPerIon() test failed for vector.");

  // Retrieve the specific electron internal energy. This result should also match the scalar value
  // calculated above.

  std::vector<double> vUe(2);
  vUe = spEospac->getSpecificElectronInternalEnergy(vtemps, vdensities);

  if (soft_equiv(vUe[0], vUe[1], tol))
    PASSMSG("getSpecificElectronInternalEnergy() test passed for vector.");
  else
    FAILMSG("getSpecificElectronInternalEnergy() test failed for vector.");

  // Retrieve the specific ion internal energy. This result should also match the scalar value
  // calculated above.

  std::vector<double> vUi(2);
  vUi = spEospac->getSpecificIonInternalEnergy(vtemps, vdensities);

  if (soft_equiv(vUe[0], vUe[1], tol))
    PASSMSG("getSpecificIonInternalEnergy() test passed for vector.");
  else
    FAILMSG("getSpecificIonInternalEnergy() test failed for vector.");

  return;
} // end of runTest()

//------------------------------------------------------------------------------------------------//
// Test the cdi_eospac exception class
//------------------------------------------------------------------------------------------------//
void cdi_eospac_except_test(rtt_dsxx::UnitTest &ut) {
  // Material: Iron
  // 2140 provides tables: 101 102 103 201 301 303 304 305 306 401
  int const Fe2140(2140);
  // 12140 provides tables: 101 201 501 502 503 504 505
  int const Fe12140(12140);

  // Create a SesameTables object for Iron and setup the tables
  rtt_cdi_eospac::SesameTables FeSt;

  // Total Pressure (GPa) from table 301
  FeSt.Pt_DT(Fe2140);

  // Vapor Density on coexistence line (Mg/m^3) from table 401
  FeSt.Dv_T(Fe2140);

  // Calculated versus Interpolated Opacity Grid Boundary (from table 501)
  FeSt.Ogb(Fe12140);

  // Table Comments:
  // EOS_Comment
  // EOS_Info

  // Generate an Eospac object

  std::shared_ptr<rtt_cdi_eospac::Eospac const> spEospac(new rtt_cdi_eospac::Eospac(FeSt));

  // Print table information for Pt_DT:
  {
    std::ostringstream msg;
    spEospac->printTableInformation(EOS_Pt_DT, msg);
    std::cout << "\nTable information for Fe2140:\n" << msg.str() << std::endl;

    // Examine the output
    std::map<std::string, unsigned> wordcount = rtt_dsxx::get_word_count(msg, false);

    if (wordcount[std::string("EOS_Pt_DT")] == 1 && wordcount[std::string("2140")] == 1)
      PASSMSG("Information table for Pt_DT for Fe2140 printed correctly.");
    else
      FAILMSG("Information table for Pt_DT for Fe2140 failed to print.");
  }

  // Try to print data for a table that is not loaded. Print table information for Pt_DT:
  {

    std::cout << "\nAttempting to access a table that should not be available ..." << std::endl;
    bool exceptionThrown(false);
    try {
      std::ostringstream msg;
      spEospac->printTableInformation(EOS_Ktc_DT, msg);
      std::cout << "Table information for Fe2140:\n" << msg.str() << std::endl;
    } catch (rtt_cdi_eospac::EospacException &err) {
      exceptionThrown = true;
      std::ostringstream msg;
      msg << "Correctly caught an exception.  The message is " << err.what() << std::endl;
    }

    if (exceptionThrown)
      PASSMSG("Attempted access to unloaded table throws exception.");
    else
      FAILMSG("Failed to throw exception for invalid table access.");
  }

  return;
}

//------------------------------------------------------------------------------------------------//
//!  Tests the Eospac constructor and access routines.
void cdi_eospac_tpack(rtt_dsxx::UnitTest &ut) {
  // Start the test.
  std::cout << "\nTest the pack() function for cdi_eospac().\n" << std::endl;

  int const Al3717 = 3717;
  int const Al23714 = 23714;
  rtt_cdi_eospac::SesameTables AlSt;
  AlSt.Uic_DT(Al3717).Ktc_DT(Al23714).Ue_DT(Al3717).Zfc_DT(Al23714);
  std::shared_ptr<rtt_cdi_eospac::Eospac const> spEospac(new rtt_cdi_eospac::Eospac(AlSt));

  {
    // Test the SesameTables packer
    std::cout << "Packing a SesameTables object." << std::endl;
    std::vector<char> packed(AlSt.pack());

    // Create a new SesameTables by unpacking the packed data
    std::cout << "Unpacking a SesameTables object." << std::endl;
    rtt_cdi_eospac::SesameTables unpacked_AlSt(packed);

    if (AlSt.matID(EOS_Uic_DT) == unpacked_AlSt.matID(EOS_Uic_DT))
      PASSMSG("unpacked AlSt looks okay.");
    else
      FAILMSG("unpacked AlSt is wrong!");
  }

  {
    // Test the Eospac packer.
    std::cout << "Packing an Eospac object." << std::endl;
    std::vector<char> packed(spEospac->pack());

    // Create a new Eospac by unpacking the packed data.
    std::cout << "Unpacking an Eospac object." << std::endl;
    std::shared_ptr<rtt_cdi_eospac::Eospac const> spUnpacked_Eospac(
        new rtt_cdi_eospac::Eospac(packed));

    // Sanity Check
    double density = 5.0;     // g/cm^3
    double temperature = 0.1; // keV

    double cvi1(spEospac->getSpecificIonInternalEnergy(temperature, density));
    double cvi2(spUnpacked_Eospac->getSpecificIonInternalEnergy(temperature, density));
    if (soft_equiv(cvi1, cvi2, 1.0e-10))
      PASSMSG("unpacked spEospac looks okay.");
    else
      FAILMSG("unpacked spEospac is wrong.");
  }

  return;
}

} // end of namespace rtt_cdi_eospac_test

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    rtt_cdi_eospac_test::cdi_eospac_test(ut);
    rtt_cdi_eospac_test::cdi_eospac_except_test(ut);
    rtt_cdi_eospac_test::cdi_eospac_tpack(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tEospac.cc
//------------------------------------------------------------------------------------------------//
