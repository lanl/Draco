//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/test/tstTabular_CP_Eloss.cc
 * \author Ben R. Ryan
 * \date   2019 Nov 18
 * \brief  Tabular_CP_Eloss test
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi/CDI.hh"
#include "cdi_CPEloss/Tabular_CP_Eloss.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/dbc.hh"
#include <sstream>

using rtt_cdi_cpeloss::Tabular_CP_Eloss;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void dedx_table_test(rtt_dsxx::UnitTest &ut) {

  // Datatable filename
  std::string filename_in = ut.getTestSourcePath() + "001-H-001";
  // Deuterium target
  int deuterium_zaid = 1002;
  double deuterium_mass = 3.34447643e-24;
  rtt_cdi::CParticle target_in(deuterium_zaid, deuterium_mass);
  // Proton projectile
  int proton_zaid = 1001;
  double proton_mass = 1.6726219e-24;
  rtt_cdi::CParticle projectile_in(proton_zaid, proton_mass);

  Tabular_CP_Eloss eloss_mod(filename_in, target_in, projectile_in,
                             rtt_cdi::CPModelAngleCutoff::NONE);

  // Model type better be tabular:
  FAIL_IF_NOT(eloss_mod.getModelType() == rtt_cdi::CPModelType::TABULAR_ETYPE);

  // Tabular data
  FAIL_IF_NOT(eloss_mod.is_data_in_tabular_form());

  // Check that grid accessors agree
  FAIL_IF_NOT(eloss_mod.getTemperatureGrid().size() == eloss_mod.getNumTemperatures());
  FAIL_IF_NOT(eloss_mod.getDensityGrid().size() == eloss_mod.getNumDensities());
  FAIL_IF_NOT(eloss_mod.getEnergyGrid().size() == eloss_mod.getNumEnergies());

  // Check that class grid dimensions match file
  FAIL_IF_NOT(eloss_mod.getNumEnergies() == 4);
  FAIL_IF_NOT(eloss_mod.getNumDensities() == 5);
  FAIL_IF_NOT(eloss_mod.getNumTemperatures() == 6);

  // Data file should not be empty
  FAIL_IF(eloss_mod.getDataFilename().empty());

  // Get eloss value for almost first (1,1,1) grid point
  {
    double energy = 1.384272;
    double density = 3.344490e-01;
    double temperature = 3.981051e-04;
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_mod.getEloss(temperature, density, energy),
                                     0.1958064213742134, 1.e-8));
  }

  // Get eloss value for almost last (2,3,4) grid point
  {
    double energy = 1.384273e+01;
    double density = 3.344495e+03;
    double temperature = 2.511868e+01;
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_mod.getEloss(temperature, density, energy),
                                     8.502486928162006370e+04, 1.e-8));
  }

  // Get eloss value for a point between grid points (1.5,2.5,3.5, i.e. requiring linear interpolation)
  {
    double energy = 4.377453e+00;
    double density = 3.344494e+02;
    double temperature = 3.981044e+00;
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_mod.getEloss(temperature, density, energy),
                                     139104.1982932578, 1.e-8));
  }

  if (ut.numFails == 0) {
    PASSMSG("Tabular_CP_Eloss test passes.");
  } else {
    FAILMSG("Tabular_CP_Eloss test fails.");
  }
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    dedx_table_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tst_Tabular_CP_Eloss.cc
//------------------------------------------------------------------------------------------------//
