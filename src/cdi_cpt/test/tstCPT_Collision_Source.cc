//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_cpt/test/tstCPT_Collision_Source.cc
 * \author Ben R. Ryan
 * \date   2024 Oct 9
 * \brief  CPT_Collision_Source test
 * \note   Copyright (C) 2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi_cpt/Collision_Source.hh"
#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "cdi/CPCommon.hh"
#include <cstdint>
#include <string>
#include <vector>

using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

// Deterministic Number Generator
class DNG {
  int i{0}; // Counter of # of calls to this DNG.
public:
  // Note: This is obviously only good for ~999 calls; plenty for a test.
  double ran() {
    i++;
    return static_cast<double>(i) / 1000.0;
  }
};

void simple_collision_test(rtt_dsxx::UnitTest &ut) {

  // Hydrogen atom (target)
  int32_t const hydrogen_zaid = 1001;
  // Alpha particle (projectile)
  int32_t const alpha_zaid = 2004;
  // Projectile energy
  double const projectile_energy = 1.0; // keV

  // Path to stable ACE data provided by Avalanche for testing
  std::string const ndatk_path_base = std::string(AVALANCHE_DIR) + std::string("/../../../") +
                                      std::string(avalanche_test_resources_DIR);
  std::vector<std::string> const ndatk_path_list = {ndatk_path_base,
                                                    ndatk_path_base + std::string("/data")};

  // Create CollisionSource object for this
  rtt_cdi_cpt::CollisionSource const cs(hydrogen_zaid, alpha_zaid, "endf71x", ndatk_path_list);

  // Create false RNG
  DNG dng;

  FAIL_IF_NOT(cs.get_target() == hydrogen_zaid);

  FAIL_IF_NOT(cs.get_projectile() == alpha_zaid);

  FAIL_IF_NOT(soft_equiv(cs.get_awr(), 9.991670e-01));

  auto temps = cs.get_temperatures();
  FAIL_IF_NOT(temps.size() == 7);
  std::vector<double> temps_ref = {8.617400e-09, 2.154300e-05, 2.530100e-05, 5.170400e-05,
                                   7.755600e-05, 1.034100e-04, 2.154300e-04};

  for (int n = 0; n < 7; n++) {
    FAIL_IF_NOT(soft_equiv(temps[n], temps_ref[n]));
  }

  FAIL_IF_NOT(soft_equiv(cs.get_cross_section(projectile_energy, 0.1) * 1.e24, 2.030482e1, 1.e-5));
  FAIL_IF_NOT(soft_equiv(cs.get_cross_section(projectile_energy, 1.0) * 1.e24, 2.030482e1, 1.e-5));
  FAIL_IF_NOT(
      soft_equiv(cs.get_cross_section(projectile_energy, 0.00001) * 1.e24, 2.030269e1, 1.e-5));

  // * Don't warn about missing include directives for avalanche, ndatk, and shacl.
  // * NOLINTBEGIN(misc-include-cleaner)

  // Retest with persistent interpolation caches at different temperatures
  auto interpCache = shacl::physics::makeInterpolationCache(cs.get_incident_energies(0.1),
                                                            projectile_energy / 1000.);
  FAIL_IF_NOT(soft_equiv(cs.get_cross_section(projectile_energy, 0.1, interpCache) * 1.e24,
                         2.030482e1, 1.e-5));
  interpCache = shacl::physics::makeInterpolationCache(cs.get_incident_energies(1.0),
                                                       projectile_energy / 1000.);
  FAIL_IF_NOT(soft_equiv(cs.get_cross_section(projectile_energy, 1.0, interpCache) * 1.e24,
                         2.030482e1, 1.e-5));
  interpCache = shacl::physics::makeInterpolationCache(cs.get_incident_energies(0.00001),
                                                       projectile_energy / 1000.);
  FAIL_IF_NOT(soft_equiv(cs.get_cross_section(projectile_energy, 0.00001, interpCache) * 1.e24,
                         2.030269e1, 1.e-5));

  // NOLINTEND(misc-include-cleaner)
  FAIL_IF_NOT(
      soft_equiv(cs.sample_target_energy_delta(projectile_energy, 0.1, dng), 9.98004e-01, 1e-5));
  FAIL_IF_NOT(
      soft_equiv(cs.sample_target_energy_delta(projectile_energy, 1.0, dng), 9.96008e-01, 1.e-5));
  FAIL_IF_NOT(soft_equiv(cs.sample_target_energy_delta(projectile_energy, 0.00001, dng),
                         9.94012e-01, 1.e-5));

  if (ut.numFails == 0) {
    PASSMSG("CPT_Collision_Source test passes.");
  } else {
    FAILMSG("CPT_Collision_Source test fails.");
  }
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    simple_collision_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstCPT_Collision_Source.cc
//------------------------------------------------------------------------------------------------//
