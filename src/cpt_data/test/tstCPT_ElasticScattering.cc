//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cpt_data/test/tstCPT_ElasticScattering.cc
 * \author Kyle S. Beling
 * \brief  Elastic scattering data test.
 * \note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cpt_data/CPT_Collision.hh"
#include "cpt_data/CPT_ElasticScattering.hh"

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"

#include <cstddef>
#include <string>
#include <vector>

using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void elastic_scattering_test(rtt_dsxx::UnitTest &ut) {
  std::string const dictionary_path = ut.getTestSourcePath();
  std::string const dataset = "Example";
  std::string const reaction_network = "Big 5";

  const std::vector<double> energy_grid{2.0e-2, 2.0e-1, 5.0e-1, 2.0e+0, 3.0e+1};

  // Build ElasticScattering object
  rtt_cpt_data::CPT_ElasticScattering const elastic(dictionary_path, dataset, reaction_network,
                                                    energy_grid);

  FAIL_IF_NOT(elastic.get_dataset() == "Example");
  FAIL_IF_NOT(elastic.reaction_network_name() == "Big 5");
  FAIL_IF_NOT(elastic.get_library() == "elastic_scattering");

  // Get list of collisions
  std::vector<std::string> list_of_collisions = elastic.list_of_collisions();
  FAIL_IF_NOT(list_of_collisions[0] == "d+t");

  // Get "d+t" collision
  const rtt_cpt_data::CPT_Collision collision = elastic.get_collision("d+t");
  FAIL_IF_NOT(collision.name() == "d+t");

  // Check scattered and recoiling particle
  FAIL_IF_NOT(collision.scattered_particle() == 1002);
  FAIL_IF_NOT(collision.recoiling_particle() == 1003);

  // Check sigma_total
  FAIL_IF_NOT(soft_equiv(collision.total_cross_section(1.0e-1), 1.6000e+03, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(collision.total_cross_section(1.0e+0), 2.0000e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(collision.total_cross_section(1.0e+1), 5.5000e-01, 1.0e-4));

  // Check pdf for deuteron
  std::vector<double> reference_pdf_deuteron{1.4414e-01, 4.0541e-01, 4.5045e-01, 0.0000e+00};
  std::vector<double> result_pdf_deuteron = collision.pdf(1002, 1.0e+0);
  FAIL_IF_NOT(result_pdf_deuteron.size() == reference_pdf_deuteron.size());
  for (size_t i = 0; i < reference_pdf_deuteron.size(); ++i) {
    FAIL_IF_NOT(soft_equiv(result_pdf_deuteron[i], reference_pdf_deuteron[i], 1.0e-4));
  }

  // Check pdf for triton
  std::vector<double> reference_pdf_triton{5.0000e-01, 2.3333e-01, 2.6667e-01, 0.0000e+00};
  std::vector<double> result_pdf_triton = collision.pdf(1003, 1.0e+0);
  FAIL_IF_NOT(result_pdf_triton.size() == reference_pdf_triton.size());
  for (size_t i = 0; i < reference_pdf_triton.size(); ++i) {
    FAIL_IF_NOT(soft_equiv(result_pdf_triton[i], reference_pdf_triton[i], 1.0e-4));
  }

  if (ut.numFails == 0)
    PASSMSG("Elastic scattering test passes.");
  else
    FAILMSG("Elastic scattering test fails.");
  return;
}

void custom_elastic_scattering_test(rtt_dsxx::UnitTest &ut) {
  std::string const dictionary_path = ut.getTestSourcePath();
  std::string const dataset = "Example";
  std::vector<std::string> const custom_reaction_network{"d+t"};

  const std::vector<double> energy_grid{2.0e-2, 2.0e-1, 5.0e-1, 2.0e+0, 3.0e+1};

  // Build ElasticScattering object
  rtt_cpt_data::CPT_ElasticScattering const elastic(dictionary_path, dataset,
                                                    custom_reaction_network, energy_grid);

  FAIL_IF_NOT(elastic.get_dataset() == "Example");
  FAIL_IF_NOT(elastic.reaction_network_name() == "Custom");
  FAIL_IF_NOT(elastic.get_library() == "elastic_scattering");

  // Get list of collisions
  std::vector<std::string> list_of_collisions = elastic.list_of_collisions();
  FAIL_IF_NOT(list_of_collisions[0] == "d+t");

  // Get "d+t" collision
  const rtt_cpt_data::CPT_Collision collision = elastic.get_collision("d+t");
  FAIL_IF_NOT(collision.name() == "d+t");

  // Check scattered and recoiling particle
  FAIL_IF_NOT(collision.scattered_particle() == 1002);
  FAIL_IF_NOT(collision.recoiling_particle() == 1003);

  // Check sigma_total
  FAIL_IF_NOT(soft_equiv(collision.total_cross_section(1.0e-1), 1.6000e+03, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(collision.total_cross_section(1.0e+0), 2.0000e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(collision.total_cross_section(1.0e+1), 5.5000e-01, 1.0e-4));

  // Check pdf for deuteron
  std::vector<double> reference_pdf_deuteron{1.4414e-01, 4.0541e-01, 4.5045e-01, 0.0000e+00};
  std::vector<double> result_pdf_deuteron = collision.pdf(1002, 1.0e+0);
  FAIL_IF_NOT(result_pdf_deuteron.size() == reference_pdf_deuteron.size());
  for (size_t i = 0; i < reference_pdf_deuteron.size(); ++i) {
    FAIL_IF_NOT(soft_equiv(result_pdf_deuteron[i], reference_pdf_deuteron[i], 1.0e-4));
  }

  // Check pdf for triton
  std::vector<double> reference_pdf_triton{5.0000e-01, 2.3333e-01, 2.6667e-01, 0.0000e+00};
  std::vector<double> result_pdf_triton = collision.pdf(1003, 1.0e+0);
  FAIL_IF_NOT(result_pdf_triton.size() == reference_pdf_triton.size());
  for (size_t i = 0; i < reference_pdf_triton.size(); ++i) {
    FAIL_IF_NOT(soft_equiv(result_pdf_triton[i], reference_pdf_triton[i], 1.0e-4));
  }

  if (ut.numFails == 0)
    PASSMSG("Custom elastic scattering test passes.");
  else
    FAILMSG("Custom elastic scattering test fails.");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    elastic_scattering_test(ut);
    custom_elastic_scattering_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstCPT_ElasticScattering.cc
//------------------------------------------------------------------------------------------------//
