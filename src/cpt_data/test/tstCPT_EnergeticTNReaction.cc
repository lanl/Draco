//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cpt_data/test/tstCPT_EnergeticTNReaction.cc
 * \author Kyle S. Beling
 * \date   Jan 15 2026
 * \brief  energetic TN reaction test.
 * \note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cpt_data/CPT_Collision.hh"
#include "cpt_data/CPT_EnergeticTNReaction.hh"

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

void rif_test(rtt_dsxx::UnitTest &ut) {
  std::string const dictionary_path = ut.getTestSourcePath();
  std::string const dataset = "Example";
  std::string const reaction_network = "Big 5";

  const std::vector<double> energy_grid{2.0e-2, 2.0e-1, 5.0e-1, 2.0e+0, 3.0e+1};

  // Build ReactionInFlight object
  rtt_cpt_data::CPT_EnergeticTNReaction const rif(dictionary_path, dataset, reaction_network,
                                                  energy_grid);

  FAIL_IF_NOT(rif.get_dataset() == "Example");
  FAIL_IF_NOT(rif.reaction_network_name() == "Big 5");
  FAIL_IF_NOT(rif.get_library() == "reaction_in_flight");

  // Get list reactions
  std::vector<std::string> list_of_reactions = rif.list_of_reactions();
  FAIL_IF_NOT(list_of_reactions[0] == "d+t->n+a");

  // Get "d+t->n+a" reaction
  const rtt_cpt_data::CPT_Collision reaction = rif.get_reaction("d+t->n+a");

  // Check reaction information
  FAIL_IF_NOT(reaction.name() == "d+t->n+a");

  // Check scattered and recoiling particle
  FAIL_IF_NOT(reaction.scattered_particle() == 1);
  FAIL_IF_NOT(reaction.recoiling_particle() == 2004);

  // Check sigma_total
  FAIL_IF_NOT(soft_equiv(reaction.total_cross_section(1.0e-1), 2.1500e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(reaction.total_cross_section(1.0e+0), 2.0000e-01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(reaction.total_cross_section(1.0e+1), 6.0000e-02, 1.0e-4));

  // Check pdf for alpha-particle
  std::vector<double> reference_pdf_alpha{0.0000e+00, 0.0000e+00, 1.6667e-01, 8.3333e-01};
  std::vector<double> result_pdf_alpha = reaction.pdf(2004, 1.0e+0);
  FAIL_IF_NOT(result_pdf_alpha.size() == reference_pdf_alpha.size());
  for (size_t i = 0; i < reference_pdf_alpha.size(); ++i) {
    FAIL_IF_NOT(soft_equiv(result_pdf_alpha[i], reference_pdf_alpha[i], 1.0e-4));
  }

  // Check pdf for neutron
  std::vector<double> reference_pdf_neutron{0.0000e+00, 0.0000e+00, 0.0000e+00, 1.0000e+00};
  std::vector<double> result_pdf_neutron = reaction.pdf(1, 1.0e+0);
  FAIL_IF_NOT(result_pdf_neutron.size() == reference_pdf_neutron.size());
  for (size_t i = 0; i < reference_pdf_neutron.size(); ++i) {
    FAIL_IF_NOT(soft_equiv(result_pdf_neutron[i], reference_pdf_neutron[i], 1.0e-4));
  }

  if (ut.numFails == 0)
    PASSMSG("energetic TN reaction test passes.");
  else
    FAILMSG("energetic TN reaction test fails.");
  return;
}

void custom_reaction_network_test(rtt_dsxx::UnitTest &ut) {
  std::string const dictionary_path = ut.getTestSourcePath();
  std::string const dataset = "Example";
  std::vector<std::string> const custom_reaction_network{"d+t->n+a"};

  const std::vector<double> energy_grid{2.0e-2, 2.0e-1, 5.0e-1, 2.0e+0, 3.0e+1};

  // Build ReactionInFlight object
  rtt_cpt_data::CPT_EnergeticTNReaction const rif(dictionary_path, dataset, custom_reaction_network,
                                                  energy_grid);

  FAIL_IF_NOT(rif.get_dataset() == "Example");
  FAIL_IF_NOT(rif.reaction_network_name() == "Custom");
  FAIL_IF_NOT(rif.get_library() == "reaction_in_flight");

  // Get list reactions
  std::vector<std::string> list_of_reactions = rif.list_of_reactions();
  FAIL_IF_NOT(list_of_reactions[0] == "d+t->n+a");

  // Get "d+t->n+a" reaction
  const rtt_cpt_data::CPT_Collision reaction = rif.get_reaction("d+t->n+a");

  // Check reaction information
  FAIL_IF_NOT(reaction.name() == "d+t->n+a");

  // Check scattered and recoiling particle
  FAIL_IF_NOT(reaction.scattered_particle() == 1);
  FAIL_IF_NOT(reaction.recoiling_particle() == 2004);

  // Check sigma_total
  FAIL_IF_NOT(soft_equiv(reaction.total_cross_section(1.0e-1), 2.1500e+00, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(reaction.total_cross_section(1.0e+0), 2.0000e-01, 1.0e-4));
  FAIL_IF_NOT(soft_equiv(reaction.total_cross_section(1.0e+1), 6.0000e-02, 1.0e-4));

  // Check pdf for alpha-particle
  std::vector<double> reference_pdf_alpha{0.0000e+00, 0.0000e+00, 1.6667e-01, 8.3333e-01};
  std::vector<double> result_pdf_alpha = reaction.pdf(2004, 1.0e+0);
  FAIL_IF_NOT(result_pdf_alpha.size() == reference_pdf_alpha.size());
  for (size_t i = 0; i < reference_pdf_alpha.size(); ++i) {
    FAIL_IF_NOT(soft_equiv(result_pdf_alpha[i], reference_pdf_alpha[i], 1.0e-4));
  }

  // Check pdf for neutron
  std::vector<double> reference_pdf_neutron{0.0000e+00, 0.0000e+00, 0.0000e+00, 1.0000e+00};
  std::vector<double> result_pdf_neutron = reaction.pdf(1, 1.0e+0);
  FAIL_IF_NOT(result_pdf_neutron.size() == reference_pdf_neutron.size());
  for (size_t i = 0; i < reference_pdf_neutron.size(); ++i) {
    FAIL_IF_NOT(soft_equiv(result_pdf_neutron[i], reference_pdf_neutron[i], 1.0e-4));
  }

  if (ut.numFails == 0)
    PASSMSG("Custom energetic TN reaction test passes.");
  else
    FAILMSG("Custom energetic TN reaction test fails.");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    rif_test(ut);
    custom_reaction_network_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstCPT_EnergeticTNReaction.cc
//------------------------------------------------------------------------------------------------//
