//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cpt_data/CPT_MultigroupCollapse.hh
 * \author Kyle S. Beling
 * \brief  CPT_MultigroupCollapse class definition.
 * \note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_CPT_multigroup_collapse
#define rtt_cdi_CPT_multigroup_collapse

#include "cpt_data/config.h" // Definition of JSON_FOUND

#include <vector>

namespace rtt_cpt_data {
//================================================================================================//
/*!
 * \class CPT_MultigroupCollapse
 * \brief Base class for collapse multigroup 
 *
 * Collapses multigroup data into the coarse group structure by finding the closest fine energy bin
 * boundaries. The fine group data should be fine enough that coarse groups do not overlap within a
 * single fine energy group.
 */
//================================================================================================//

#ifdef JSON_FOUND

class CPT_MultigroupCollapse {
  //! Coarse energy bounds
  std::vector<double> coarseBounds_{};

  //! Fine energy bounds
  std::vector<double> fineBounds_{};

  //! Mapping coarse onto fine indecies
  std::vector<int> mapCoarseToFine_{};

protected:
  //! Number of energy groups
  int numberEnergyGroups_{};

  //! Norm between requested and fit coarse energy bounds
  double normBetweenOldAndNew_{};

  //! Mapped energy bounds
  std::vector<double> mappedEnergyBounds_{};

public:
  CPT_MultigroupCollapse(const std::vector<double> &coarseBounds,
                         const std::vector<double> &fineBounds);

private:
  //! Check energy bounds
  void check_energy_bounds(const std::vector<double> &coarseGrid,
                           const std::vector<double> &fineGrid);

  //! Map coarse groups onto fine groups
  void remap_coarse_onto_fine(const std::vector<double> &coarseGrid,
                              const std::vector<double> &fineGrid);

  //! Compute norm between requested and fit coarse energy bounds
  void compute_norm_between_requested_and_fit_bounds();

protected:
  //! Elimante last element in tensor (only isotropic scattering for now)
  [[nodiscard]] std::vector<std::vector<double>> tensor_to_matrix(const nlohmann::json &tensor,
                                                                  int index) const;

  //! Collapse group average values
  [[nodiscard]] std::vector<double>
  collapse_group_average_values(const std::vector<double> &unCollapsed,
                                const std::vector<double> &weights) const;

  //! Collapse transfer matrices
  [[nodiscard]] std::vector<std::vector<double>>
  collapse_transfer_matrices(const std::vector<std::vector<double>> &unCollapsed,
                             const std::vector<double> &weights) const;

  //! Returns bin of a value
  [[nodiscard]] int get_index(const double &value) const;
};

#else

class CPT_MultigroupCollapse {};

#endif

} // namespace rtt_cpt_data

#endif // rtt_cdi_CPT_multigroup_collapse

//------------------------------------------------------------------------------------------------//
// End cpt_data/CPT_MultigroupCollapse.hh
//------------------------------------------------------------------------------------------------//
