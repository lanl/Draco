//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cpt_data/CPT_MultigroupCollapse.cc
 * \author Kyle S. Beling
 * \brief  CPT_MultigroupCollapse member definitions.
 * \note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

// Definition of JSON_FOUND
#include "cpt_data/config.h" // NOLINT [misc-include-cleaner]

#ifdef JSON_FOUND

#include "cpt_data/CPT_MultigroupCollapse.hh"

#include "dsxx/Assert.hh"
#include "dsxx/dbc.hh"
#include <string>
#include <utility>

namespace rtt_cpt_data {

//================================================================================================//
// CONSTRUCTOR
//================================================================================================//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for JSON reader.
 *
 * \param[in] coarseBounds User defined energy group bounds that the fine bounds will be 
 *                         collapsed to in this class
 * \param[in] fineBounds Fine energy bounds that will be collapsed
 */
CPT_MultigroupCollapse::CPT_MultigroupCollapse(const std::vector<double> &coarseBounds,
                                               const std::vector<double> &fineBounds)
    : coarseBounds_(coarseBounds), fineBounds_(fineBounds) {
  // Ensure the coarse and fine bounds contain data
  Require(coarseBounds_.size() > 0);
  Require(fineBounds_.size() > 0);

  // Check that coarseBounds_ and fineBounds_ is monotonically increasing
  Insist(rtt_dsxx::is_strict_monotonic_increasing(coarseBounds_.begin(), coarseBounds_.end()),
         "Requested multigroup bounds are not strictly monotonic increasing!");
  Insist(rtt_dsxx::is_strict_monotonic_increasing(fineBounds_.begin(), fineBounds_.end()),
         "Multigroup bounds in the collision file are not strictly monotonic increasing!");

  // Ensure coarse bounds are within fine bounds
  Require(coarseBounds_.front() >= fineBounds_.front());
  Require(coarseBounds_.back() <= fineBounds_.back());

  // Check that coarse bounds are more coarse than the fine bounds
  Require(coarseBounds_.size() < fineBounds_.size());

  // Check that there is not more than one bound within any fine group
  check_energy_bounds(coarseBounds_, fineBounds_);

  // Remap coarse bounds onto fine bounds
  remap_coarse_onto_fine(coarseBounds, fineBounds);

  // Compute difference between requested and fit coarse energy bounds
  compute_norm_between_requested_and_fit_bounds();
}
//------------------------------------------------------------------------------------------------//

//================================================================================================//
// MEMBER FUNCTIONS
//================================================================================================//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Check coarse and fine energy grids to enure that only one coarse grid point falls within 
 *        any single fine energy group.
 *
 * \param[in] coarseGrid User defined energy group bounds that the fine bounds will be 
 *                       collapsed to in this class
 * \param[in] fineGrid Fine energy bounds of the multigroup data
 */
void CPT_MultigroupCollapse::check_energy_bounds(const std::vector<double> &coarseGrid,
                                                 const std::vector<double> &fineGrid) {
  // Loop through coarse bounds to see if more than one falls within a fineGrid bin
  for (int index = 0; index < static_cast<int>(coarseGrid.size() - 1); ++index) {
    // Current energy bounds
    double firstEnergyBound = coarseGrid[index];
    double secondEnergyBound = coarseGrid[index + 1];

    // Find locations of the closest bounds
    auto firstIterator = std::lower_bound(fineGrid.begin(), fineGrid.end(), firstEnergyBound);
    auto secondIterator = std::lower_bound(fineGrid.begin(), fineGrid.end(), secondEnergyBound);

    // If locations are the same we cannot collapse to the coarse group
    Insist(std::abs(*firstIterator - *secondIterator) > std::numeric_limits<double>::epsilon(),
           "Cannot collapse to coarse group structure: More than one coarse bound within a single "
           "fine group.");
  }
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Here the coarseGrid is remapped onto the fineGrid by selecting the closest fineGrid 
 *        point for each coarseGrid point. Note: the coarseGrid is not necessarily preserved,  
 *        but is instead preserved to the best of the fineGrid's points.
 *
 * \param[in] coarseGrid User defined energy group bounds that the fine bounds will be 
 *                       collapsed to in this class
 * \param[in] fineGrid Fine energy bounds of the multigroup data
 */
void CPT_MultigroupCollapse::remap_coarse_onto_fine(const std::vector<double> &coarseGrid,
                                                    const std::vector<double> &fineGrid) {
  // Update number of energy groups
  this->numberEnergyGroups_ = static_cast<int>(coarseGrid.size() - 1);

  // Remap coasre grid onto fine grid
  for (size_t g = 0; g < coarseGrid.size(); ++g) {
    auto value =
        std::min_element(fineGrid.begin(), fineGrid.end(), [&coarseGrid, g](double a, double b) {
          return std::abs(a - coarseGrid[g]) < std::abs(b - coarseGrid[g]);
        });
    this->mappedEnergyBounds_.emplace_back(*value);

    int index = static_cast<int>(std::distance(fineGrid.begin(), value));
    this->mapCoarseToFine_.emplace_back(index);
  }
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Determine how well the requested and fit coarse energy bounds match by using an L2 norm
 */
void CPT_MultigroupCollapse::compute_norm_between_requested_and_fit_bounds() {
  double result = 0.0;
  for (int g = 0; g < this->numberEnergyGroups_; ++g) {
    result += std::pow(this->coarseBounds_[g] - this->mappedEnergyBounds_[g], 2);
  }
  result = std::sqrt(result) / this->numberEnergyGroups_;
  this->normBetweenOldAndNew_ = result;
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Elimante last element in a nlohmann::json tensor (only isotropic scattering for now) 
 *        and recast as a std::vector<std::vector<double>> by evaluating the tensor as 
 *        tensor[:][:][index].
 *
 * \param[in] tensor A nlohmann::json tensor object
 * \param[in] index Integer of where to evalute the tensor
 *
 * \return A matrix consisting of tensor[:][:][index]
 */
[[nodiscard]] std::vector<std::vector<double>>
CPT_MultigroupCollapse::tensor_to_matrix(const nlohmann::json &tensor, int index) const {
  std::vector<std::vector<double>> matrix(fineBounds_.size() - 1,
                                          std::vector<double>(fineBounds_.size() - 1, 0.0));
  for (size_t i = 0; i < fineBounds_.size() - 1; ++i) {
    for (size_t j = 0; j < fineBounds_.size() - 1; ++j) {
      matrix[i][j] = tensor[i][j][index].get<double>();
    }
  }
  return matrix;
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Collapsed group averaged values on the fineBounds grid to the collapsed bounds 
 *        using the weights. The collapse algorithm is outlined in Reference:
 *        D. E. Cullen, “Nuclear data preparation,” in Handbook of Nuclear Engineering
 *
 * \param[in] unCollapsed Vector of un-collapsed values to be collapsed
 * \param[in] weights Fine bounds weights that are used to collapsed
 *
 * \return Vector of collapsed values
 */
[[nodiscard]] std::vector<double>
CPT_MultigroupCollapse::collapse_group_average_values(const std::vector<double> &unCollapsed,
                                                      const std::vector<double> &weights) const {
  std::vector<double> collapsed(this->numberEnergyGroups_, 0.0);
  for (int g = 0; g < this->numberEnergyGroups_; ++g) {
    int const dnBound = this->mapCoarseToFine_[g + 0];
    int const upBound = this->mapCoarseToFine_[g + 1];
    if (dnBound == upBound) {
      collapsed[g] = 0.0;
    } else {
      double weight_sum = 0.0;
      for (int i = dnBound; i < upBound; ++i) {
        double const EjMStar = std::max(this->mappedEnergyBounds_[g + 0], this->fineBounds_[i + 0]);
        double const EjPStar = std::min(this->mappedEnergyBounds_[g + 1], this->fineBounds_[i + 1]);
        double const pj = (EjPStar - EjMStar) / (this->fineBounds_[i + 1] - this->fineBounds_[i]);
        collapsed[g] += pj * unCollapsed[i] * weights[i];
        weight_sum += pj * weights[i];
      }

      Insist(weight_sum > std::numeric_limits<double>::epsilon(),
             "Cannot divide by zero: Weights summed to zero.");
      collapsed[g] /= weight_sum;
    }
  }
  return collapsed;
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Collapses a fine matrix of group averaged values to the coarse bounds 
 *        using the weights. The collapse algorithm is outlined in Reference:
 *        D. E. Cullen, “Nuclear data preparation,” in Handbook of Nuclear Engineering
 *
 * \param[in] unCollapsed Matrix of un-collapsed values to be collapsed
 * \param[in] weights Fine bounds weights that are used to collapsed
 *
 * \return Matrix of collapsed values
 */
[[nodiscard]] std::vector<std::vector<double>> CPT_MultigroupCollapse::collapse_transfer_matrices(
    const std::vector<std::vector<double>> &unCollapsed, const std::vector<double> &weights) const {
  std::vector<std::vector<double>> collapsed(this->numberEnergyGroups_,
                                             std::vector<double>(this->numberEnergyGroups_, 0.0));
  for (int g = 0; g < this->numberEnergyGroups_; ++g) {
    double const Eg_d = mappedEnergyBounds_[g];
    double const Eg_u = mappedEnergyBounds_[g + 1];
    int const dnBound_g = mapCoarseToFine_[g];
    int const upBound_g = mapCoarseToFine_[g + 1];
    for (int k = 0; k < this->numberEnergyGroups_; ++k) {
      double const Ek_d = mappedEnergyBounds_[k];
      double const Ek_u = mappedEnergyBounds_[k + 1];
      int const dnBound_k = mapCoarseToFine_[k];
      int const upBound_k = mapCoarseToFine_[k + 1];

      double weight_sum = 0.0;
      for (int i = dnBound_g; i < upBound_g; ++i) {
        double const Ei_m_star = std::max(Eg_d, fineBounds_[i]);
        double const Ei_p_star = std::min(Eg_u, fineBounds_[i + 1]);
        double const pi = (Ei_p_star - Ei_m_star) / (fineBounds_[i + 1] - fineBounds_[i]);

        // Sum weights only on last time through "j" iteration
        weight_sum = 0.0;
        for (int j = dnBound_k; j < upBound_k; ++j) {
          double const Ej_m_star = std::max(Ek_d, fineBounds_[j]);
          double const Ej_p_star = std::min(Ek_u, fineBounds_[j + 1]);
          double const pj = (Ej_p_star - Ej_m_star) / (fineBounds_[j + 1] - fineBounds_[j]);

          collapsed[g][k] += pi * pj * unCollapsed[i][j] * weights[j];
          weight_sum += pj * weights[j];
        }
      }
      Insist(weight_sum > std::numeric_limits<double>::epsilon(),
             "Cannot divide by zero: Weights summed to zero.");
      collapsed[g][k] /= weight_sum;
    }
  }
  return collapsed;
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Find energy bin where the energyIn value is located
 *
 * \param[in] energyIn Energy value [MeV]
 *
 * \return Energy bin where the energyIn value is located
 */
[[nodiscard]] int CPT_MultigroupCollapse::get_index(const double &energyIn) const {
  // Check that energyIn is in the energy grid
  Require(energyIn >= this->mappedEnergyBounds_.front());
  Require(energyIn <= this->mappedEnergyBounds_.back());

  // Locate the energy bin that energyIn is in
  auto it = std::lower_bound(this->mappedEnergyBounds_.begin(), this->mappedEnergyBounds_.end(),
                             energyIn);
  return static_cast<int>(std::distance(this->mappedEnergyBounds_.begin(), it) - 1);
}
//------------------------------------------------------------------------------------------------//

} // namespace rtt_cpt_data

#else

#endif

//------------------------------------------------------------------------------------------------//
// End cpt_data/Reaction.cc
//------------------------------------------------------------------------------------------------//
