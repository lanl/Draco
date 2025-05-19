//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/NDI_MultigroupNeutron.hh
 * \author Kyle S. Beling
 * \date   2025 March 18
 * \brief  NDI_MultigroupNeutron class definition.
 * \note   Copyright (C) 2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef cdi_ndi_NDI_MultigroupNeutron_hh
#define cdi_ndi_NDI_MultigroupNeutron_hh

#include "NDI_Base.hh"
#include "dsxx/Assert.hh"
#include "dsxx/DracoStrings.hh"
#include "dsxx/Soft_Equivalence.hh"
#include <algorithm>
#include <memory>
#include <numeric>
#include <string>

namespace rtt_cdi_ndi {

//================================================================================================//
/*!
 * \class NDI_MultigroupNeutron
 *
 * \brief Class for wrapping NDI calls to NDI MultigroupNeutron data. Reads data, constructs internal storage
 *        amenable to radiation calculations, and provides accessors.
 * \example cdi_ndi/test/tstNDI_MultigroupNeutron.cc
 */
//================================================================================================//

class NDI_MultigroupNeutron : public NDI_Base {

public:
  //! Constructor
  NDI_MultigroupNeutron(const std::string &library_in, int target_zaid_in,
                        const std::vector<double> &mg_e_bounds_in);
  //! Constructor
  NDI_MultigroupNeutron(const std::string &gendir_in, const std::string &library_in,
                        int target_zaid_in, const std::vector<double> &mg_e_bounds_in);

  //! Default destructor
  ~NDI_MultigroupNeutron() = default;

  //! Disable default constructor
  NDI_MultigroupNeutron() = delete;

  //! Disable copy constructor
  NDI_MultigroupNeutron(const NDI_MultigroupNeutron &rhs) = delete;

  //! Disable move constructor
  NDI_MultigroupNeutron(NDI_MultigroupNeutron &&rhs) noexcept = delete;

  //! Disable assignment
  NDI_MultigroupNeutron &operator=(const NDI_MultigroupNeutron &rhs) = delete;

  //! Disable move assignment
  NDI_MultigroupNeutron &operator=(NDI_MultigroupNeutron &&rhs) noexcept = delete;

  // >>>>>>>>>> TOTAL REACTION DATA ACCESSORS >>>>>>>>>>

  //! Multigroup total cross section
  [[nodiscard]] std::vector<double> totalCrossSection() const { return this->sigma_tot; }

  //! Energy-dependent total cross section
  [[nodiscard]] double totalCrossSection(double incident_energy) const;

  //! Multigroup total kermas
  [[nodiscard]] std::vector<double> kermas() const { return this->group_kermas; }

  //! Transfer matrix
  [[nodiscard]] std::vector<std::vector<std::vector<double>>> scatteringMatrix() const {
    return this->total_scattering_matrix;
  }

  //! Total PDF
  [[nodiscard]] std::vector<double> pdf(double incident_energy) const;

  // >>>>>>>>>> REACTION SPECIFIC DATA ACCESSORS >>>>>>>>>>

  //! List of products for a specific reaction
  [[nodiscard]] std::vector<int> listOfProducts(const int mt_number) const {
    return this->reaction_products[find_reaction_index(mt_number)];
  }

  //! List of product multiplicites for a specific reaction
  [[nodiscard]] std::vector<int> productMultiplicities(const int mt_number) const {
    return this->reaction_multiplicities[find_reaction_index(mt_number)];
  }

  //! Q-value for a specific reaction
  [[nodiscard]] double qValue(const int mt_number) const {
    return this->reaction_q_val[find_reaction_index(mt_number)];
  }

  //! Total cross section for a specific reaction
  [[nodiscard]] std::vector<double> totalCrossSection(const int mt_number) const {
    return this->reaction_total_sigmas[find_reaction_index(mt_number)];
  }

  //! Total cross section for a specific reaction and energy
  [[nodiscard]] double totalCrossSection(const int mt_number, double incident_energy) const;

  // >>>>>>>>>> PRODUCT SPECIFIC DATA ACCESSORS >>>>>>>>>>

  //! Multigroup kermas for a specific product
  [[nodiscard]] std::vector<double> kermas(const int product_zaid) const {
    return this->product_kermas[find_product_index(product_zaid)];
  }

  //! Transfer matrix for a specific product
  [[nodiscard]] std::vector<std::vector<std::vector<double>>>
  scatteringMatrix(const int product_zaid) const {
    return this->product_scattering_matrix[find_product_index(product_zaid)];
  }

  //! Return scattering PDF for a specific product and incident energy
  [[nodiscard]] std::vector<double> pdf(const int product_zaid, const double incident_energy) const;

  // >>>>>>>>>> OTHER USEFUL ACCESSORS >>>>>>>>>>

  //! Target zaid
  [[nodiscard]] int targetZAID() const { return this->target_zaid; }

  //! List of reactions
  [[nodiscard]] std::vector<int> listOfReactions() const { return this->reaction_mts; }

  //! List of all products
  [[nodiscard]] std::vector<int> listOfAllProducts() const { return this->product_list; }

  //! Group bounds
  [[nodiscard]] std::vector<double> groupBounds() const { return this->coarse_group_bounds; }

  //! Group velocities
  [[nodiscard]] std::vector<double> groupVelocities() const { return this->group_velocities; }

private:
  //! Target zaid
  int target_zaid = 0;

  //! Coarse group bounds (MeV)
  std::vector<double> coarse_group_bounds{};

  //! Group averaged sigmaT (barns), velocities, kermas
  std::vector<double> sigma_tot{}, group_velocities{}, group_kermas{};

  //! Total scattering matrix
  std::vector<std::vector<std::vector<double>>> total_scattering_matrix{};

  //! Reaction MT number
  std::vector<int> reaction_mts{};

  //! Reaction Q-value (MeV)
  std::vector<double> reaction_q_val{};

  //! Reaction products, and multiplicities
  std::vector<std::vector<int>> reaction_products{}, reaction_multiplicities{};

  //! Group averaged sigma reaction (barns)
  std::vector<std::vector<double>> reaction_total_sigmas{};

  //! Labels (ZAIDs) for reaction products
  std::vector<int> product_list{};

  //! Kermas and total cross sections (barns) for each product
  std::vector<std::vector<double>> product_kermas{};

  //! Scattering matrix for each product (barns)
  std::vector<std::vector<std::vector<std::vector<double>>>> product_scattering_matrix{};

  // >>> MEMBER FUNCTIONS >>>

  //! Find reaction index
  [[nodiscard]] int find_reaction_index(const int mt_number) const {
    Require(std::count(this->reaction_mts.begin(), this->reaction_mts.end(), mt_number) == 1);
    auto it = std::find(this->reaction_mts.begin(), this->reaction_mts.end(), mt_number);
    return static_cast<int>(std::distance(this->reaction_mts.begin(), it));
  }

  //! Find product index
  [[nodiscard]] int find_product_index(const int product_zaid) const {
    Require(std::count(this->product_list.begin(), this->product_list.end(), product_zaid) == 1);
    auto it = std::find(this->product_list.begin(), this->product_list.end(), product_zaid);
    return static_cast<int>(std::distance(this->product_list.begin(), it));
  }

  //! Find bin
  [[nodiscard]] int find_bin(const std::vector<double> &grid, const double value) const {
    Require(value < grid.front());
    Require(value > grid.back());
    const int bin =
        static_cast<int>(std::upper_bound(grid.rbegin(), grid.rend(), value) - grid.rbegin());
    return static_cast<int>(grid.size()) - 1 - bin;
  }

  //! NDI
  void ndi_throw_rt_error(int const ndi_err) {
    if (ndi_err != 0) {
      std::string errmsg("NDI error code: ");
      errmsg += rtt_dsxx::to_string(ndi_err);
      throw std::runtime_error(errmsg);
    }
  }

  //! Map course group onto fine group
  [[nodiscard]] std::vector<double> remapCoarseOntoFine(const std::vector<double> &coarse_bounds,
                                                        const std::vector<double> &fine_bounds) {
    std::vector<double> new_energy_grid(coarse_bounds.size());

    // Find closest fine group value and index
    for (size_t g = 0; g < coarse_bounds.size(); ++g) {
      new_energy_grid[g] =
          *std::min_element(fine_bounds.begin(), fine_bounds.end(), [&](double a, double b) {
            return std::abs(a - coarse_bounds[g]) < std::abs(b - coarse_bounds[g]);
          });
    }
    return new_energy_grid;
  }

  //! Remap vector into 2D matrix
  template <typename T>
  [[nodiscard]] std::vector<std::vector<T>> vectorToSquareMatrix(const std::vector<T> &vec,
                                                                 size_t size) const {
    std::vector<std::vector<T>> matrix(size, std::vector<T>(size, 0.0));
    for (size_t i = 0; i < size; ++i) {
      for (size_t j = 0; j < size; ++j) {
        matrix[i][j] = vec[i * size + j];
      }
    }
    return matrix;
  }

  // Only implemented if NDI is found
#ifdef NDI_FOUND
  void load_ndi();
#endif
};

} // namespace rtt_cdi_ndi

#endif // cdi_ndi_NDI_MultigroupNeutron_hh

//------------------------------------------------------------------------------------------------//
// End cdi_ndi/NDI_MultigroupNeutron.hh
//------------------------------------------------------------------------------------------------//
