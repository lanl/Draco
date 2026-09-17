//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/NDI_TNReaction.hh
 * \author Ben R. Ryan
 * \brief  NDI_TNReaction class definition.
 * \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef cdi_ndi_NDI_TNReaction_hh
#define cdi_ndi_NDI_TNReaction_hh

#include "NDI_Base.hh"
#include "dsxx/Assert.hh"
#include "dsxx/Soft_Equivalence.hh"
#include <algorithm>
#include <memory>
#include <string>
#include <utility>

namespace rtt_cdi_ndi {

//================================================================================================//
/*!
 * \class NDI_TNReaction
 *
 * \brief Class for wrapping NDI calls to NDI tn data. Reads data, constructs internal storage
 *        amenable to radiation calculations, and provides accessors.
 * \example cdi_ndi/test/tstNDI_TNReaction.cc
 */
//================================================================================================//

class NDI_TNReaction : public NDI_Base {

public:
  //! Constructor (default gendir path)
  NDI_TNReaction(const std::string &library_in, std::string &reaction_in,
                 const std::vector<double> &mg_e_bounds_in = {});

  //! Constructor (overridden gendir path)
  NDI_TNReaction(const std::string &gendir_in, const std::string &library_in,
                 std::string reaction_in, const std::vector<double> &mg_e_bounds_in = {});

  //! Default destructor
  ~NDI_TNReaction() = default;

  //! Disable default constructor
  NDI_TNReaction() = delete;

  //! Disable copy constructor
  NDI_TNReaction(const NDI_TNReaction &rhs) = delete;

  //! Disable move constructor
  NDI_TNReaction(NDI_TNReaction &&rhs) noexcept = delete;

  //! Disable assignment
  NDI_TNReaction &operator=(const NDI_TNReaction &rhs) = delete;

  //! Disable move assignment
  NDI_TNReaction &operator=(NDI_TNReaction &&rhs) noexcept = delete;

  //! Return spectrum PDF at a given temperature.
  [[nodiscard]] std::vector<double> get_PDF(const int product_zaid, const double temperature) const;
  [[nodiscard]] std::pair<const std::vector<double> &, std::vector<double>>
  get_ebins_and_PDF(const int product_zaid, const double temperature) const;

  //! Get the reaction
  [[nodiscard]] inline std::string get_reaction() const & { return reaction; }

  //! Get the name of the reaction from the NDI file
  [[nodiscard]] inline std::string get_reaction_name() const & { return reaction_name; }

  //! Get number of reaction products
  [[nodiscard]] inline uint32_t get_num_products() const {
    return static_cast<uint32_t>(products.size());
  }

  //! Get vector of reaction products
  [[nodiscard]] inline std::vector<int> get_products() const & { return products; }

  //! Get vector of reaction product multiplicities
  [[nodiscard]] inline std::vector<int> get_product_multiplicities() const & {
    return product_multiplicities;
  }

  //! Get vector of reaction temperature grid support points
  [[nodiscard]] inline std::vector<double> get_reaction_temperature() const & {
    return reaction_temperature;
  }

  //! Get vector of incident energy support points
  //! TODO: Rm? Is this actually needed for our use case?
  [[nodiscard]] inline std::vector<double> get_einbar() const & { return einbar; }

  //! Get vector of cross section support points
  //! TODO: Rm? Is this actually needed for our use case?
  [[nodiscard]] inline std::vector<double> get_sigvbar() const & { return sigvbar; }

  //! Get change in energy due to reaction
  [[nodiscard]] inline double get_reaction_q() const { return q_reaction; }

  //! Get number of groups for product
  [[nodiscard]] inline int get_num_groups(const int pzaid) const {
    Require(product_zaid_to_index.count(pzaid) == 1);
    const int pidx = product_zaid_to_index.at(pzaid);
    Check(std::cmp_less(pidx, prod_mg_e_bounds.size()));
    return num_groups[pidx];
  }

  //! Get group boundaries for product (MeV)
  [[nodiscard]] inline std::vector<double> get_mg_e_bounds(const int pzaid) const & {
    Require(product_zaid_to_index.count(pzaid) == 1);
    const int pidx = product_zaid_to_index.at(pzaid);
    Check(std::cmp_less(pidx, prod_mg_e_bounds.size()));
    return prod_mg_e_bounds[pidx];
  }

  //! Get group energies for product (MeV)
  [[nodiscard]] inline std::vector<double> get_group_energies(const int pzaid) const & {
    Require(product_zaid_to_index.count(pzaid) == 1);
    const int pidx = product_zaid_to_index.at(pzaid);
    Check(std::cmp_less(pidx, group_energies.size()));
    return group_energies[pidx];
  }

private:
  //! Flag to indicate if MG bounds were passed in on construction.
  bool bounds_passed{false};

  //! Parameters used to generate group bounds (if not passed to the constructor). These determine
  //! the portion of the outgoing energy CDF that we want to resolve for the lowest and central
  //! temperature distributions. The remainder of the bounds will come from the highest temperature
  //! data (widest in outgoing energy)
  static constexpr double lowT_CDF_fr{
      0.5}; // Ensure bounds resolve the middle 50% of the low-temp CDF
  static constexpr double midT_CDF_fr{
      0.8}; // Ensure bounds resolve the middle 80% of the mid-temp CDF

  //! Name of reaction to read
  const std::string reaction;

  //! Name of reaction as found in NDI data
  std::string reaction_name;

  //! Labels (ZAIDs) for reaction products
  std::vector<int> products;

  //! Map from reaction product ZAID to index
  std::map<int, int> product_zaid_to_index;

  //! Multiplicities for each reaction product
  std::vector<int> product_multiplicities;

  //! Temperature support point grid for reaction (keV)
  std::vector<double> reaction_temperature;

  //! Incident energy support point grid for reaction (MeV)
  std::vector<double> einbar;

  //! Incident cross section support point grid for reaction (cm^3 sh^-1)
  std::vector<double> sigvbar;

  //! Temperature support point grids for each reaction product (keV)
  std::vector<std::vector<double>> product_temperatures;

  //! Distribution support point grids for each reaction product
  std::vector<std::vector<std::vector<double>>> product_distributions;

  //! Reaction Q value i.e. change in energy
  double q_reaction = 0.0;

  //! Number of groups per TN product
  std::vector<uint32_t> num_groups;

  //! Group average energies (MeV)
  std::vector<std::vector<double>> group_energies;

  //! Energy bounds of multigroup data (MeV) per product, user-input or determined internally
  std::vector<std::vector<double>> prod_mg_e_bounds;

// Only implemented if NDI is found
#ifdef NDI_FOUND
  void load_ndi(const std::vector<double> &input_kev_bounds);
  std::vector<double> compute_mg_bounds_for_product(const int &ndi_handle, const size_t idx);
#endif
};

} // namespace rtt_cdi_ndi

#endif // cdi_ndi_NDI_TNReaction_hh

//------------------------------------------------------------------------------------------------//
// End cdi_ndi/NDI_TNReaction.hh
//------------------------------------------------------------------------------------------------//
