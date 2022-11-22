//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/NDI_TNReaction.hh
 * \author Ben R. Ryan
 * \date   2019 Nov 4
 * \brief  NDI_TNReaction class definition.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef cdi_ndi_NDI_TNReaction_hh
#define cdi_ndi_NDI_TNReaction_hh

#include "NDI_Base.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include <algorithm>
#include <memory>
#include <string>

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
                 std::vector<double> mg_e_bounds_in);

  //! Constructor (overridden gendir path)
  NDI_TNReaction(const std::string &gendir_in, const std::string &library_in,
                 std::string reaction_in, std::vector<double> mg_e_bounds_in);

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

  //! Return spectrum PDF at a given temperature
  std::vector<double> get_PDF(const int product_zaid, const double temperature) const;

  //! Get the reaction
  inline std::string get_reaction() const & { return reaction; }

  //! Get the name of the reaction from the NDI file
  inline std::string get_reaction_name() const & { return reaction_name; }

  //! Get number of reaction products
  inline uint32_t get_num_products() const { return static_cast<uint32_t>(products.size()); }

  //! Get vector of reaction products
  inline std::vector<int> get_products() const & { return products; }

  //! Get vector of reaction product multiplicities
  inline std::vector<int> get_product_multiplicities() const & { return product_multiplicities; }

  //! Get vector of reaction temperature grid support points
  inline std::vector<double> get_reaction_temperature() const & { return reaction_temperature; }

  //! Get vector of incident energy support points
  inline std::vector<double> get_einbar() const & { return einbar; }

  //! Get vector of cross section support points
  inline std::vector<double> get_sigvbar() const & { return sigvbar; }

  //! Get change in energy due to reaction
  inline double get_reaction_q() const { return q_reaction; }

  //! Get number of groups
  inline int get_num_groups() const { return num_groups; }

  //! Get group boundaries
  inline std::vector<double> get_group_bounds() const & { return group_bounds; }

  //! Get group energies
  inline std::vector<double> get_group_energies() const & { return group_energies; }

private:
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

  //! Incident energy support point grid for reaction (keV)
  std::vector<double> einbar;

  //! Incident cross section support point grid for reaction (cm^3 sh^-1)
  std::vector<double> sigvbar;

  //! Temperature support point grids for each reaction product (keV)
  std::vector<std::vector<double>> product_temperatures;

  //! Distribution support point grids for each reaction product
  std::vector<std::vector<std::vector<double>>> product_distributions;

  //! Reaction Q value i.e. change in energy
  double q_reaction = 0.0;

  //! Number of groups
  uint32_t num_groups = 0;

  //! Group boundaries (keV)
  std::vector<double> group_bounds;

  //! Group average energies (keV)
  std::vector<double> group_energies;

  //! Energy bounds of multigroup data (MeV) to be passed to NDI
  std::vector<double> mg_e_bounds;

// Only implemented if NDI is found
#ifdef NDI_FOUND
  void load_ndi();
#endif
};

} // namespace rtt_cdi_ndi

#endif // cdi_ndi_NDI_TNReaction_hh

//------------------------------------------------------------------------------------------------//
// End cdi_ndi/NDI_TNReaction.hh
//------------------------------------------------------------------------------------------------//
