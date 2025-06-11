//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/NDI_MultigroupNeutron.cc
 * \author Kyle S. Beling
 * \date   2025 March 18
 * \brief  NDI_MultigroupNeutron member definitions.
 * \note   Copyright (C) 2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi_ndi/NDI_MultigroupNeutron.hh"
#include "dsxx/Assert.hh"
#include "dsxx/dbc.hh"
#include "cdi_ndi/NDI_Base.hh"
#include "cdi_ndi/config.h"
#include <array>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace rtt_cdi_ndi {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for NDI reader specific to ZAID data with CMake-set gendir file.
 *
 * \param[in] library_in name of requested NDI data library
 * \param[in] target_zaid_in name of requested target ZAID
 * \param[in] mg_e_bounds_in energy boundaries of multigroup bins (MeV)
 */
#ifdef NDI_CP2011_FOUND
NDI_MultigroupNeutron::NDI_MultigroupNeutron(const std::string &library_in, int target_zaid_in,
                                             const std::vector<double> &mg_e_bounds_in)
    : NDI_Base(NDI_CP2011_GENDIR, "multigroup_neutron", library_in), target_zaid(target_zaid_in),
      coarse_group_bounds(mg_e_bounds_in) {

  Require(mg_e_bounds_in.size() > 0);

  // Check that mg_e_bounds_in is monotonically decreasing (NDI requirement)
  Insist(rtt_dsxx::is_strict_monotonic_decreasing(mg_e_bounds_in.begin(), mg_e_bounds_in.end()),
         "Product multigroup bounds not strictly monotonic decreasing!");
  Insist(mg_e_bounds_in.back() > 0, "Negative product multigroup bounds!");

  load_ndi();
}
#else
NDI_MultigroupNeutron::NDI_MultigroupNeutron(const std::string &library_in, int target_zaid_in,
                                             const std::vector<double> &mg_e_bounds_in)
    : NDI_Base("multigroup_neutron", library_in), target_zaid(target_zaid_in),
      coarse_group_bounds(mg_e_bounds_in) {
  Insist(false, "NDI_MultigroupNeutron objects require CP2011 data!");
}
#endif

/*!
 * \brief Constructor for NDI reader specific to ZAID data with provided path to gendir file.
 *
 * \param[in] gendir_in path to gendir file
 * \param[in] library_in name of requested NDI data library
 * \param[in] target_zaid_in name of requested target ZAID
 * \param[in] mg_e_bounds_in energy boundaries of multigroup bins (MeV)
 */
NDI_MultigroupNeutron::NDI_MultigroupNeutron(const std::string &gendir_in,
                                             const std::string &library_in, int target_zaid_in,
                                             const std::vector<double> &mg_e_bounds_in)
    : NDI_Base(gendir_in, "multigroup_neutron", library_in), target_zaid(target_zaid_in),
      coarse_group_bounds(mg_e_bounds_in) {

#ifdef NDI_FOUND
  Require(mg_e_bounds_in.size() > 0);

  // Check that mg_e_bounds_in is monotonically decreasing (NDI requirement)
  Insist(rtt_dsxx::is_strict_monotonic_decreasing(mg_e_bounds_in.begin(), mg_e_bounds_in.end()),
         "Product multigroup bounds not strictly monotonic decreasing!");
  Insist(mg_e_bounds_in.back() > 0, "Negative product multigroup bounds!");

  load_ndi();
#endif
}

// Protect actual NDI calls with NDI_FOUND macro:
#ifdef NDI_FOUND

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Load NDI dataset. Split off from constructor to allow for both default and overridden
 *        gendir paths.
 *
 * This function opens an NDI file, navigates to the appropriate data, reads the data into internal
 * buffers, and closes the file. For more details on NDI, see
 * https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html
 *
 * Don't warn about ndi specific variables used in this function.
 * NOLINTBEGIN(misc-include-cleaner)
 */
void NDI_MultigroupNeutron::load_ndi() {
  int gendir_handle = -1;
  int dataset_handle = -1;
  int ndi_error = -9999;
  constexpr int c_str_len = 4096;
  std::array<char, c_str_len> c_str_buf{' '};

  // Open gendir file (index of a complete NDI dataset)
  ndi_error = NDI2_open_gendir(&gendir_handle, gendir.c_str());
  Insist(ndi_error == 0, "Error when opening gendir file" + gendir);

  // Set dataset option by changing default value for this handle
  ndi_error = NDI2_set_option_gendir(gendir_handle, NDI_LIB_TYPE_DEFAULT, dataset.c_str());
  ndi_throw_rt_error(ndi_error);

  //! Set library option by changing default value for this handle
  ndi_error = NDI2_set_option_gendir(gendir_handle, NDI_LIBRARY_DEFAULT, library.c_str());
  ndi_throw_rt_error(ndi_error);

  //! Get dataset handle
  ndi_error = NDI2_make_handle(gendir_handle, dataset.c_str(), &dataset_handle);
  ndi_throw_rt_error(ndi_error);
  Insist(dataset_handle != -1, "dataset_handle still has default value!");

  //! Set ZAID
  ndi_error = NDI2_set_zaid(dataset_handle, std::to_string(target_zaid).c_str());
  ndi_throw_rt_error(ndi_error);

  //! Get group bounds and average energies
  // =========================================================================================== //
  int num_fine_groups = 0;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_GRPS, &num_fine_groups);
  ndi_throw_rt_error(ndi_error);
  Require(num_fine_groups > 0);

  std::vector<double> fine_group_bounds(num_fine_groups + 1);
  ndi_error =
      NDI2_get_float64_vec(dataset_handle, NDI_E_BOUNDS, fine_group_bounds.data(), num_fine_groups);
  ndi_throw_rt_error(ndi_error);

  // Remap coarse groups to fine groups
  coarse_group_bounds = remapCoarseOntoFine(coarse_group_bounds, fine_group_bounds);

  // Group collapse
  ndi_error = NDI2_set_float64_vec_option(dataset_handle, NDI_COLLAPSE, coarse_group_bounds.data(),
                                          static_cast<int>(coarse_group_bounds.size()));
  ndi_throw_rt_error(ndi_error);

  //! Get total cross section
  // =========================================================================================== //
  const int num_sigtot = NDI2_get_size(dataset_handle, NDI_SIG_TOT, &ndi_error);
  ndi_throw_rt_error(ndi_error);

  sigma_tot.resize(num_sigtot);
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_SIG_TOT, sigma_tot.data(), num_sigtot);
  ndi_throw_rt_error(ndi_error);

  //! Get group velocities
  // =========================================================================================== //
  const int num_vel = NDI2_get_size(dataset_handle, NDI_VEL, &ndi_error);
  ndi_throw_rt_error(ndi_error);

  group_velocities.resize(num_vel);
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_VEL, group_velocities.data(), num_vel);
  ndi_throw_rt_error(ndi_error);

  //! Get total kermas
  // =========================================================================================== //
  int num_kermas = NDI2_get_size(dataset_handle, NDI_KERMA, &ndi_error);
  ndi_throw_rt_error(ndi_error);

  group_kermas.resize(num_kermas);
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_KERMA, group_kermas.data(), num_kermas);
  ndi_throw_rt_error(ndi_error);

  //! Get total scattering matrix
  // =========================================================================================== //
  int num_legendre = 0;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_PN_ORDER, &num_legendre);
  ndi_throw_rt_error(ndi_error);

  std::vector<double> vector_temp;
  total_scattering_matrix.resize(num_legendre);
  for (int n = 0; n < num_legendre; ++n) {
    const int pn_size = NDI2_get_size_n(dataset_handle, NDI_PN_FULL, n, &ndi_error);
    ndi_throw_rt_error(ndi_error);

    vector_temp.resize(pn_size);
    ndi_error = NDI2_get_float64_vec_n(dataset_handle, NDI_PN_FULL, n, vector_temp.data(), pn_size);
    ndi_throw_rt_error(ndi_error);

    total_scattering_matrix[n] = vectorToSquareMatrix(vector_temp, coarse_group_bounds.size() - 1);
  }

  //! Get reaction info
  // =========================================================================================== //
  int num_reac = 0;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_REAC, &num_reac);
  ndi_throw_rt_error(ndi_error);

  // Resize reaction vectors
  reaction_mts.resize(num_reac);
  reaction_q_val.resize(num_reac);
  reaction_total_sigmas.resize(num_reac);
  reaction_products.resize(num_reac);
  reaction_multiplicities.resize(num_reac);

  // Loop through reactions
  for (int r = 0; r < num_reac; ++r) {
    ndi_error = NDI2_get_int_val_n(dataset_handle, NDI_REAC_MT, r, &reaction_mts[r]);
    ndi_throw_rt_error(ndi_error);

    ndi_error =
        NDI2_get_float64_val_n(dataset_handle, NDI_REAC_Q, reaction_mts[r], &reaction_q_val[r]);
    ndi_throw_rt_error(ndi_error);

    // Get products and their multiplicites
    const int prod_size =
        NDI2_get_size_n(dataset_handle, NDI_RPROD_ID_ALL, reaction_mts[r], &ndi_error);
    ndi_throw_rt_error(ndi_error);

    reaction_products[r].resize(prod_size);
    ndi_error = NDI2_get_int_vec_n(dataset_handle, NDI_RPROD_ID_ALL, reaction_mts[r],
                                   reaction_products[r].data(), prod_size);
    ndi_throw_rt_error(ndi_error);

    reaction_multiplicities[r].resize(prod_size);
    ndi_error = NDI2_get_int_vec_n(dataset_handle, NDI_RPROD_MLT_ALL, reaction_mts[r],
                                   reaction_multiplicities[r].data(), prod_size);
    ndi_throw_rt_error(ndi_error);

    // Get total cross section
    const int vec_size = NDI2_get_size_n(dataset_handle, NDI_SIG_REAC, reaction_mts[r], &ndi_error);
    ndi_throw_rt_error(ndi_error);

    reaction_total_sigmas[r].resize(vec_size);
    ndi_error = NDI2_get_float64_vec_n(dataset_handle, NDI_SIG_REAC, reaction_mts[r],
                                       reaction_total_sigmas[r].data(), vec_size);
    ndi_throw_rt_error(ndi_error);
  }

  //! Get number of reaction products
  int num_products = 0;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_SEC_PARTS, &num_products);
  ndi_throw_rt_error(ndi_error);

  product_list.resize(num_products);
  product_kermas.resize(num_products);
  product_scattering_matrix.resize(num_products);

  //! Loop over reaction products
  // =========================================================================================== //
  for (int n = 0; n < num_products; n++) {
    //! Get ZAID of reaction product
    ndi_error =
        NDI2_get_string_val_n(dataset_handle, NDI_SEC_PART_TYPES, n, c_str_buf.data(), c_str_len);
    ndi_throw_rt_error(ndi_error);
    const std::string product_zaid = std::string(c_str_buf.data());
    product_list[n] = std::stoi(product_zaid);

    //! Set NDI to reaction product
    ndi_error = NDI2_set_option(dataset_handle, NDI_CURR_PART, product_zaid.c_str());
    ndi_throw_rt_error(ndi_error);

    //! Get energy grid for reaction product
    int num_fine_groups_prod = 0;
    ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_GRPS_CURR_PART, &num_fine_groups_prod);
    ndi_throw_rt_error(ndi_error);

    std::vector<double> fine_bounds_prod(num_fine_groups_prod + 1);
    ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_E_BOUNDS_CURR_PART,
                                     fine_bounds_prod.data(), num_fine_groups_prod + 1);
    ndi_throw_rt_error(ndi_error);

    // Remap coarse groups to fine groups for product
    std::vector<double> new_group_bounds =
        remapCoarseOntoFine(coarse_group_bounds, fine_bounds_prod);

    //! Group collapse product
    ndi_error =
        NDI2_set_float64_vec_option(dataset_handle, NDI_E_BOUNDS_CURR_PART, new_group_bounds.data(),
                                    static_cast<int>(new_group_bounds.size()));
    ndi_throw_rt_error(ndi_error);

    //! Get number of kermas
    num_kermas = NDI2_get_size(dataset_handle, NDI_KERMA_CURR_PART, &ndi_error);
    ndi_throw_rt_error(ndi_error);

    product_kermas[n].resize(num_kermas);
    ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_KERMA_CURR_PART, product_kermas[n].data(),
                                     num_kermas);
    ndi_throw_rt_error(ndi_error);

    //! Get number of Legendre moments
    ndi_error = NDI2_get_int_val(dataset_handle, NDI_PN_ORDER_CURR_PART, &num_legendre);
    ndi_throw_rt_error(ndi_error);

    //! Get scattering matrix
    product_scattering_matrix[n].resize(num_legendre);
    for (int l = 0; l < num_legendre; ++l) {
      const int data_size =
          NDI2_get_size_n(dataset_handle, NDI_PN_PROD_FULL_CURR_PART, l, &ndi_error);
      ndi_throw_rt_error(ndi_error);

      std::vector<double> temp_data(data_size);
      ndi_error = NDI2_get_float64_vec_n(dataset_handle, NDI_PN_PROD_FULL_CURR_PART, l,
                                         temp_data.data(), data_size);
      ndi_throw_rt_error(ndi_error);

      product_scattering_matrix[n][l] =
          vectorToSquareMatrix(temp_data, new_group_bounds.size() - 1);
    }
  }

  //! Close datafile
  // =========================================================================================== //
  ndi_error = NDI2_close_gendir(gendir_handle);
  ndi_throw_rt_error(ndi_error);
}

// NOLINTEND(misc-include-cleaner)

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the total cross section for a neutron incident energy.
 * \param[in] incident_energy of neutron (MeV)
 * \return Total cross section value
 */
double NDI_MultigroupNeutron::totalCrossSection(const double incident_energy) const {
  const int index = find_bin(this->coarse_group_bounds, incident_energy);
  return this->sigma_tot[index];
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return normalized probability distribution function for a incident neutron energy
 * \param[in] incident_energy of neutron (MeV)
 * \return Normalized PDF of reaction product energy
 */
std::vector<double> NDI_MultigroupNeutron::pdf(const double incident_energy) const {
  const int index = find_bin(this->coarse_group_bounds, incident_energy);
  std::vector<double> distribution = total_scattering_matrix[0][index];
  const double sum = std::accumulate(distribution.begin(), distribution.end(), 0.0);
  if (sum > std::numeric_limits<double>::epsilon()) {
    for (double &value : distribution)
      value /= sum;
  }
  return distribution;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the total cross section for a specific reaction.
 * \param[in] mt_number MT number of the reaction
 * \param[in] incident_energy of neutron (MeV)
 * \return Total cross section value
 */
double NDI_MultigroupNeutron::totalCrossSection(const int mt_number,
                                                const double incident_energy) const {
  const int mt_index = find_reaction_index(mt_number);
  const int index = find_bin(this->coarse_group_bounds, incident_energy);
  return this->reaction_total_sigmas[mt_index][index];
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return normalized probability distribution function of a reaction product for an incident 
 *        neutron energy.
 * \param[in] product_zaid ZAID of reaction product to sample
 * \param[in] incident_energy of neutron (MeV)
 * \return Normalized PDF of reaction product energy
 */
std::vector<double> NDI_MultigroupNeutron::pdf(const int product_zaid,
                                               const double incident_energy) const {
  const int prod_index = find_product_index(product_zaid);
  const int index = find_bin(this->coarse_group_bounds, incident_energy);
  std::vector<double> distribution = product_scattering_matrix[prod_index][0][index];
  const double sum = std::accumulate(distribution.begin(), distribution.end(), 0.0);
  if (sum > std::numeric_limits<double>::epsilon()) {
    for (double &value : distribution)
      value /= sum;
  }
  return distribution;
}

#else

double NDI_MultigroupNeutron::totalCrossSection(const double /*incident_energy*/) const {
  Insist(0, "totalCrossSection() only available when NDI library is found");
}

double NDI_MultigroupNeutron::totalCrossSection(const int /*product_zaid*/,
                                                const double /*incident_energy*/) const {
  Insist(0, "totalCrossSection() only available when NDI library is found");
}

std::vector<double> NDI_MultigroupNeutron::pdf(const double /*incident_energy*/) const {
  Insist(0, "pdf() only available when NDI library is found");
}

std::vector<double> NDI_MultigroupNeutron::pdf(const int /*product_zaid*/,
                                               const double /*incident_energy*/) const {
  Insist(0, "pdf() only available when NDI library is found");
}

#endif

} // namespace rtt_cdi_ndi

//------------------------------------------------------------------------------------------------//
// End cdi_ndi/NDI_MultigroupNeutron.cc
//------------------------------------------------------------------------------------------------//
