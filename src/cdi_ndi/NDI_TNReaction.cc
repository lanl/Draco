//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/NDI_TNReaction.cc
 * \author Ben R. Ryan
 * \date   2020 Feb 4
 * \brief  NDI_TNReaction member definitions.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "NDI_TNReaction.hh"
#include "ds++/Query_Env.hh"
#include "ds++/dbc.hh"
#include <array>
#include <cmath>

namespace rtt_cdi_ndi {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//------------------------------------------------------------------------------------------------//

/*!
 * \brief Constructor for NDI reader specific to TN reaction data with provided path to gendir file.
 *
 * \param[in] gendir_in path to gendir file
 * \param[in] library_in name of requested NDI data library
 * \param[in] reaction_in name of requested reaction
 * \param[in] mg_e_bounds_in energy boundaries of multigroup bins (keV)
 */
NDI_TNReaction::NDI_TNReaction(const std::string &gendir_in, const std::string &library_in,
                               std::string reaction_in, std::vector<double> mg_e_bounds_in)
    : NDI_Base(gendir_in, "tn", library_in), reaction(std::move(reaction_in)),
      mg_e_bounds(std::move(mg_e_bounds_in)) {

#ifdef NDI_FOUND
  Require(reaction.length() > 0);
  Require(mg_e_bounds.size() > 0);

  for (double &eb : mg_e_bounds)
    eb /= 1000.0; // keV -> MeV

  // Check that mg_e_bounds is monotonically decreasing (NDI requirement)
  Insist(rtt_dsxx::is_strict_monotonic_decreasing(mg_e_bounds.begin(), mg_e_bounds.end()),
         "Product multigroup bounds not strictly monotonic decreasing!");
  Insist(mg_e_bounds.back() > 0, "Negative product multigroup bounds!");

  load_ndi();
#endif
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for NDI reader specific to TN reaction data using default gendir file.
 *
 * \param[in] library_in name of requested NDI data library
 * \param[in] reaction_in name of requested reaction
 * \param[in] mg_e_bounds_in energy boundaries of multigroup bins (keV)
 *
 * This is a delegating constructor. It defaults the gendir filepath to the value found in the
 * environment.
 */
NDI_TNReaction::NDI_TNReaction(const std::string &library_in, std::string &reaction_in,
                               std::vector<double> mg_e_bounds_in)
#ifdef NDI_FOUND
    : NDI_TNReaction(rtt_dsxx::get_env_val<std::string>("NDI_GENDIR_PATH").second, library_in,
                     reaction_in, std::move(mg_e_bounds_in))
#else
    : NDI_Base("tn", library_in), reaction(std::move(reaction_in)),
      mg_e_bounds(std::move(mg_e_bounds_in))
#endif
{
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
 */
void NDI_TNReaction::load_ndi() {
  int gendir_handle = -1;
  int dataset_handle = -1;
  int ndi_error = -9999;
  constexpr int c_str_len = 4096;
  std::array<char, c_str_len> c_str_buf{' '};

  // Open gendir file (index of a complete NDI dataset)
  ndi_error = NDI2_open_gendir(&gendir_handle, gendir.c_str());
  Require(ndi_error == 0);
  Insist(gendir_handle != -1, "gendir_handle still has default value!");

  // Set dataset option by changing default value for this handle
  ndi_error = NDI2_set_option_gendir(gendir_handle, NDI_LIB_TYPE_DEFAULT, dataset.c_str());
  Require(ndi_error == 0);

  //! Set library option by changing default value for this handle
  ndi_error = NDI2_set_option_gendir(gendir_handle, NDI_LIBRARY_DEFAULT, library.c_str());
  Require(ndi_error == 0);

  //! Get dataset handle
  ndi_error = NDI2_make_handle(gendir_handle, dataset.c_str(), &dataset_handle);
  Require(ndi_error == 0);
  Insist(dataset_handle != -1, "dataset_handle still has default value!");

  //! Set reaction
  ndi_error = NDI2_set_reaction(dataset_handle, reaction.c_str());
  Require(ndi_error == 0);

  //! Store reaction name from NDI file
  ndi_error = NDI2_get_string_val(dataset_handle, NDI_ZAID, c_str_buf.data(), c_str_len);
  Require(ndi_error == 0);
  reaction_name = std::string(c_str_buf.data());

  //! Get number of temperature support points for reaction
  int num_temps;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_TEMPS, &num_temps);
  Require(ndi_error == 0);
  Require(num_temps > 1);
  reaction_temperature.resize(num_temps);

  //! Get temperature support points for reaction
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_TEMPS, reaction_temperature.data(),
                                   static_cast<int>(reaction_temperature.size()));
  Require(ndi_error == 0);
  // MeV -> keV
  for (auto &temperature : reaction_temperature) {
    temperature *= 1000.;
  }

  //! Get number of incident energy support points for reaction
  int num_einbar = NDI2_get_size(dataset_handle, NDI_EINBAR, &ndi_error);
  Require(ndi_error == 0);
  einbar.resize(num_einbar);

  //! Get incident energy support points for reaction
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_EINBAR, einbar.data(),
                                   static_cast<int>(einbar.size()));
  Require(ndi_error == 0);
  // MeV -> keV
  for (auto &energy : einbar) {
    energy *= 1000.;
  }

  //! Get the number of interp regions... for now just throw an exception if this is not equal to 1
  int num_einbar_interp_regions = NDI2_get_size(dataset_handle, NDI_EINBAR_INTERP_REG, &ndi_error);
  Require(ndi_error == 0);
  Insist(num_einbar_interp_regions == 1, "Only 1 einbar interp region supported!");

  //! Get number of cross section support points for reaction
  int num_sigvbar = NDI2_get_size(dataset_handle, NDI_SIGVBARS, &ndi_error);
  Require(ndi_error == 0);
  sigvbar.resize(num_sigvbar);

  //! Get cross section support points (cm^3 s-1 in NDI) for reaction
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_SIGVBARS, sigvbar.data(),
                                   static_cast<int>(sigvbar.size()));
  Require(ndi_error == 0);
  // cm^3 s^-1 -> cm^3 sh^-1
  for (auto &sigma : sigvbar) {
    sigma *= 1.e-8;
  }

  //! Get the number of interp regions... for now just throw an exception if this is not equal to 1
  int num_sigvbar_interp_regions =
      NDI2_get_size(dataset_handle, NDI_SIGVBAR_INTERP_REG, &ndi_error);
  Require(ndi_error == 0);
  Insist(num_sigvbar_interp_regions == 1, "Only 1 sigvbar interp region supported!");

  //! Get number of reaction products
  int num_products;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_SEC_PARTS, &num_products);
  Require(ndi_error == 0);
  Require(num_products > 0);
  products.resize(num_products);
  product_temperatures.resize(num_products);
  product_distributions.resize(num_products);
  product_multiplicities.resize(num_products);
  Require(num_products > 0);

  //! Get reaction product multiplicity
  ndi_error = NDI2_get_int_vec(dataset_handle, NDI_RPRODS_MLT, product_multiplicities.data(),
                               static_cast<int>(product_multiplicities.size()));
  Require(ndi_error == 0);

  //! Get change in energy due to reaction
  ndi_error = NDI2_get_float64_val(dataset_handle, NDI_REAC_Q, &q_reaction);
  Require(ndi_error == 0);
  q_reaction *= 1000.; // MeV -> keV

  //! Specify multigroup option
  ndi_error = NDI2_set_float64_vec_option(dataset_handle, NDI_COLLAPSE, mg_e_bounds.data(),
                                          static_cast<int>(mg_e_bounds.size()));
  Require(ndi_error == 0);

  //! Get number of groups
  int num_groups_int;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_GRPS, &num_groups_int);
  Require(ndi_error == 0);
  Require(num_groups_int > 0);
  num_groups = static_cast<uint32_t>(num_groups_int);
  group_bounds.resize(num_groups + 1);
  group_energies.resize(num_groups);

  //! Get boundaries of energy groups
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_E_BOUNDS, group_bounds.data(),
                                   static_cast<int>(group_bounds.size()));
  Require(ndi_error == 0);
  // MeV -> keV
  for (auto &bound : group_bounds) {
    bound *= 1000;
  }

  //! Get average energies of energy groups
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_E_AVG, group_energies.data(),
                                   static_cast<int>(group_energies.size()));
  Require(ndi_error == 0);
  // MeV -> keV
  for (auto &energy : group_energies) {
    energy *= 1000;
  }

  //! Loop over reaction products
  for (int n = 0; n < num_products; n++) {
    //! Get ZAID of reaction product
    ndi_error =
        NDI2_get_string_val_n(dataset_handle, NDI_SEC_PART_TYPES, n, c_str_buf.data(), c_str_len);
    Require(ndi_error == 0);
    const std::string product_zaid = std::string(c_str_buf.data());

    // Ensure no duplicate products
    Require(std::count(products.begin(), products.end(), std::stoi(product_zaid)) == 0);

    products[n] = std::stoi(product_zaid);
    product_zaid_to_index.insert(std::pair<int, int>(products[n], n));

    //! Set NDI to reaction product
    ndi_error = NDI2_set_option(dataset_handle, NDI_CURR_PART, product_zaid.c_str());
    Require(ndi_error == 0);

    //! Get number of temperature support points (this can depend on reaction product)
    const int num_product_temps = NDI2_get_size(dataset_handle, NDI_EDIST_TEMPS, &ndi_error);
    Require(ndi_error == 0);
    Require(num_product_temps > 1);
    product_temperatures[n].resize(num_product_temps);
    product_distributions[n].resize(num_product_temps);

    //! Get temperature support points
    ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_TEMPS, product_temperatures[n].data(),
                                     static_cast<int>(product_temperatures[n].size()));
    Require(ndi_error == 0);
    // MeV -> keV
    for (auto &temperature : product_temperatures[n]) {
      temperature *= 1000.;
    }

    //! Get the number of interp regions... for now just throw an exception if this is not equal to
    //! 1
    int num_edist_interp_regions = NDI2_get_size(dataset_handle, NDI_EDIST_INTERP_REG, &ndi_error);
    Require(ndi_error == 0);
    Insist(num_edist_interp_regions == 1, "Only 1 edist interp region supported!");

    // Loop over temperatures
    for (size_t m = 0; m < product_temperatures[n].size(); m++) {
      std::ostringstream temp_stream;
      temp_stream << product_temperatures[n][m] / 1000; // keV -> MeV
      ndi_error = NDI2_set_option(dataset_handle, NDI_TEMP, temp_stream.str().c_str());
      Require(ndi_error == 0);

      product_distributions[n][m].resize(num_groups);
      ndi_error =
          NDI2_get_float64_vec(dataset_handle, NDI_EDIST, product_distributions[n][m].data(),
                               static_cast<int>(product_distributions[n][m].size()));
      Require(ndi_error == 0);
    }
  }

  //! Close datafile
  ndi_error = NDI2_close_gendir(gendir_handle);
  Require(ndi_error == 0);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return normalized probability distribution function for energy of a reaction product at a
 *        given temperature.
 * \param[in] product_zaid ZAID of reaction product to sample
 * \param[in] temperature of plasma (keV)
 * \return Normalized PDF of reaction product energy
 */
std::vector<double> NDI_TNReaction::get_PDF(const int product_zaid,
                                            const double temperature) const {
  std::vector<double> pdf(num_groups);

  Require(std::count(products.begin(), products.end(), product_zaid) == 1);

  const int product_index = (*(product_zaid_to_index.find(product_zaid))).second;

  Require(temperature > product_temperatures[product_index].front());
  Require(temperature < product_temperatures[product_index].back());

  auto temp_1 = std::upper_bound(product_temperatures[product_index].begin(),
                                 product_temperatures[product_index].end(), temperature);
  uint32_t index_1 = static_cast<uint32_t>(temp_1 - product_temperatures[product_index].begin());
  uint32_t index_0 = index_1 - 1;
  double temp_0 = product_temperatures[product_index][index_0];
  Check(*temp_1 - temp_0 > std::numeric_limits<double>::min());
  double fac = 1. - (temperature - temp_0) / (*temp_1 - temp_0);
  Check(fac >= 0 && fac <= 1);

  Remember(double dist_sum = 0.);
  for (uint32_t n = 0; n < num_groups; n++) {
    double pdf_0 = product_distributions[product_index][index_0][n];
    double pdf_1 = product_distributions[product_index][index_1][n];
    pdf[n] = pdf_0 * fac + pdf_1 * (1. - fac);
    Remember(dist_sum += pdf[n]);
  }

  Require(rtt_dsxx::soft_equiv(dist_sum, 1., 1.e-10));

  return pdf;
}

#else

std::vector<double> NDI_TNReaction::get_PDF(const int /*product_zaid*/,
                                            const double /*temperature*/) const {
  Insist(0, "get_PDF() only available when NDI library is found");
}

#endif // NDI_FOUND

} // namespace rtt_cdi_ndi

//------------------------------------------------------------------------------------------------//
// End cdi_ndi/NDI_TNReaction.cc
//------------------------------------------------------------------------------------------------//
