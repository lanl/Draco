//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/NDI_TNReaction.cc
 * \author Ben R. Ryan
 * \brief  NDI_TNReaction member definitions.
 * \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi_ndi/NDI_TNReaction.hh"
#include "dsxx/Assert.hh"
#include "cdi_ndi/NDI_Base.hh"
#include <string>
#include <utility> // move
#include <vector>
#ifdef NDI_FOUND
#include "dsxx/Query_Env.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/dbc.hh"
#include <algorithm> // upper_bound, count
#include <array>
#include <cstddef> // size_t
#include <cstdint> // uint32_t
#include <limits>
#include <numeric>
#include <sstream>
#endif
#ifdef DRACO_DUMP_TN_DATA
#include <fstream>
#endif

// cSpell:ignore edist, einbar, mlt, reac, rprods, sigvbars

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
                               std::string reaction_in, const std::vector<double> &mg_e_bounds_in)
    : NDI_Base(gendir_in, "tn", library_in), reaction(std::move(reaction_in)) {

#ifdef NDI_FOUND
  Require(reaction.length() > 0);

  if (mg_e_bounds_in.size() > 0) {
    bounds_passed = true;
  }

  load_ndi(mg_e_bounds_in);
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
                               const std::vector<double> &mg_e_bounds_in)
#ifdef NDI_FOUND
    : NDI_TNReaction(rtt_dsxx::get_env_val<std::string>("NDI_GENDIR_PATH").second, library_in,
                     reaction_in, mg_e_bounds_in)
#else
    : NDI_Base("tn", library_in), reaction(std::move(reaction_in))
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
 * \param[in] input_kev_bounds (optional) input MG bounds for NDI data collapse 
 *
 * This function opens an NDI file, navigates to the appropriate data, reads the data into internal
 * buffers, and closes the file. For more details on NDI, see
 * https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html
 *
 * Don't warn about ndi specific variables used in this function.
 * NOLINTBEGIN(misc-include-cleaner)
 */
void NDI_TNReaction::load_ndi(const std::vector<double> &input_kev_bounds) {
  int gendir_handle = -1;
  int dataset_handle = -1;
  int ndi_error = -9999;
  constexpr int c_str_len = 4096;
  std::array<char, c_str_len> c_str_buf{' '};

  // Open gendir file (index of a complete NDI dataset)
  ndi_error = NDI2_open_gendir(&gendir_handle, gendir.c_str());
  ndi_throw_rt_error(ndi_error);
  Insist(gendir_handle != -1, "gendir_handle still has default value!");

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

  //! Set reaction
  ndi_error = NDI2_set_reaction(dataset_handle, reaction.c_str());
  ndi_throw_rt_error(ndi_error);

  //! Store reaction name from NDI file
  ndi_error = NDI2_get_string_val(dataset_handle, NDI_ZAID, c_str_buf.data(), c_str_len);
  ndi_throw_rt_error(ndi_error);
  reaction_name = std::string(c_str_buf.data());

  //! Get number of temperature support points for reaction
  int num_temps;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_TEMPS, &num_temps);
  ndi_throw_rt_error(ndi_error);
  Require(num_temps > 1);
  reaction_temperature.resize(num_temps);

  //! Get temperature support points for reaction
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_TEMPS, reaction_temperature.data(),
                                   static_cast<int>(reaction_temperature.size()));
  ndi_throw_rt_error(ndi_error);
  // MeV -> keV
  for (auto &temperature : reaction_temperature) {
    temperature *= 1000.;
  }

  //! Get number of incident energy support points for reaction
  int const num_einbar = NDI2_get_size(dataset_handle, NDI_EINBAR, &ndi_error);
  ndi_throw_rt_error(ndi_error);
  einbar.resize(num_einbar);

  //! Get incident energy support points for reaction
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_EINBAR, einbar.data(),
                                   static_cast<int>(einbar.size()));
  ndi_throw_rt_error(ndi_error);

  //! Get the number of interp regions... for now just throw an exception if this is not equal to 1
  int const num_einbar_interp_regions =
      NDI2_get_size(dataset_handle, NDI_EINBAR_INTERP_REG, &ndi_error);
  ndi_throw_rt_error(ndi_error);
  Insist(num_einbar_interp_regions == 1, "Only 1 einbar interp region supported!");

  //! Get number of cross section support points for reaction
  int const num_sigvbar = NDI2_get_size(dataset_handle, NDI_SIGVBARS, &ndi_error);
  ndi_throw_rt_error(ndi_error);
  sigvbar.resize(num_sigvbar);

  //! Get cross section support points (cm^3 s-1 in NDI) for reaction
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_SIGVBARS, sigvbar.data(),
                                   static_cast<int>(sigvbar.size()));
  ndi_throw_rt_error(ndi_error);

  //! Get the number of interp regions... for now just throw an exception if this is not equal to 1
  int const num_sigvbar_interp_regions =
      NDI2_get_size(dataset_handle, NDI_SIGVBAR_INTERP_REG, &ndi_error);
  ndi_throw_rt_error(ndi_error);
  Insist(num_sigvbar_interp_regions == 1, "Only 1 sigvbar interp region supported!");

  //! Get number of reaction products
  int num_products;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_SEC_PARTS, &num_products);
  ndi_throw_rt_error(ndi_error);
  Require(num_products > 0);
  products.resize(num_products);
  product_temperatures.resize(num_products);
  prod_mg_e_bounds.resize(num_products);
  num_groups.resize(num_products);
  group_energies.resize(num_products);
  product_distributions.resize(num_products);
  product_multiplicities.resize(num_products);
  Require(num_products > 0);

  //! Get reaction product multiplicity
  ndi_error = NDI2_get_int_vec(dataset_handle, NDI_RPRODS_MLT, product_multiplicities.data(),
                               static_cast<int>(product_multiplicities.size()));
  ndi_throw_rt_error(ndi_error);

  //! Get change in energy due to reaction
  ndi_error = NDI2_get_float64_val(dataset_handle, NDI_REAC_Q, &q_reaction);
  ndi_throw_rt_error(ndi_error);
  q_reaction *= 1000.; // MeV -> keV

  std::vector<double> input_mev_bounds;
  if (bounds_passed) {
    Check(input_kev_bounds.size() > 0);
    // Make a copy of the input bounds so you can reverse/convert units.
    input_mev_bounds.resize(input_kev_bounds.size());
    std::ranges::transform(input_kev_bounds, input_mev_bounds.begin(),
                           [](double x) noexcept { return x / 1000.0; });

    // If passed, check that mg_e_bounds is monotonically decreasing (NDI requirement)
    Insist(
        rtt_dsxx::is_strict_monotonic_decreasing(input_mev_bounds.begin(), input_mev_bounds.end()),
        "Product multigroup bounds not strictly monotonic decreasing!");
    Insist(input_mev_bounds.back() > 0, "Negative product multigroup bounds!");

  } // else, bounds will be dynamically determined per-product later.

  //! Loop over reaction products
  for (int n = 0; n < num_products; n++) {
    //! Get ZAID of reaction product
    ndi_error =
        NDI2_get_string_val_n(dataset_handle, NDI_SEC_PART_TYPES, n, c_str_buf.data(), c_str_len);
    ndi_throw_rt_error(ndi_error);
    const std::string product_zaid = std::string(c_str_buf.data());

    // Ensure no duplicate products
    Require(std::count(products.begin(), products.end(), std::stoi(product_zaid)) == 0);

    products[n] = std::stoi(product_zaid);
    product_zaid_to_index.insert(std::pair<int, int>(products[n], n));

    //! Set NDI to reaction product
    ndi_error = NDI2_set_option(dataset_handle, NDI_CURR_PART, product_zaid.c_str());
    ndi_throw_rt_error(ndi_error);

    //! Get number of temperature support points (this can depend on reaction product)
    const int num_product_temps = NDI2_get_size(dataset_handle, NDI_EDIST_TEMPS, &ndi_error);
    ndi_throw_rt_error(ndi_error);
    Require(num_product_temps > 1);
    product_temperatures[n].resize(num_product_temps);
    product_distributions[n].resize(num_product_temps);

    //! Get temperature support points
    ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_TEMPS, product_temperatures[n].data(),
                                     static_cast<int>(product_temperatures[n].size()));
    ndi_throw_rt_error(ndi_error);

    //! Get the number of interp regions... for now just throw an exception if this is not equal to
    //! 1
    int const num_edist_interp_regions =
        NDI2_get_size(dataset_handle, NDI_EDIST_INTERP_REG, &ndi_error);
    ndi_throw_rt_error(ndi_error);
    Insist(num_edist_interp_regions == 1, "Only 1 edist interp region supported!");

    // If bounds aren't passed, determine them for this product; otherwise simply copy in passed bounds.
    if (bounds_passed) {
      prod_mg_e_bounds[n] = input_mev_bounds;
    } else {
      prod_mg_e_bounds[n] = compute_mg_bounds_for_product(dataset_handle, n);
    }
    num_groups[n] = static_cast<int>(prod_mg_e_bounds[n].size() - 1);

    // NOW use the dynamically-determined bounds to compute collapsed PDFs
    ndi_error =
        NDI2_set_float64_vec_option(dataset_handle, NDI_COLLAPSE, prod_mg_e_bounds[n].data(),
                                    static_cast<int>(prod_mg_e_bounds[n].size()));
    ndi_throw_rt_error(ndi_error);
    //! Check that number of groups is consistent
    int num_groups_int;
    ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_GRPS, &num_groups_int);
    Check(std::cmp_equal(num_groups_int, prod_mg_e_bounds[n].size() - 1));
    ndi_throw_rt_error(ndi_error);
    Require(num_groups_int > 0);
    group_energies[n].resize(num_groups[n]);

    //! Get average energies of energy groups (returned in MeV)
    ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_E_AVG, group_energies[n].data(),
                                     static_cast<int>(group_energies[n].size()));
    ndi_throw_rt_error(ndi_error);

    // Loop over temperatures
    for (size_t m = 0; m < product_temperatures[n].size(); m++) {
      std::ostringstream temp_stream;

      temp_stream << product_temperatures[n][m];
      ndi_error = NDI2_set_option(dataset_handle, NDI_TEMP, temp_stream.str().c_str());
      ndi_throw_rt_error(ndi_error);
      product_distributions[n][m].resize(num_groups[n]);
      ndi_error =
          NDI2_get_float64_vec(dataset_handle, NDI_EDIST, product_distributions[n][m].data(),
                               static_cast<int>(product_distributions[n][m].size()));
      ndi_throw_rt_error(ndi_error);
    } // END OF TEMPERATURE LOOP

#ifdef DRACO_DUMP_TN_DATA
    std::cout << "!!! Hey! We made it in here! " << std::endl;
    // Option to print out the product-dependent dynamically determined bounds
    // This data shares a single set of MG bounds, so the format is:
    //      Temp_1  Temp_2    Temp_3  ...
    // G_1  PDF_1,1  PDF_2,1  PDF_3,1 ...
    // G_2  PDF_1,2  PDF_2,2  PDF_3,2 ...
    // [...]
    // G_N  PDF_1,N  PDF_2,N  PDF_3,N ...
    std::ofstream datafile;
    datafile.open(std::string("fixed_group_bound_data_") + product_zaid);
    datafile << std::scientific << std::setprecision(4);
    datafile << "# T(MeV): ";
    for (size_t t = 0; t < product_temperatures[n].size(); t++) {
      datafile << product_temperatures[n][t] << "\t";
    }
    datafile << std::endl;

    for (size_t a = 0; a < prod_mg_e_bounds[n].size() - 1; a++) {
      datafile << 0.5 * (prod_mg_e_bounds[n][a] + prod_mg_e_bounds[n][a + 1]) << "\t";
      for (size_t t = 0; t < product_temperatures[n].size(); t++) {
        datafile << product_distributions[n][t][a] /
                        (prod_mg_e_bounds[n][a] - prod_mg_e_bounds[n][a + 1])
                 << "\t";
      }
      datafile << std::endl;
    }
    datafile.close();
#endif

    // Convert this product's temperatures from MeV -> keV
    std::transform(product_temperatures[n].begin(), product_temperatures[n].end(),
                   product_temperatures[n].begin(), [](double x) { return x * 1000.0; });
  }

  // Convert xs cm^3 s^-1 -> cm^3 sh^-1
  for (auto &sigma : sigvbar) {
    sigma *= 1.e-8;
  }

  //! Close datafile
  ndi_error = NDI2_close_gendir(gendir_handle);
  ndi_throw_rt_error(ndi_error);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Use a simple scheme to determine a single set of MG bounds across all temperature evals
 * \param[in] ndi_handle Handle to the current working NDI dataset
 * \param[in] prod_idx Index of product in product_temperatures array
 * \return Bounds computed to capture [most] salient features across temperature range
 */
std::vector<double> NDI_TNReaction::compute_mg_bounds_for_product(const int &ndi_handle,
                                                                  const size_t prod_idx) {

  Ensure(prod_idx < product_temperatures.size());
  Ensure(prod_idx < products.size());

  std::vector<double> final_bounds;

  const size_t num_temps(product_temperatures[prod_idx].size());
  //! Temp storage for energy support point grids for this reaction product
  std::vector<std::vector<double>> cont_energy_pts(num_temps);
  std::vector<std::vector<double>> cont_dist_pts(num_temps);
  // Gather all pointwise data for this product:
  int ndi_error = -9999;
  for (size_t m = 0; m < num_temps; m++) {
    std::ostringstream temp_stream;

    temp_stream << product_temperatures[prod_idx][m];
    ndi_error = NDI2_set_option(ndi_handle, NDI_TEMP, temp_stream.str().c_str());
    ndi_throw_rt_error(ndi_error);
    // Get number of energy out points for reaction
    const int num_eout = NDI2_get_size(ndi_handle, NDI_EDIST_CONT_E, &ndi_error);
    cont_energy_pts[m].resize(num_eout);
    cont_dist_pts[m].resize(num_eout);
    ndi_throw_rt_error(ndi_error);
    ndi_error = NDI2_get_float64_vec(ndi_handle, NDI_EDIST_CONT_E, cont_energy_pts[m].data(),
                                     static_cast<int>(cont_energy_pts[m].size()));
    ndi_throw_rt_error(ndi_error);
    ndi_error = NDI2_get_float64_vec(ndi_handle, NDI_EDIST_CONT_P, cont_dist_pts[m].data(),
                                     static_cast<int>(cont_dist_pts[m].size()));
    ndi_throw_rt_error(ndi_error);
  }

#ifdef DRACO_DUMP_TN_DATA
  // Option to print out full pointwise reference data from NDI file.
  // Because each temperature point has its own set of bounds, the format is:
  // Temp_1        Temp_2        Temp_3        ...
  // E_1,1 PDF_1,1 E_2,1 PDF_2,1 E_3,1 PDF_3,1 ...
  // E_1,2 PDF_1,2 E_2,2 PDF_2,2 E_3,2 PDF_3,2 ...
  // [...]
  // E_1,N PDF_1,N E_2,N PDF_2,N E_3,N PDF_3,N ...
  const auto zaid = products[prod_idx];
  std::ofstream datafile;
  datafile.open("full_tdep_bounds_reference_" + std::to_string(zaid));
  // Write the temperature values
  datafile << std::scientific << std::setprecision(4);
  datafile << "# T(MeV): ";
  for (size_t n = 0; n < num_temps; n++) {
    datafile << product_temperatures[prod_idx][n] << "\t\t";
  }
  datafile << std::endl;
  // Write the data in a gnuplottable format (alternating energy points and PDF
  // values for each temperature point)
  for (size_t d = 0; d < cont_energy_pts[0].size(); d++) {
    for (size_t n = 0; n < num_temps; n++) {
      Check(cont_energy_pts[0].size() == cont_energy_pts[n].size());
      datafile << cont_energy_pts[n][d] << " " << cont_dist_pts[n][d] << "\t";
    }
    datafile << std::endl;
  }
  datafile.close();
#endif

  // The goal here is to find a single (temperature-independent) set of bounds for this
  // reaction product. We do this by understanding that the doppler broadening with
  // increasing temperature always acts to widen and flatten the distribution. So, we
  // capture some fine bounds around the center of the "most-peaked" (lowest-temperature)
  // distribution, tack on some wider bounds for the shoulders using mid-temperature data
  // and finally round out the set by including high/low-E tails from the highest-temperature data.
  // This is simple, product-dependent, and easy to fine-tune, but it is likely overkill -
  // you could almost certainly coarsen the group structure at the tails for most parts of the
  // temperature grid.  You could also fix a max temperature well below 1 MeV for many problems, which
  // would improve sampling efficiency *and* save memory.
  final_bounds.reserve(2 * cont_energy_pts[0].size());

  // First, determine the energy bounds for the "middle 50%" of the lowest temperature
  {
    bool done{false};
    double cdf{cont_dist_pts[0][0]};
    size_t idx{1};
    const double cdf_norm = std::accumulate(cont_dist_pts[0].begin(), cont_dist_pts[0].end(), 0.0);
    const double low_cdf_val{cdf_norm * (0.5 - 0.5 * lowT_CDF_fr)};
    const double hi_cdf_val{cdf_norm * (0.5 + 0.5 * lowT_CDF_fr)};
    while (!done) {
      Insist(idx < cont_energy_pts[0].size() - 1,
             "Hit index outside of bounds in dynamic bounds determination");
      const double this_bound{0.5 * (cont_energy_pts[0][idx - 1] + cont_energy_pts[0][idx])};
      cdf += cont_dist_pts[0][idx];
      idx++;
      if (cdf >= low_cdf_val && cdf <= hi_cdf_val) {
        final_bounds.push_back(this_bound);
      } else if (cdf > 0.75 * cdf_norm) {
        done = true;
      }
    }
  }

  // Now, we've built in most of the fidelity we need for the peaked center of the distribution.
  // Next we'll take roughly the middle temperature point to capture the "shoulders"
  {
    bool done{false};
    const auto mid_t{static_cast<size_t>(std::floor(num_temps / 2))};
    double cdf{cont_dist_pts[mid_t][0]};
    size_t idx{1};
    const double cdf_norm =
        std::accumulate(cont_dist_pts[mid_t].begin(), cont_dist_pts[mid_t].end(), 0.0);
    const double low_cdf_val{cdf_norm * (0.5 - 0.5 * midT_CDF_fr)};
    const double hi_cdf_val{cdf_norm * (0.5 + 0.5 * midT_CDF_fr)};
    const auto [minit, maxit] = std::ranges::minmax_element(final_bounds);
    const double min{*minit};
    const double max{*maxit};
    while (!done) {
      Insist(idx < cont_energy_pts[mid_t].size() - 1,
             "Hit index outside of bounds in dynamic bounds determination");
      const double this_bound{0.5 *
                              (cont_energy_pts[mid_t][idx - 1] + cont_energy_pts[mid_t][idx])};
      cdf += cont_dist_pts[mid_t][idx];

      idx++;
      if ((cdf >= low_cdf_val && cdf <= hi_cdf_val) && (this_bound < min || this_bound > max)) {
        final_bounds.push_back(this_bound);
      } else if (cdf > 0.9 * cdf_norm) {
        done = true;
      }
    }
  }

  // Now,we've built in most of the fidelity we need for the center/shoulder of the distribution.
  // Finally, we consider the highest temperature point to build the higher/lower energy groups
  {
    const size_t last_t{cont_dist_pts.size() - 1};
    const auto [minit, maxit] = std::ranges::minmax_element(final_bounds);
    const double min{*minit};
    const double max{*maxit};
    // Push back the first (lowest) bound
    final_bounds.push_back(0.5 * cont_energy_pts[last_t][0]);
    for (size_t k = 1; k < cont_energy_pts[last_t].size(); k++) {
      const double this_bound{0.5 * (cont_energy_pts[last_t][k - 1] + cont_energy_pts[last_t][k])};
      if (this_bound < min || this_bound > max) {
        final_bounds.push_back(this_bound);
      } // else ignore it - we don't add bounds between the existing ones.
    }
    // Push back the last (highest) bound
    final_bounds.push_back(2.0 * cont_energy_pts[last_t][cont_energy_pts[last_t].size() - 1]);
    // Sort these into descending order:
    std::ranges::sort(final_bounds, std::greater<>());
  }
  Insist(rtt_dsxx::is_strict_monotonic_decreasing(final_bounds.begin(), final_bounds.end()),
         "Product multigroup bounds not strictly monotonic decreasing!");
  final_bounds.shrink_to_fit();

  return final_bounds;
}

// NOLINTEND(misc-include-cleaner)

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return normalized probability distribution function for energy of a reaction product at a
 *        given temperature, along with energy bin edges
 * \param[in] product_zaid ZAID of reaction product to sample
 * \param[in] temperature of plasma (keV)
 * \return Normalized PDF of reaction product energy
 */
std::pair<const std::vector<double> &, std::vector<double>>
NDI_TNReaction::get_ebins_and_PDF(const int product_zaid, const double temperature) const {

  Require(product_zaid_to_index.count(product_zaid) == 1);
  const int product_index = product_zaid_to_index.at(product_zaid);

  std::vector<double> pdf(num_groups[product_index]);

  Require(std::count(products.begin(), products.end(), product_zaid) == 1);

  Require(temperature > product_temperatures[product_index].front());
  Require(temperature < product_temperatures[product_index].back());

  auto temp_1 = std::upper_bound(product_temperatures[product_index].begin(),
                                 product_temperatures[product_index].end(), temperature);
  uint32_t const index_1 =
      static_cast<uint32_t>(temp_1 - product_temperatures[product_index].begin());
  uint32_t const index_0 = index_1 - 1;
  double const temp_0 = product_temperatures[product_index][index_0];
  Check(*temp_1 - temp_0 > std::numeric_limits<double>::min());
  double const fac = 1. - (temperature - temp_0) / (*temp_1 - temp_0);
  Check(fac >= 0 && fac <= 1);

  Remember(double dist_sum = 0.);
  for (uint32_t n = 0; n < num_groups[product_index]; n++) {
    double const pdf_0 = product_distributions[product_index][index_0][n];
    double const pdf_1 = product_distributions[product_index][index_1][n];
    pdf[n] = pdf_0 * fac + pdf_1 * (1. - fac);
    Remember(dist_sum += pdf[n]);
  }

  Require(rtt_dsxx::soft_equiv(dist_sum, 1., 1.e-10));

  // Careful here - std::pair and std::make_pair both use type decay
  // (discarding the ref on the first return value)
  return {prod_mg_e_bounds[product_index], pdf};
}

#else
std::pair<const std::vector<double> &, std::vector<double>>
NDI_TNReaction::get_ebins_and_PDF(const int /*product_zaid*/, const double /*temperature*/) const {
  Insist(0, "get_ebins_and_PDF() only available when NDI library is found");
}

#endif // NDI_FOUND

} // namespace rtt_cdi_ndi

//------------------------------------------------------------------------------------------------//
// End cdi_ndi/NDI_TNReaction.cc
//------------------------------------------------------------------------------------------------//
