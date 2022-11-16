//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/NDI_CP_Eloss.cc
 * \author Ben R. Ryan
 * \date   2020 Jun 3
 * \brief  NDI_CP_Eloss member definitions.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "NDI_CP_Eloss.hh"
#include "ds++/DracoStrings.hh"
#include "ds++/Interpolate.hh"
#include "ds++/Soft_Equivalence.hh"

#include <cmath>

namespace rtt_cdi_ndi {

// Protect actual NDI calls with NDI_FOUND macro:
#ifdef NDI_FOUND
//----------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------//
/*!
 * \brief Constructor for NDI reader specific to TN DEDX data with provided path to gendir file.
 *
 * \param[in] gendir_in path to gendir file
 * \param[in] library_in name of requested NDI data library
 * \param[in] target_in Target CParticle object
 * \param[in] projectile_in Projectile CParticle object
 */
NDI_CP_Eloss::NDI_CP_Eloss(const std::string &gendir_in, const std::string &library_in,
                           rtt_cdi::CParticle target_in, rtt_cdi::CParticle projectile_in)
    : rtt_cdi::CPEloss(target_in, projectile_in, rtt_cdi::CPModelType::TABULAR_ETYPE,
                       rtt_cdi::CPModelAngleCutoff::NONE),
      NDI_Base(gendir_in, "dedx", library_in) {
  // TODO BRR set angle cutoff by parsing library_in

  load_ndi();
}
/*!
 * \brief Constructor for NDI reader specific to TN DEDX data using default gendir file.
 *
 * \param[in] library_in name of requested NDI data library
 * \param[in] target_in Target CParticle object
 * \param[in] projectile_in Projectile CParticle object
 */
NDI_CP_Eloss::NDI_CP_Eloss(const std::string &library_in, rtt_cdi::CParticle target_in,
                           rtt_cdi::CParticle projectile_in)
    : rtt_cdi::CPEloss(target_in, projectile_in, rtt_cdi::CPModelType::TABULAR_ETYPE,
                       rtt_cdi::CPModelAngleCutoff::NONE),
      NDI_Base("dedx", library_in) {
  // TODO BRR set angle cutoff by parsing library_in

  load_ndi();
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Load NDI dataset
 *
 * This function opens an NDI file, navigates to the appropriate dataset, reads the data into
 * internal buffers, and closes the file. For more details on NDI, see
 * https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html
 */
void NDI_CP_Eloss::load_ndi() {
#if NDI_DEDX_SUPPORT
  int gendir_handle = -1;
  int dataset_handle = -1;
  int ndi_error = -9999;

  // Open gendir file (index of a complete NDI dataset)
  ndi_error = NDI2_open_gendir(&gendir_handle, gendir.c_str());
  Insist(ndi_error == 0, "Error when opening gendir file" + gendir);

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

  //! Set projectile isotope
  ndi_error = NDI2_set_isotope(dataset_handle, std::to_string(projectile.get_zaid()).c_str());
  Require(ndi_error == 0);

  int num_targets = 0;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_TARGET, &num_targets);
  Require(ndi_error == 0);

  std::vector<int> target_zaids(num_targets);
  ndi_error = NDI2_get_int_vec(dataset_handle, NDI_TARGET_ZAID, target_zaids.data(),
                               static_cast<int>(target_zaids.size()));
  Require(ndi_error == 0);

  int num_grps = 0;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_GRPS, &num_grps);
  Require(ndi_error == 0);
  n_energy = static_cast<uint32_t>(num_grps);

  energies.resize(n_energy);
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_ENERGIES, energies.data(),
                                   static_cast<int>(energies.size()));
  Require(ndi_error == 0);
  min_log_energy = energies.front();
  d_log_energy = energies[1] - energies[0];
  min_energy = exp(min_log_energy);
  max_energy = exp(min_log_energy + d_log_energy * n_energy);

  int num_densities = 0;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_DENSITIES, &num_densities);
  Require(ndi_error == 0);
  n_density = static_cast<uint32_t>(num_densities);

  densities.resize(n_density);
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_DENSITIES, densities.data(),
                                   static_cast<int>(densities.size()));
  Require(ndi_error == 0);
  min_log_density = densities.front();
  d_log_density = densities[1] - densities[0];
  min_density = target.get_mass() * exp(min_log_density);
  max_density = target.get_mass() * exp(min_log_density + d_log_density * n_density);

  int num_temperatures = 0;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_TEMPS, &num_temperatures);
  Require(ndi_error == 0);
  n_temperature = static_cast<uint32_t>(num_temperatures);

  temperatures.resize(n_temperature);
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_TEMPS, temperatures.data(),
                                   static_cast<int>(temperatures.size()));
  Require(ndi_error == 0);
  min_log_temperature = temperatures.front();
  d_log_temperature = temperatures[1] - temperatures[0];
  min_temperature = exp(min_log_temperature);
  max_temperature = exp(min_log_temperature + d_log_temperature * n_temperature);

  stopping_data_1d.resize(n_energy * n_density * n_temperature);
  if (target.get_zaid() == -1) {
    ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_DEDX, stopping_data_1d.data(),
                                     static_cast<int>(stopping_data_1d.size()));
  } else {
    ndi_error = NDI2_get_float64_vec_x(
        dataset_handle, NDI_TARGET_DEDX, std::to_string(target.get_zaid()).c_str(),
        stopping_data_1d.data(), static_cast<int>(stopping_data_1d.size()));
  }
  Require(ndi_error == 0);

  // Check for uniform log spacing
  for (uint32_t n = 1; n < n_energy; n++) {
    Require(rtt_dsxx::soft_equiv(d_log_energy, energies[n] - energies[n - 1], 1.e-4));
  }
  for (uint32_t n = 1; n < n_density; n++) {
    Require(rtt_dsxx::soft_equiv(d_log_density, densities[n] - densities[n - 1], 1.e-4));
  }
  for (uint32_t n = 1; n < n_temperature; n++) {
    Require(rtt_dsxx::soft_equiv(d_log_temperature, temperatures[n] - temperatures[n - 1], 1.e-4));
  }

  // Convert units on table to match those of getEloss:
  //   energy:      MeV -> cm/shk (using target particle mass)
  const double energy_cgs = exp(min_log_energy) * (1.e6 * pc.electronVolt());
  min_log_energy = log(sqrt(2. * energy_cgs / target.get_mass()) * 1.e-8);
  d_log_energy = d_log_energy / 2.;
  //   density:     cm^-3 -> g cm^-3
  min_log_density = log(exp(min_log_density) * target.get_mass());
  //   temperature: keV -> keV
  // Note that d log x = dx / x is not affected by unit conversion factors
  for (auto &energy : energies) {
    energy = sqrt(2. * (exp(energy) * 1.e6 * pc.electronVolt()) / target.get_mass()) * 1.e-8;
  }
  for (auto &density : densities) {
    density = exp(density) * target.get_mass();
  }
  for (auto &temperature : temperatures) {
    temperature = exp(temperature);
  }

  //! Close datafile
  ndi_error = NDI2_close_gendir(gendir_handle);
  Require(ndi_error == 0);
#else
  Insist(0,
         "NDI version " + std::string(NDI_VERSION_STRING) + " does not support stopping powers!");
#endif // NDI_DEDX_SUPPORT
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate the tabulated stopping power for a given material and projectile state.
 * \param[in] temperature Material temperature [keV]
 * \param[in] density Material density [g cm^-3]
 * \param[in] partSpeed Particle speed [cm shk^-1]
 */
double NDI_CP_Eloss::getEloss(const double temperature, const double density,
                              const double partSpeed) const {
  if (temperature <= min_temperature || temperature >= max_temperature || density <= min_density ||
      density >= max_density || partSpeed <= min_energy || partSpeed >= max_energy) {
    // Outside of the table
    return 0.;
  }

  const int pt0_energy =
      static_cast<int>(std::floor((log(partSpeed) - min_log_energy) / d_log_energy));
  const int pt1_energy = pt0_energy + 1;
  const int pt0_density =
      static_cast<int>(std::floor((log(density) - min_log_density) / d_log_density));
  const int pt1_density = pt0_density + 1;
  const int pt0_temperature =
      static_cast<int>(std::floor((log(temperature) - min_log_temperature) / d_log_temperature));
  const int pt1_temperature = pt0_temperature + 1;

  const double x0 = exp(min_log_energy + pt0_energy * d_log_energy);
  const double x1 = exp(min_log_energy + pt1_energy * d_log_energy);
  const double y0 = exp(min_log_density + pt0_density * d_log_density);
  const double y1 = exp(min_log_density + pt1_density * d_log_density);
  const double z0 = exp(min_log_temperature + pt0_temperature * d_log_temperature);
  const double z1 = exp(min_log_temperature + pt1_temperature * d_log_temperature);

  const double f000 = get_stopping_data(pt0_energy, pt0_density, pt0_temperature);
  const double f100 = get_stopping_data(pt1_energy, pt0_density, pt0_temperature);
  const double f001 = get_stopping_data(pt0_energy, pt0_density, pt1_temperature);
  const double f101 = get_stopping_data(pt1_energy, pt0_density, pt1_temperature);
  const double f010 = get_stopping_data(pt0_energy, pt1_density, pt0_temperature);
  const double f110 = get_stopping_data(pt1_energy, pt1_density, pt0_temperature);
  const double f011 = get_stopping_data(pt0_energy, pt1_density, pt1_temperature);
  const double f111 = get_stopping_data(pt1_energy, pt1_density, pt1_temperature);
  const double dedx =
      exp(rtt_dsxx::interpolate::linear_3d(x0, x1, y0, y1, z0, z1, f000, f100, f001, f101, f010,
                                           f110, f011, f111, partSpeed, density, temperature));
  const double number_density = density / target.get_mass();
  return dedx * 1000. * number_density * partSpeed; // MeV cm^2 -> keV shk^-1
}

#else

//----------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------//

/*!
 * \brief Constructor for NDI reader specific to TN DEDX data with provided path to gendir file.
 *
 * \param[in] gendir_in path to gendir file
 * \param[in] library_in name of requested NDI data library
 * \param[in] target_in name of requested reaction
 * \param[in] projectile_in energy boundaries of multigroup bins (keV)
 */
NDI_CP_Eloss::NDI_CP_Eloss(const std::string &gendir_in, const std::string &library_in,
                           rtt_cdi::CParticle target_in, rtt_cdi::CParticle projectile_in)
    : rtt_cdi::CPEloss(target_in, projectile_in, rtt_cdi::CPModelType::TABULAR_ETYPE,
                       rtt_cdi::CPModelAngleCutoff::NONE),
      NDI_Base(gendir_in, "dedx", library_in) {}
/*!
 * \brief Constructor for NDI reader specific to TN DEDX data using default gendir file.
 *
 * \param[in] library_in name of requested NDI data library
 * \param[in] target_in name of requested reaction
 * \param[in] projectile_in energy boundaries of multigroup bins (keV)
 */
NDI_CP_Eloss::NDI_CP_Eloss(const std::string &library_in, rtt_cdi::CParticle target_in,
                           rtt_cdi::CParticle projectile_in)
    : rtt_cdi::CPEloss(target_in, projectile_in, rtt_cdi::CPModelType::TABULAR_ETYPE,
                       rtt_cdi::CPModelAngleCutoff::NONE),
      NDI_Base("dedx", library_in) {}

//----------------------------------------------------------------------------//
/*!
 * \brief Interpolate the tabulated stopping power for a given material and projectile state.
 *
 * \param[in] temperature Material temperature [keV]
 * \param[in] density Material density [g cm^-3]
 * \param[in] partSpeed Particle speed [cm shk^-1]
 */
double NDI_CP_Eloss::getEloss(const double /*temperature*/, const double /*density*/,
                              const double /*partSpeed*/) const {
  Insist(0, "getEloss only defined if NDI_FOUND!");
}

#endif // NDI_FOUND
} // namespace rtt_cdi_ndi

//----------------------------------------------------------------------------//
// End cdi_ndi/NDI_CP_Eloss.cc
//----------------------------------------------------------------------------//
