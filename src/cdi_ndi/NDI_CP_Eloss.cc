//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ndi/NDI_CP_Eloss.cc
 * \author Ben R. Ryan
 * \date   2020 Jun 3
 * \brief  NDI_CP_Eloss member definitions.
 * \note   Copyright (C) 2020-2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi_ndi/NDI_CP_Eloss.hh"
#include "dsxx/Assert.hh"
#include "dsxx/Interpolate.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "cdi/CPCommon.hh"
#include "cdi/CPEloss.hh"
#include "cdi_ndi/NDI_Base.hh"
#include <algorithm>
#include <cmath>
#include <string>

// cSpell:ignore Eloss, dedx, zaids

namespace rtt_cdi_ndi {

// Protect actual NDI calls with NDI_FOUND macro:
#ifdef NDI_FOUND

//------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for NDI reader specific to TN DEDX data with provided path to gendir file.
 *
 * \param[in] gendir_in path to gendir file
 * \param[in] library_in name of requested NDI data library
 * \param[in] zaidt_ target ZAID
 * \param[in] mt_ target mass (g)
 * \param[in] zaidp_ projectile ZAID
 * \param[in] mp_ projectile mass (g)
 */
NDI_CP_Eloss::NDI_CP_Eloss(const std::string &gendir_in, const std::string &library_in,
                           const int zaidt_, const double mt_, const int zaidp_, const double mp_)
    : rtt_cdi::CPEloss(rtt_cdi::CPModelType::TABULAR_ETYPE, rtt_cdi::CPModelAngleCutoff::NONE),
      NDI_Base(gendir_in, "dedx", library_in), zaidt(zaidt_), mt(mt_), zaidp(zaidp_), mp(mp_) {
  // TODO BRR set angle cutoff by parsing library_in

  load_ndi();
}
/*!
 * \brief Constructor for NDI reader specific to TN DEDX data using default gendir file.
 *
 * \param[in] library_in name of requested NDI data library
 * \param[in] zaidt_ target ZAID
 * \param[in] mt_ target mass (g)
 * \param[in] zaidp_ projectile ZAID
 * \param[in] mp_ projectile mass (g)
 */
NDI_CP_Eloss::NDI_CP_Eloss(const std::string &library_in, const int zaidt_, const double mt_,
                           const int zaidp_, const double mp_)
    : rtt_cdi::CPEloss(rtt_cdi::CPModelType::TABULAR_ETYPE, rtt_cdi::CPModelAngleCutoff::NONE),
      NDI_Base("dedx", library_in), zaidt(zaidt_), mt(mt_), zaidp(zaidp_), mp(mp_) {
  // TODO BRR set angle cutoff by parsing library_in

  load_ndi();
}

// Don't warn about directly including NDI header files.

// NOLINTBEGIN(misc-include-cleaner)

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
  ndi_error = NDI2_set_isotope(dataset_handle, std::to_string(zaidp).c_str());
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
  d_log_energy = energies[1] - energies[0];

  int num_densities = 0;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_DENSITIES, &num_densities);
  Require(ndi_error == 0);
  n_density = static_cast<uint32_t>(num_densities);

  densities.resize(n_density);
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_DENSITIES, densities.data(),
                                   static_cast<int>(densities.size()));
  Require(ndi_error == 0);
  d_log_density = densities[1] - densities[0];

  int num_temperatures = 0;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_TEMPS, &num_temperatures);
  Require(ndi_error == 0);
  n_temperature = static_cast<uint32_t>(num_temperatures);

  temperatures.resize(n_temperature);
  ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_TEMPS, temperatures.data(),
                                   static_cast<int>(temperatures.size()));
  Require(ndi_error == 0);
  d_log_temperature = temperatures[1] - temperatures[0];

  stopping_data_1d.resize(n_energy * n_density * n_temperature);
  if (zaidt == -1) {
    ndi_error = NDI2_get_float64_vec(dataset_handle, NDI_DEDX, stopping_data_1d.data(),
                                     static_cast<int>(stopping_data_1d.size()));
  } else {
    ndi_error =
        NDI2_get_float64_vec_x(dataset_handle, NDI_TARGET_DEDX, std::to_string(zaidt).c_str(),
                               stopping_data_1d.data(), static_cast<int>(stopping_data_1d.size()));
  }
  Require(ndi_error == 0);

  // Check for uniform log spacing
  // =========================================================================================== //
  for (uint32_t n = 1; n < n_energy; n++) {
    Require(rtt_dsxx::soft_equiv(d_log_energy, energies[n] - energies[n - 1], 1.e-4));
  }
  for (uint32_t n = 1; n < n_density; n++) {
    Require(rtt_dsxx::soft_equiv(d_log_density, densities[n] - densities[n - 1], 1.e-4));
  }
  for (uint32_t n = 1; n < n_temperature; n++) {
    Require(rtt_dsxx::soft_equiv(d_log_temperature, temperatures[n] - temperatures[n - 1], 1.e-4));
  }

  // Set minimum and maximum grid values
  // =========================================================================================== //
  min_energy = energies.front();
  max_energy = energies.back();

  min_density = densities.front();
  max_density = densities.back();

  min_temperature = temperatures.front();
  max_temperature = temperatures.back();

  // Close datafile
  // =========================================================================================== //
  ndi_error = NDI2_close_gendir(gendir_handle);
  Require(ndi_error == 0);
#else
  Insist(0,
         "NDI version " + std::string(NDI_VERSION_STRING) + " does not support stopping powers!");
#endif // NDI_DEDX_SUPPORT
}

// NOLINTEND(misc-include-cleaner)

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate the tabulated stopping power for a given material and projectile state.
 *
 * \todo Share code with Tabular_CP_Eloss::getEloss (appears to be identical interp)
 *
 * \param[in] temperature Ion temperature [keV]
 * \param[in] density Target number density [atoms cm^-3]
 * \param[in] partSpeed Particle speed [cm shk^-1]
 * \param[in] debye_length Plasma Debye length [cm], unused
 */
double NDI_CP_Eloss::getEloss(const double temperature, const double density,
                              const double partSpeed,
                              [[maybe_unused]] const double debye_length) const {

  // temperature: keV -> keV
  const double log_temperature = std::log(temperature);

  // density: cm^-3 -> cm^-3
  const double log_density = std::log(density);

  // energy: cm/shk -> MeV  (using projectile particle's mass)
  const double log_energy =
      std::log(std::pow(partSpeed / 1.e-8, 2) * mp / (2.0 * 1.e6 * pc.electronVolt()));

  if (log_temperature <= min_temperature || log_temperature >= max_temperature ||
      log_density <= min_density || log_density >= max_density || log_energy <= min_energy ||
      log_energy >= max_energy) {
    // Outside of the table
    return 0.;
  }

  // Get grid indices
  const int pt1_temperature =
      static_cast<int>(std::lower_bound(temperatures.begin(), temperatures.end(), log_temperature) -
                       temperatures.begin());
  const int pt0_temperature = pt1_temperature - 1;

  const int pt1_density = static_cast<int>(
      std::lower_bound(densities.begin(), densities.end(), log_density) - densities.begin());
  const int pt0_density = pt1_density - 1;

  const int pt1_energy = static_cast<int>(
      std::lower_bound(energies.begin(), energies.end(), log_energy) - energies.begin());
  const int pt0_energy = pt1_energy - 1;

  // Get grid values
  const double x0 = temperatures[pt0_temperature];
  const double x1 = temperatures[pt1_temperature];
  const double y0 = densities[pt0_density];
  const double y1 = densities[pt1_density];
  const double z0 = energies[pt0_energy];
  const double z1 = energies[pt1_energy];

  // Get function values
  const double f000 = get_stopping_data(pt0_temperature, pt0_density, pt0_energy);
  const double f100 = get_stopping_data(pt1_temperature, pt0_density, pt0_energy);
  const double f001 = get_stopping_data(pt0_temperature, pt0_density, pt1_energy);
  const double f101 = get_stopping_data(pt1_temperature, pt0_density, pt1_energy);
  const double f010 = get_stopping_data(pt0_temperature, pt1_density, pt0_energy);
  const double f110 = get_stopping_data(pt1_temperature, pt1_density, pt0_energy);
  const double f011 = get_stopping_data(pt0_temperature, pt1_density, pt1_energy);
  const double f111 = get_stopping_data(pt1_temperature, pt1_density, pt1_energy);
  const double dedx = std::exp(
      rtt_dsxx::interpolate::linear_3d(x0, x1, y0, y1, z0, z1, f000, f100, f001, f101, f010, f110,
                                       f011, f111, log_temperature, log_density, log_energy));
  return dedx * 1000. * density * partSpeed; // MeV cm^2 -> keV shk^-1
}

#else

//------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//------------------------------------------------------------------------------------------------//

/*!
 * \brief Constructor for NDI reader specific to TN DEDX data with provided path to gendir file.
 *
 * \param[in] gendir_in path to gendir file
 * \param[in] library_in name of requested NDI data library
 * \param[in] zaidt_ target ZAID
 * \param[in] mt_ target mass (g)
 * \param[in] zaidp_ projectile ZAID
 * \param[in] mp_ projectile mass (g)
 */
NDI_CP_Eloss::NDI_CP_Eloss(const std::string &gendir_in, const std::string &library_in,
                           const int /*zaidt_*/, const double /*mt_*/, const int /*zaidp_*/,
                           const double /*mp_*/)
    : rtt_cdi::CPEloss(rtt_cdi::CPModelType::TABULAR_ETYPE, rtt_cdi::CPModelAngleCutoff::NONE),
      NDI_Base(gendir_in, "dedx", library_in) {}
/*!
 * \brief Constructor for NDI reader specific to TN DEDX data using default gendir file.
 *
 * \param[in] library_in name of requested NDI data library
 * \param[in] zaidt_ target ZAID
 * \param[in] mt_ target mass (g)
 * \param[in] zaidp_ projectile ZAID
 * \param[in] mp_ projectile mass (g)
 */
NDI_CP_Eloss::NDI_CP_Eloss(const std::string &library_in, const int /*zaidt_*/,
                           const double /*mt_*/, const int /*zaidp_*/, const double /*mp_*/)
    : rtt_cdi::CPEloss(rtt_cdi::CPModelType::TABULAR_ETYPE, rtt_cdi::CPModelAngleCutoff::NONE),
      NDI_Base("dedx", library_in) {}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate the tabulated stopping power for a given material and projectile state.
 *
 * \param[in] temperature Material temperature [keV]
 * \param[in] density Material density [g cm^-3]
 * \param[in] partSpeed Particle speed [cm shk^-1]
 * \param[in] debye_length Plasma Debye length [cm], unused
 */
double NDI_CP_Eloss::getEloss([[maybe_unused]] const double temperature,
                              [[maybe_unused]] const double density,
                              [[maybe_unused]] const double partSpeed,
                              [[maybe_unused]] const double debye_length) const {
  Insist(0, "getEloss only defined if NDI_FOUND!");
}

#endif // NDI_FOUND
} // namespace rtt_cdi_ndi

//------------------------------------------------------------------------------------------------//
// End cdi_ndi/NDI_CP_Eloss.cc
//------------------------------------------------------------------------------------------------//
