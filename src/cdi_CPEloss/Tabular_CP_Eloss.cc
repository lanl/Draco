//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/Tabular_CP_Eloss.cc
 * \author Ben R. Ryan
 * \date   2019 Nov 4
 * \brief  Tabular_CP_Eloss member definitions.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Tabular_CP_Eloss.hh"
#include "ds++/DracoStrings.hh"
#include "ds++/Interpolate.hh"
#include <utility>

namespace rtt_cdi_cpeloss {

using std::experimental::basic_mdspan;
using std::experimental::dynamic_extent;
using std::experimental::extents;
using std::experimental::layout_left;

//------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for an analytic tabular eloss model.
 *
 * This constructor builds an eloss model defined by the rtt_cdi_cpeloss::Tabular_Eloss_Model
 * derived class argument.
 *
 * The path to an eloss datafile is passed to the constructor, which opens and parses the file. The
 * file format is the usual LANL format for stopping powers.
 *
 * \param[in] filename_in path to eloss file
 * \param[in] target_in target particle zaid
 * \param[in] projectile_in transporting particle zaid
 * \param[in] model_angle_cutoff_in rtt_cdi::CPModelAngleCutoff the angle separating the stopping
 *                 power approximation from analog scattering
 */
Tabular_CP_Eloss::Tabular_CP_Eloss(std::string filename_in, rtt_cdi::CParticle target_in,
                                   rtt_cdi::CParticle projectile_in,
                                   rtt_cdi::CPModelAngleCutoff model_angle_cutoff_in)
    : rtt_cdi::CPEloss(target_in, projectile_in, rtt_cdi::CPModelType::TABULAR_ETYPE,
                       model_angle_cutoff_in),
      filename(std::move(filename_in)) {
  using std::stod;
  using std::stoi;

  file.open(filename);
  Insist(file.is_open(), "Error opening DEDX file \"" + filename + "\"");

  constexpr uint32_t max_entries =
      6; // This is a statement about the file format, maximum of six entries per row

  std::vector<std::string> line_entries = read_line(); // ZAID
  int32_t projectile_zaid_file = stoi(line_entries[0]);
  Insist(projectile.get_zaid() == projectile_zaid_file,
         "File projectile ZAID does not match constructor argument!");

  skip_lines(1); // Z, A, mass

  line_entries = read_line(); // Number of bins for energy, density, temperature
  Check(line_entries.size() >= 3);
  n_energy = stoi(line_entries[0]);
  n_density = stoi(line_entries[1]);
  n_temperature = stoi(line_entries[2]);

  line_entries = read_line(); // Bin spacing for energy, density, temperature (log)
  Check(line_entries.size() >= 3);
  d_log_energy = stod(line_entries[0]);
  d_log_density = stod(line_entries[1]);
  d_log_temperature = stod(line_entries[2]);
  Require(d_log_energy > 0);
  Require(d_log_density > 0);
  Require(d_log_temperature > 0);
  d_log_energy = 1. / d_log_energy;
  d_log_density = 1. / d_log_density;
  d_log_temperature = 1. / d_log_temperature;

  // Get first energy support point
  uint32_t nlines = (n_energy + max_entries - 1) / max_entries;
  line_entries = read_line();
  min_log_energy = stod(line_entries[0]);
  skip_lines(nlines - 1);
  energies.resize(n_energy);
  for (uint32_t n = 0; n < n_energy; n++) {
    energies[n] = exp(min_log_energy + n * d_log_energy);
  }

  // Get first density support point
  nlines = (n_density + max_entries - 1) / max_entries;
  line_entries = read_line();
  min_log_density = stod(line_entries[0]);
  skip_lines(nlines - 1);
  densities.resize(n_density);
  for (uint32_t n = 0; n < n_density; n++) {
    densities[n] = exp(min_log_density + n * d_log_density);
  }

  // Get first temperature support point
  nlines = (n_temperature + max_entries - 1) / max_entries;
  line_entries = read_line();
  min_log_temperature = stod(line_entries[0]);
  skip_lines(nlines - 1);
  temperatures.resize(n_temperature);
  for (uint32_t n = 0; n < n_temperature; n++) {
    temperatures[n] = exp(min_log_temperature + n * d_log_temperature);
  }

  stopping_data_1d.resize(n_energy * n_density * n_temperature);

  bool target_found = false;
  nlines = (n_energy * n_density * n_temperature + max_entries - 1) /
           max_entries; // The number of lines taken up by stopping power data for one target
  if (target.get_zaid() == -1) {
    // Target is free electrons
    target_found = true;
    uint32_t nentry = 0;
    for (uint32_t n = 0; n < nlines; n++) {
      line_entries = read_line();
      for (std::string const &entry : line_entries) {
        stopping_data_1d[nentry++] = stod(entry);
      }
    }
  } else {
    // Skip electrons
    skip_lines(nlines);

    // Find ion target, if it exists
    const uint32_t n_target_ions = stoi(read_line()[0]); // Number of target ions in file
    for (uint32_t n_target_ion = 0; n_target_ion < n_target_ions; n_target_ion++) {
      int zaid_target_ion = stoi(read_line()[0]); // ZAID
      read_line();                                // Z, A, mass
      if (zaid_target_ion == target.get_zaid()) {
        // This is the requested target ion
        target_found = true;
        uint32_t nentry = 0;
        for (uint32_t n = 0; n < nlines; n++) {
          line_entries = read_line();
          for (std::string const &entry : line_entries) {
            stopping_data_1d[nentry] = stod(entry);
            nentry++;
          }
        }
        break;
      } else {
        // This is not the requested target ion
        skip_lines(nlines);
      }
    }
  }
  file.close();

  Insist(target_found, "Error finding target ZAID \"" + std::to_string(target.get_zaid()) +
                           "\" in DEDX file \"" + filename + "\"");

  stopping_data =
      basic_mdspan<double, extents<dynamic_extent, dynamic_extent, dynamic_extent>, layout_left>(
          stopping_data_1d.data(), n_energy, n_density, n_temperature);

  // Convert units on table to match those of getEloss:
  //
  // energy:      MeV -> cm/shk (using target particle mass)
  double energy_cgs = exp(min_log_energy) * (1.e6 * pc.electronVolt());
  min_log_energy = log(sqrt(2. * energy_cgs / target.get_mass()) * 1.e-8);
  d_log_energy = d_log_energy / 2.;
  // density:     cm^-3 -> g cm^-3
  min_log_density = log(exp(min_log_density) * target.get_mass());
  // temperature: keV -> keV
  //
  // Note that d log x = dx / x is not affected by unit conversion factors
  for (auto &energy : energies) {
    energy = sqrt(2. * (energy * 1.e6 * pc.electronVolt()) / target.get_mass()) * 1.e-8;
  }
  for (auto &density : densities) {
    density *= target.get_mass();
  }

  // Initialize table bounds
  min_energy = exp(min_log_energy);
  max_energy = exp(min_log_energy + d_log_energy * (n_energy));
  min_density = exp(min_log_density);
  max_density = exp(min_log_density + d_log_density * (n_density));
  min_temperature = exp(min_log_temperature);
  max_temperature = exp(min_log_temperature + d_log_temperature * (n_temperature));
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Read through the next nlines lines and ignore them.
 * \param[in] nlines number of lines in file to skip
 */
void Tabular_CP_Eloss::skip_lines(uint32_t nlines) {
  std::string line;
  for (uint32_t nline = 0; nline < nlines; nline++) {
    std::getline(file, line);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Read a line from an eloss datafile and return as a vector of strings.
 * \return entries the resulting vector of entries in the datafile line.
 */
std::vector<std::string> const Tabular_CP_Eloss::read_line() {
  std::string line;
  std::getline(file, line);
  return rtt_dsxx::tokenize(line);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate the tabulated stopping power for a given material and projectile state.
 * \param[in] temperature Material temperature [keV]
 * \param[in] density Material density [g cm^-3]
 * \param[in] partSpeed Particle speed [cm shk^-1]
 */
double Tabular_CP_Eloss::getEloss(const double temperature, const double density,
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

  const double f000 = stopping_data(pt0_energy, pt0_density, pt0_temperature);
  const double f100 = stopping_data(pt1_energy, pt0_density, pt0_temperature);
  const double f001 = stopping_data(pt0_energy, pt0_density, pt1_temperature);
  const double f101 = stopping_data(pt1_energy, pt0_density, pt1_temperature);
  const double f010 = stopping_data(pt0_energy, pt1_density, pt0_temperature);
  const double f110 = stopping_data(pt1_energy, pt1_density, pt0_temperature);
  const double f011 = stopping_data(pt0_energy, pt1_density, pt1_temperature);
  const double f111 = stopping_data(pt1_energy, pt1_density, pt1_temperature);
  const double dedx =
      exp(rtt_dsxx::interpolate::linear_3d(x0, x1, y0, y1, z0, z1, f000, f100, f001, f101, f010,
                                           f110, f011, f111, partSpeed, density, temperature));
  const double number_density = density / target.get_mass();
  return dedx * 1000. * number_density * partSpeed; // MeV cm^2 -> keV shk^-1
}

} // namespace rtt_cdi_cpeloss

//------------------------------------------------------------------------------------------------//
// End cdi_CPEloss/Tabular_CP_Eloss.cc
//------------------------------------------------------------------------------------------------//
