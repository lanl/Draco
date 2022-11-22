//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   compton_tools/Compton_Edep.cc
 * \author Kendra Long
 * \date   10 October 2022
 * \brief  Implementation file for native Compton binary-read and temperature interpolation
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Compton_Edep.hh"
#include "ds++/dbc.hh"
#include "units/PhysicalConstexprs.hh"
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

// Empty namespace contains local helper functions
namespace {

// Define a small nudge to keep frequencies off of exact interp points:
constexpr double fuzz = std::numeric_limits<double>::epsilon();

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Extract a fixed number of data points from an fstream object
 *
 * \param[in] n_entries The number of data values to extract
 * \param[in] datafile An fstream to extract from
 * \return An n_entries-length vector of double data
 */
std::vector<double> get_data_line(size_t n_entries, std::ifstream &datafile) {

  std::vector<double> values(n_entries);

  std::string line;
  std::getline(datafile, line);
  std::stringstream datavals(line);
  for (size_t a = 0; a < n_entries; a++) {
    datavals >> values[a];
    // check for failure
    if (datavals.fail() || values[a] < 0) {
      Insist(0, "Failure during Edep data read!");
    }
  }

  return values;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Extract a single size_t value from a stream
 *
 * \param[in] sizeline A stream with size_t values in it
 * \return The next size_t value in the stream
 */
size_t get_data_size(std::stringstream &sizeline) {
  size_t val;
  sizeline >> val;
  // check for failure
  if (sizeline.fail() || val <= 0) {
    Insist(0, "Failure during Edep header data read!");
  }
  return val;
}
} // namespace

namespace compton_edep {

// read all multigroup LLNL library data (multiple libfiles)
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Read in an average energy deposition dataset and populate/return a simple data container
 *
 * \param[in] infile the (CSK-generated) base library filename
 * \return A data container with the raw sigma_c / EREC values
 *
 * \post the local reconstruction of the original data is returned.
 */
Edep_Container read_edep_file(const std::string &infile) {

  std::ifstream ascii_lib1;
  std::ifstream ascii_lib2;

  std::string sigc_file(infile + "_sigc");
  std::string erec_file(infile + "_erec");

  // try to open the base file with no additional qualifiers
  ascii_lib1.open(sigc_file);
  ascii_lib2.open(erec_file);

  // check file for "openness":
  if (!ascii_lib1.is_open()) {
    std::ostringstream msg;
    msg << "Failed to open " << sigc_file << "! ";
    Insist(0, msg.str());
  }
  // check file for "openness":
  if (!ascii_lib2.is_open()) {
    std::ostringstream msg;
    msg << "Failed to open " << erec_file << "! ";
    Insist(0, msg.str());
  }

  // Get the data from the first line, because it lets us size the return vector
  std::string line;
  getline(ascii_lib1, line);
  std::stringstream data_sizes(line);

  const size_t netempbp = get_data_size(data_sizes);
  const size_t netemp = get_data_size(data_sizes);
  const size_t nginbp = get_data_size(data_sizes);
  const size_t ngin = get_data_size(data_sizes);
  // Create and size the return container
  Edep_Container edepdata(netempbp, nginbp, netemp, ngin);

  // next line contains a list of electron temp breakpoints
  edepdata.T_break_pts = get_data_line(netempbp, ascii_lib1);
  edepdata.freq_break_pts = get_data_line(nginbp, ascii_lib1);
  edepdata.T_eval_pts = get_data_line(netemp, ascii_lib1);
  // discard corresponding header lines from erec file...
  std::string discard_line;
  getline(ascii_lib2, discard_line);
  getline(ascii_lib2, discard_line);
  getline(ascii_lib2, discard_line);
  getline(ascii_lib2, discard_line);

  std::string line1, line2;
  std::vector<int> curr_position((netempbp - 1) * (nginbp - 1), 0);
  for (size_t g = 0; g < ngin; g++) {
    // get region index for this frequency:
    const size_t f_index = static_cast<int>(g / edepdata.f_ppr);

    // grab a line from each file
    getline(ascii_lib1, line1);
    std::stringstream dataline1(line1);
    getline(ascii_lib2, line2);
    std::stringstream dataline2(line2);

    // stream in the gin value(s):
    double this_gin2;
    dataline1 >> edepdata.freq_eval_pts[g];
    dataline2 >> this_gin2;

    Insist(rtt_dsxx::soft_equiv(edepdata.freq_eval_pts[g], this_gin2),
           "Frequency evaluation points must be identical!");

    // stream in a value of erec/sigc for each electron temp:
    for (size_t k = 0; k < netemp; k++) {
      // get region index for this temperature:
      const size_t T_index = static_cast<int>(k / edepdata.T_ppr);
      const size_t region_index = f_index * (netempbp - 1) + T_index;
      dataline1 >> edepdata.sigc[region_index][curr_position[region_index]];
      dataline2 >> edepdata.erec[region_index][curr_position[region_index]];
      // increment the index for the region:
      curr_position[region_index]++;
    }
  }
  // Finally, init the temperature interpolation multipliers, using the data in the container:
  edepdata.compute_multipliers();

  return edepdata;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate Compton opacity (sigc) and Expected Relative Energy Change (EREC) data to a
 *        given frequency, using data preinterpolated in cell electron temperature.
 *
 * \param[in] cell  Cell index
 * \param[in] freq  Incident frequency (keV)
 * \param[in] data  Container of temperature-interpolated Lagrange poly coefficients by cell
 * \return    The interpolated sigma_c and EREC values
 */
std::pair<double, double> interpolate_sigc_erec(const size_t cell, const double freq,
                                                const Preinterp_Container &data) {
  Ensure(cell > 0 && cell <= data.erec_coeffs.size());
  rtt_units::PhysicalConstexprs<rtt_units::CGSH> consts;

  // Calculate electron rest-mass in keV:
  // units: me (g) * c^2 (cm^2 / shk^2) * (shk^2 / s^2) * (m^2 / cm^2) * (kg / g) * (keV / J)
  const double e_mocsq =
      consts.me() * consts.c() * consts.c() * (1e8) * (1e8) * (1e-3) * (1e-2) * (1e-2) * 6.2415e15;
  // Scale frequency by electron rest-mass
  double mygamma = freq / e_mocsq;

  // check it against the library (which should always have 2+ break points):
  Insist(data.freq_break_pts.size() > 1, "Missing CSK frequency break points");
  if (mygamma < data.freq_break_pts[0] || mygamma > data.freq_break_pts.back()) {
#if DRACO_DIAGNOSTICS & 3
    std::cout << "WARNING: Frequency " << mygamma << " is outside the bounds of the library ["
              << data.freq_break_pts[0] << ", " << data.freq_break_pts.back() << "]" << std::endl;
#endif
    if (mygamma < data.freq_break_pts[0]) {
      mygamma = data.freq_break_pts[0] * (1.0 + fuzz);
    }
    if (mygamma > data.freq_break_pts.back()) {
      mygamma = data.freq_break_pts.back() * (1.0 - fuzz);
    }
  }
  // compute the index of the frequency interpolation region:
  const auto f_index =
      std::lower_bound(data.freq_break_pts.begin(), data.freq_break_pts.end(), mygamma) -
      data.freq_break_pts.begin() - 1;

  // Collect all necessary data for this interpolation region.
  // TODO: Pass iterators directly?
  const auto ffirst = data.freq_eval_pts.begin() + f_index * data.f_ppr;
  const auto flast = ffirst + data.f_ppr;
  const std::vector<double> local_gamma(ffirst, flast);

  const auto mfirst = data.lag_freq_mults.begin() + f_index * data.f_ppr;
  const auto mlast = mfirst + data.f_ppr;
  const std::vector<double> local_mult(mfirst, mlast);

  const auto efirst = data.erec_coeffs[cell - 1].begin() + f_index * data.f_ppr;
  const auto elast = efirst + data.f_ppr;
  const std::vector<double> local_erec(efirst, elast);

  const auto sfirst = data.sigc_coeffs[cell - 1].begin() + f_index * data.f_ppr;
  const auto slast = sfirst + data.f_ppr;
  const std::vector<double> local_sigc(sfirst, slast);

  // return a pair with the interpolated sigma_c and EREC values:
  return std::make_pair<double, double>(
      rtt_dsxx::interpolate::lagrange_1d(local_sigc, local_gamma, local_mult, mygamma),
      rtt_dsxx::interpolate::lagrange_1d(local_erec, local_gamma, local_mult, mygamma));
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Pre-interpolate raw data in temperature and scale sigma_c by density /
 *        opacity normalization factor. This method belongs to the Preinterp_Container struct,
 *        and operates directly on members of that struct.
 *
 * \param[in] Tms     Cell-wise electron temperatures
 * \param[in] density Cell-wise densities
 * \param[in] data    Container of pointwise sigma_c / EREC
 */
void Preinterp_Container::preinterp_sigc_erec(const std::vector<double> &Tms,
                                              const std::vector<double> &density,
                                              const Edep_Container &data) {

  rtt_units::PhysicalConstexprs<rtt_units::CGSH> consts;

  // Calculate electron rest-mass in keV:
  // units: me (g) * c^2 (cm^2 / shk^2) * (shk^2 / s^2) * (m^2 / cm^2) * (kg / g) * (keV / J)
  const double e_mocsq =
      consts.me() * consts.c() * consts.c() * (1e8) * (1e8) * (1e-3) * (1e-2) * (1e-2) * 6.2415e15;
  const double norm_factor = 0.5 * consts.pi() * consts.re() * consts.re() * 0.5 * consts.Na();

  // Loop over all cell temperatures provided:
  for (size_t b = 0; b < Tms.size(); b++) {
    double myetemp = Tms[b] / e_mocsq;

    // find local intervals on which to do interpolation
    // check it against the library (which should always have 2+ break points):
    Insist(data.T_break_pts.size() > 1, "Missing CSK temperature break points");
    if (myetemp < data.T_break_pts[0] || myetemp > data.T_break_pts.back()) {
#if DRACO_DIAGNOSTICS & 3
      std::cout << "WARNING: Material temp " << myetemp << " is outside the bounds of the library ["
                << data.T_break_pts[0] << ", " << data.T_break_pts.back() << "]" << std::endl;
      ;
#endif
      if (myetemp < data.T_break_pts[0]) {
        myetemp = data.T_break_pts[0] * (1.0 + fuzz);
      }
      if (myetemp > data.T_break_pts.back()) {
        myetemp = data.T_break_pts.back() * (1.0 - fuzz);
      }
    }
    const auto T_index =
        std::lower_bound(data.T_break_pts.begin(), data.T_break_pts.end(), myetemp) -
        data.T_break_pts.begin() - 1;
    // Get Lagrange multipliers for this temperature breakpoint region
    const auto Tmfirst = data.lag_T_mults.begin() + T_index * data.T_ppr;
    const auto Tmlast = Tmfirst + data.T_ppr;
    const std::vector<double> local_Tmult(Tmfirst, Tmlast);

    // Get eval points for this temperature breakpoint region
    const auto Tfirst = data.T_eval_pts.begin() + T_index * data.T_ppr;
    const auto Tlast = Tfirst + data.T_ppr;
    const std::vector<double> local_etemp(Tfirst, Tlast);

    // get the data for this "stripe" in etemp
    for (size_t a = 0; a < data.nfe; a++) {
      const auto f_index = std::lower_bound(data.freq_break_pts.begin(), data.freq_break_pts.end(),
                                            data.freq_eval_pts[a]) -
                           data.freq_break_pts.begin() - 1;
      const size_t region_index = f_index * (data.ntbp - 1) + T_index;
      // Get a start iterator into the data:
      const auto sstart = (a % data.f_ppr) * data.T_ppr;
      const auto send = sstart + data.T_ppr;
      const std::vector<double> oned_sigc_data(data.sigc[region_index].begin() + sstart,
                                               data.sigc[region_index].begin() + send);
      const std::vector<double> oned_erec_data(data.erec[region_index].begin() + sstart,
                                               data.erec[region_index].begin() + send);

      sigc_coeffs[b][a] =
          norm_factor * density[b] *
          rtt_dsxx::interpolate::lagrange_1d(oned_sigc_data, local_etemp, local_Tmult, myetemp);
      erec_coeffs[b][a] =
          rtt_dsxx::interpolate::lagrange_1d(oned_erec_data, local_etemp, local_Tmult, myetemp);
    }
  }
}
} // namespace compton_edep

//------------------------------------------------------------------------------------------------//
// End of Compton_Edep.cc
//------------------------------------------------------------------------------------------------//
