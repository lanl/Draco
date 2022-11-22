//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   compton_tools/Compton_Edep.hh
 * \author Kendra Long
 * \date   10 October 2022
 * \brief  Header file for Compton energy deposition data handling
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_compton_tools_Compton_Edep_hh
#define rtt_compton_tools_Compton_Edep_hh

#include "ds++/Interpolate.hh"
#include <string>

namespace compton_edep {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief This is a simple container to hold pointwise sigma_c and EREC values read from a CSK-
 *        generated file. The data is designed to be interpolated with Lagrange polynomials in
 *        both temperature and frequency. It is most efficient to ``pre-interpolate'' the data in
 *        electron temperature at beginning of cycle; the resulting data can then be interpolated
 *        in frequency ``on the fly'' during the cycle.
 */
class Edep_Container {

public:
  // Convenient sizes, initialized to nonsense
  const size_t ntbp; //! (N)umber of (T)emperature (B)reak(P)oints
  const size_t nte;  //! (N)umber of (T)emperature (E)vals
  const size_t nfbp; //! (N)umber of (F)requency (B)reak(P)oints
  const size_t nfe;  //! (N)umber of (F)requency (E)vals

  // Points per interpolation region:
  const size_t T_ppr; //! Temperature
  const size_t f_ppr; //! Frequency

  // Evaluation points and break points:
  std::vector<double> T_eval_pts;
  std::vector<double> T_break_pts;
  std::vector<double> freq_eval_pts;
  std::vector<double> freq_break_pts;
  std::vector<std::vector<double>> sigc; // this data is laid out BY INTERPOLATION REGION
  std::vector<std::vector<double>> erec; // this data is laid out BY INTERPOLATION REGION

  // This is the Lagrange multiplier data in temperature:
  std::vector<double> lag_T_mults;

  Edep_Container(size_t n_temp_breakpts, size_t n_freq_breakpts, size_t n_temp_evals,
                 size_t n_freq_evals)
      : ntbp(n_temp_breakpts), nte(n_temp_evals), nfbp(n_freq_breakpts), nfe(n_freq_evals),
        T_ppr(nte / (ntbp - 1)), f_ppr(nfe / (nfbp - 1)), T_eval_pts(n_temp_evals, 0.0),
        T_break_pts(n_temp_breakpts, 0.0), freq_eval_pts(n_freq_evals, 0.0),
        freq_break_pts(n_freq_breakpts, 0.0),
        sigc((ntbp - 1) * (nfbp - 1), std::vector<double>(T_ppr * f_ppr, 0.0)),
        erec((ntbp - 1) * (nfbp - 1), std::vector<double>(T_ppr * f_ppr, 0.0)),
        lag_T_mults(n_temp_evals, 0.0) {}

  void compute_multipliers() {
    lag_T_mults = rtt_dsxx::interpolate::lagrange_multipliers(ntbp, T_ppr, T_eval_pts);
  }
};
//------------------------------------------------------------------------------------------------//
/*!
 * \brief This is a simple container to hold coefficients for one-dimensional Lagrange polynomial
 *        interpolation in frequency. It is constructed using the data in an Edep_Container object,
 *        along with cell temperatures and densities. This data is valid for a single cycle.
 */
class Preinterp_Container {

public:
  // Convenient sizes:
  const size_t nfbp; //! (N)umber of (F)requency (B)reak(P)oints
  const size_t nfe;  //! (N)umber of (F)requency (E)vals

  // Points per interpolation region:
  const size_t f_ppr; //! Frequency

  // Evaluation points and break points:
  std::vector<double> freq_eval_pts;
  std::vector<double> freq_break_pts;
  std::vector<std::vector<double>> sigc_coeffs; // this data is laid out BY SPATIAL CELL
  std::vector<std::vector<double>> erec_coeffs; // this data is laid out BY SPATIAL CELL

  // This is the Lagrange multiplier data in frequency:
  std::vector<double> lag_freq_mults;

  Preinterp_Container(const std::vector<double> &cell_temps, const std::vector<double> &density,
                      const Edep_Container &raw_edep)
      : nfbp(raw_edep.nfbp), nfe(raw_edep.nfe), f_ppr(raw_edep.f_ppr),
        freq_eval_pts(raw_edep.freq_eval_pts), freq_break_pts(raw_edep.freq_break_pts),
        sigc_coeffs(cell_temps.size(), std::vector<double>(nfe, 0.0)),
        erec_coeffs(cell_temps.size(), std::vector<double>(nfe, 0.0)) {

    // Compute the multipliers in frequency
    lag_freq_mults = rtt_dsxx::interpolate::lagrange_multipliers(nfbp, f_ppr, freq_eval_pts);
    // Pre-interpolate coefficients in cell temperature
    preinterp_sigc_erec(cell_temps, density, raw_edep);
  }
  //! Use pointwise data + temperatures and densities to pre-compute interpolation coefficients
  //! for sigma_c and EREC.
  void preinterp_sigc_erec(const std::vector<double> &Tms, const std::vector<double> &density,
                           const Edep_Container &data);
};

//! Public interface to read an energy deposition library and return a filled container
Edep_Container read_edep_file(const std::string &infile);
//! Public interface to interpolate to a given frequency in a given cell
std::pair<double, double> interpolate_sigc_erec(const size_t cell, const double freq,
                                                const Preinterp_Container &data);

} // namespace compton_edep

#endif // rtt_compton_tools_Compton_Edep_hh

//------------------------------------------------------------------------------------------------//
// End of Compton_Edep.hh
//------------------------------------------------------------------------------------------------//
