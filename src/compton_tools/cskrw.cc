//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   compton_tools/cskrw.cc
 * \author Andrew Till
 * \date   11 May 2020
 * \brief  Converter of ASCII to binary CSK Compton files. Intended for internal use.
 * \note   Copyright (C) 2020 Triad National Security, LLC. All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "cdi/CDI.hh"
#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/XGetopt.hh"
#include "units/PhysicalConstants.hh"
#include "units/UnitSystem.hh"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using rtt_dsxx::soft_equiv;
using std::cout;
using std::endl;
using std::string;
using UINT64 = uint64_t;
using FP = double;

//------------------------------------------------------------------------------------------------//
/*!
 * \class Sparse_Compton_Data
 *
 * Passed around internally in Dense_Compton_Data; limited metadata
 */
//------------------------------------------------------------------------------------------------//
struct Sparse_Compton_Data {

  //! Construct with zeros and known sizes
  Sparse_Compton_Data() = default;

  //! First group-to with nonzero value. index with T and gfrom; applies to all points
  std::vector<UINT64> first_groups;

  //! cumulative sum of row offsets into data and derivatives. index with T and gfrom; applies to
  //! all points
  std::vector<UINT64> indexes;

  //! sparse version of Dense_Compton_Data's data and derivatives, [point, T, gfrom, gto]
  std::vector<FP> data;
  std::vector<FP> derivatives;
};

//------------------------------------------------------------------------------------------------//
//! \class Dense_Compton_Data
//------------------------------------------------------------------------------------------------//
struct Dense_Compton_Data {
  UINT64 numEvals;
  UINT64 numTs;
  UINT64 numGroups;
  UINT64 numLegMoments;
  std::vector<FP> groupBdrs;
  std::vector<FP> Ts;
  // [eval, moment, T, gfrom, gto]
  std::vector<FP> data;
  // [eval, moment, T, gfrom, gto]
  std::vector<FP> derivatives;

  void resize(UINT64 numfiles, std::string filename);
  void read_from_file(UINT64 eval, std::string filename, bool isnonlin);
  void compute_nonlinear_difference();
  void compute_temperature_derivatives();
  void write_sparse_binary(std::string fileout);
  void print_contents(int verbosity, int precision);

private:
  Sparse_Compton_Data copy_to_sparse();
  void print_sparse(const Sparse_Compton_Data &sd);
  void write_binary(std::string fileout, Sparse_Compton_Data &sd);
};

//------------------------------------------------------------------------------------------------//
/*
 * \brief resize data and set sizes
 * \param[in] numfiles is number of Compton files
 * \param[in] filename is one such file
 */
void Dense_Compton_Data::resize(UINT64 numfiles, std::string filename) {
  // Reserve space for nl difference
  UINT64 numderived = (numfiles >= 4) ? 1 : 0;
  numEvals = numfiles + numderived;

  // Read first line of filename to get sizes
  std::ifstream f(filename);
  Insist(f.is_open(), "Unable to open " + filename);

  // Line 1: sizes
  UINT64 numTbreakpoints = 0; // not used
  numTs = 0;
  numGroups = 0;
  numLegMoments = 0;
  f >> numTbreakpoints >> numTs >> numGroups >> numLegMoments;

  // Set vector lengths
  groupBdrs.resize(numGroups + 1, 0.0);
  Ts.resize(numTs, 0.0);
  UINT64 sz = numEvals * numLegMoments * numTs * numGroups * numGroups;
  data.resize(sz, 0.0);
  derivatives.resize(sz, 0.0);

  f.close();
}

//------------------------------------------------------------------------------------------------//
// read the entire contents of one file
void Dense_Compton_Data::read_from_file(UINT64 eval, std::string filename, bool isnonlin) {
  Insist(eval < numEvals, "eval must be < numEvals");
  std::ifstream f(filename);
  Insist(f.is_open(), "Unable to open " + filename);

  // Normalization / unit change
  const FP mtocm = 100.0;
  const FP classical_electron_radius = mtocm * rtt_units::classicalElectronRadiusSI; // cm

  // Normalization constants for raw CSK data:
  // CSK to cross section: 2 * pi * classicalElectronRadius^2 / 4
  // cross section to opacity: Zbar_over_A * avogadrosNumber
  // opacity to micro xs: 1/avogadrosNumber
  //
  // convert from CSK data to opacity (cm^2/mole)
  const FP csk_opac_norm = 0.25 * 2 * rtt_units::PI * classical_electron_radius *
                           classical_electron_radius * rtt_units::AVOGADRO;
  Ensure(soft_equiv(csk_opac_norm, 2 * 0.037558, 1e-4));
  //
  // convert from opacity (cm^2/mole) to micro xs (cm^2)
  //const FP csk_xs_conv = 1. / rtt_units::AVOGADRO;
  //Ensure(soft_equiv(2 * 0.2003102 * csk_xs_conv, 6.6524587e-25, 1e-4));

  // VALUES USED IN CSK'S PHYSICAL_CONSTANTS.HH
  // Normalization for energy and nonlinear terms
  const FP mec2 = 510.998;
  const FP hplanck = 8.6173303e-8;
  const FP cspeed = 299.79245800;

  // set effective renorm base multipliers
  const FP basescale = csk_opac_norm; // cm^2/mole
  const FP nlbase =
      (4.0 / 9.0) * 2.0 / (hplanck * hplanck * hplanck * cspeed * cspeed) * (mec2 * mec2 * mec2);

  // Line 1: sizes
  UINT64 numTbreakpoints = 0;
  UINT64 numTs_check = 0;
  UINT64 numGroups_check = 0;
  UINT64 numLegMoments_check = 0;
  f >> numTbreakpoints >> numTs_check >> numGroups_check >> numLegMoments_check;
  Check(numTs == numTs_check);
  Check(numGroups == numGroups_check);
  Check(numLegMoments == numLegMoments_check);

  // Line 2: T breakpoints (unused)
  FP throwaway;
  for (UINT64 i = 0; i < numTbreakpoints; ++i) {
    f >> throwaway;
  }

  // Line 3: Group bounds
  for (UINT64 g = 0; g < numGroups + 1U; ++g) {
    f >> groupBdrs[g];
    // scale to keV
    groupBdrs[g] *= mec2;
  }

  // Remaining lines: Temperatures and MG data
  // Format:
  // T
  // gto gfrom moment0 [moment1 moment2 ...]
  // <blankline>
  for (UINT64 iT = 0; iT < numTs; ++iT) {
    f >> Ts[iT];
    // scale to keV
    Ts[iT] *= mec2;
    const FP T4 = Ts[iT] * Ts[iT] * Ts[iT] * Ts[iT];

    const FP linscale = isnonlin ? nlbase * T4 : 1.0;
    const FP renorm = basescale * linscale;

    UINT64 gfrom = 0U;
    UINT64 gto = 0U;
    bool finished = false;
    while (!finished) {
      // Read one line
      f >> gfrom >> gto;
      // 1-based to 0-based
      gfrom -= 1U;
      gto -= 1U;

      // Read xs
      for (UINT64 iL = 0; iL < numLegMoments; ++iL) {
        // Read value
        FP val;
        f >> val;
        val *= renorm;

        // Put in 1D data vector
        const UINT64 loc =
            gto + numGroups * (gfrom + numGroups * (iT + numTs * (iL + numLegMoments * eval)));
        Check(loc < data.size());
        data[loc] = val;
      }

      // Get EOL char
      char c;
      f.get(c);
      Check(c == '\n');

      // Look-ahead to next line to see if a blank line, signalling end of temperature data
      f.get(c);
      if (c == '\n') { // [kt] || c == EOF) { // EOF is not a char. This is alwasy false.
        finished = true;
      } else {
        f.putback(c);
      }
    }
  }
  f.close();
}

//------------------------------------------------------------------------------------------------//
// use 4 evaluations to determine nonlinear difference
// implicit at low E/T; explicit at high E/T
void Dense_Compton_Data::compute_nonlinear_difference() {
  Require(numEvals == 5);

  // Physically, ON - IN.T == IL.T B^-1 - B^-1 OL.
  // Due to quadrature choice and finite-precision arithmetic,
  // the above equality does NOT discretely hold.
  // The LHS suffers from small differences of large numbers at low E/T.
  // The RHS suffers from the same problem when E/T is large and so B is small.
  // To mitigate numerical error and avoid dividing by a B of 0.,
  // we use the RHS at low E/T and use the LHS at high E/T.

  // Caveat Emptor: fN, ON, and IN are scaled to a bg that sums to unity! If a downstream data
  // consumer uses a phi or bg that sums to a T_r^4, then fN, ON, and/or IN should be rescaled by
  // 1/(a * T_e^4)

  // Determine the cutoff between low and high energies (Ecutoff = N * T)
  // Based on the wgt fxn CSK uses in its MG avg, N <= 25.0
  // Using 25.0 is better for equilibrium stimulated scat
  // Using N>=9.0 is better for non-equil stimulated scat
  const FP Ncutoff = 9.0;

  // Planck MG integral
  std::vector<FP> bg(numGroups, 0.0);

  // indexes for the evaluations
  // I for inscattering, O for outscattering, f for difference
  // L for linear, N for nonlinear
  const UINT64 e_IL = 0;
  const UINT64 e_OL = 1;
  const UINT64 e_IN = 2;
  const UINT64 e_ON = 3;
  const UINT64 e_fN = 4;
  // 0th Legendre moment

  for (UINT64 iT = 0; iT < numTs; ++iT) {
    const FP T = Ts[iT];
    const FP Ecutoff = Ncutoff * T;

    // Compute bg[T]
    FP bgsum = 0.0;
    for (UINT64 g = 0; g < numGroups; ++g) {
      const FP Elow = groupBdrs[g];
      const FP Ehigh = groupBdrs[g + 1];
      bg[g] = rtt_cdi::integratePlanckSpectrum(Elow, Ehigh, T);
      bgsum += bg[g];
    }
    bgsum = bgsum > 0.0 ? bgsum : 1.0;
    //Normalize bg[T] (needed when T is near first or last group bounds)
    for (UINT64 g = 0; g < numGroups; ++g) {
      bg[g] /= bgsum;
    }
    // First pass on nldiff
    FP sumlin = 0.0;
    FP sumnonlin = 0.0;
    for (UINT64 iL = 0; iL < numLegMoments; ++iL) {
      for (UINT64 gfrom = 0; gfrom < numGroups; ++gfrom) {
        for (UINT64 gto = 0; gto < numGroups; ++gto) {
          // Look at left side of group bounds (use less implicit)
          const FP Eto = groupBdrs[gto];
          const FP Efrom = groupBdrs[gfrom];
          const bool lowE = Eto <= Ecutoff && Efrom <= Ecutoff;

          // Planck is equilibrium distribution
          const FP bgto = bg[gto];
          const FP bgfrom = bg[gfrom];

          std::array<FP, 4> vals;
          // use scattering matrix (no transpose) for outscattering
          for (UINT64 eval : {e_OL, e_ON}) {
            const UINT64 loc =
                gto + numGroups * (gfrom + numGroups * (iT + numTs * (iL + numLegMoments * eval)));
            vals[eval] = data[loc];
          }
          // use transpose of scattering matrix for inscattering
          for (UINT64 eval : {e_IL, e_IN}) {
            const UINT64 loc =
                gfrom + numGroups * (gto + numGroups * (iT + numTs * (iL + numLegMoments * eval)));
            vals[eval] = data[loc];
          }

          // Avoid dividing by zero
          constexpr FP eps = 100 * std::numeric_limits<FP>::min();
          const bool bzero = bgto <= eps || bgfrom <= eps;

          // Take differences of spontaneous and induced rates at equilibrium
          const FP impldiff = (bzero) ? 0.0 : vals[e_IL] / bgfrom - vals[e_OL] / bgto;
          const FP expldiff = vals[e_ON] - vals[e_IN];

          // For low E/T, store impldiff; for high E/T, store expldiff
          const UINT64 loc_fN =
              gto + numGroups * (gfrom + numGroups * (iT + numTs * (iL + numLegMoments * e_fN)));
          data[loc_fN] = (lowE) ? impldiff : expldiff;

          // Keep sums of 0th moment rates for later ratio
          if (iL == 0U) {
            sumlin += bgto * impldiff * bgfrom;
            sumnonlin += bgto * data[loc_fN] * bgfrom;
          }
        }
      }
    }

    // Rescale the nonlin diff to get exact detailed balance at equilibrium
    const FP scalenl = sumlin / sumnonlin;
    // we hope scalenl is within a percent or less of 1
    Check(scalenl < 1.2 && scalenl > 0.8);
    for (UINT64 iL = 0; iL < numLegMoments; ++iL) {
      for (UINT64 gfrom = 0; gfrom < numGroups; ++gfrom) {
        for (UINT64 gto = 0; gto < numGroups; ++gto) {
          const UINT64 loc_fN =
              gto + numGroups * (gfrom + numGroups * (iT + numTs * (iL + numLegMoments * e_fN)));
          data[loc_fN] *= scalenl;
        }
      }
    }
  }
}

//------------------------------------------------------------------------------------------------//
// Compute and limit temperature derivatives of data
void Dense_Compton_Data::compute_temperature_derivatives() {
  // Check that temperature grid is valid (part 1/2)
  if (numTs < 2) {
    std::fill(derivatives.begin(), derivatives.end(), 0.0);
    std::cerr << "WARNING: Cannot construct derivatives with only one "
                 "temperature. Aborting routine.";
    std::cerr << std::endl;
    return;
  }

  // Check that temperature grid is valid (part 2/2)
  bool validTs = true;
  for (UINT64 it = 0; it < (numTs - 1U); ++it) {
    FP T1 = Ts[it];
    FP T2 = Ts[it + 1];
    validTs = validTs && (T1 < T2);
  }
  Insist(validTs, "Temperatures are not monotonically increasing and unique.\n");

  // Temporary array for finite differences
  const UINT64 fd_sz = numGroups * numGroups * (numTs - 1U);
  std::vector<FP> finite_diffs(fd_sz, 0.0);

  for (UINT64 eval = 0; eval < numEvals; ++eval) {
    for (UINT64 iL = 0; iL < numLegMoments; ++iL) {
      std::fill(finite_diffs.begin(), finite_diffs.end(), 0.0);

      // Step 1: Compute finite difference in each temperature interval
      for (UINT64 iT = 0; iT < (numTs - 1U); ++iT) {
        const FP inv_dT = 1.0 / (Ts[iT + 1U] - Ts[iT]);
        for (UINT64 gfrom = 0; gfrom < numGroups; ++gfrom) {
          for (UINT64 gto = 0; gto < numGroups; ++gto) {

            const UINT64 loc_m =
                gto + numGroups * (gfrom + numGroups * (iT + numTs * (iL + numLegMoments * eval)));
            const UINT64 loc_p =
                gto +
                numGroups * (gfrom + numGroups * ((iT + 1U) + numTs * (iL + numLegMoments * eval)));
            const FP fd = (data[loc_p] - data[loc_m]) * inv_dT;

            const UINT64 loc_fd = gto + numGroups * (gfrom + numGroups * iT);
            finite_diffs[loc_fd] = fd;
          }
        }
      }

      // Step 2: Compute derivatives at each temperature point with a
      // limiter based on the finite differences in surrounding intervals

      // Step 2a: First temperature (no limiter; first-order estimate)
      {
        UINT64 iT = 0;
        const UINT64 offset = numGroups * (iT + numTs * (iL + numLegMoments * eval));
        const UINT64 offset_fd = numGroups * iT;
        for (UINT64 gfrom = 0; gfrom < numGroups; ++gfrom) {
          for (UINT64 gto = 0; gto < numGroups; ++gto) {
            const UINT64 loc = gto + numGroups * (gfrom + offset);
            const UINT64 loc_fd = gto + numGroups * (gfrom + offset_fd);
            derivatives[loc] = finite_diffs[loc_fd];
          }
        }
      }

      // Step 2b: Last temperature (no limiter; first-order estimate)
      {
        UINT64 iT = numTs - 1U;
        const UINT64 offset = numGroups * (iT + numTs * (iL + numLegMoments * eval));
        const UINT64 offset_fd = numGroups * (iT - 1U);
        for (UINT64 gfrom = 0; gfrom < numGroups; ++gfrom) {
          for (UINT64 gto = 0; gto < numGroups; ++gto) {
            const UINT64 loc = gto + numGroups * (gfrom + offset);
            const UINT64 loc_fd = gto + numGroups * (gfrom + offset_fd);
            derivatives[loc] = finite_diffs[loc_fd];
          }
        }
      }

      // Step 2c: Interior temperatures (with limiters)
      // m is left of iT; p is right of iT; p and iT have same indices
      for (UINT64 iT = 1U; iT < (numTs - 1U); ++iT) {
        const FP dT_m = Ts[iT] - Ts[iT - 1U];
        const FP dT_p = Ts[iT + 1U] - Ts[iT];
        const FP f_m = (2. * dT_m + dT_p) / (3. * (dT_m + dT_p));
        const FP f_p = 1.0 - f_m;

        const UINT64 offset = numGroups * (iT + numTs * (iL + numLegMoments * eval));

        const UINT64 offset_fd_m = numGroups * (iT - 1U);
        const UINT64 offset_fd_p = numGroups * iT;

        for (UINT64 gfrom = 0; gfrom < numGroups; ++gfrom) {
          for (UINT64 gto = 0; gto < numGroups; ++gto) {
            const UINT64 loc_fd_m = gto + numGroups * (gfrom + offset_fd_m);
            const UINT64 loc_fd_p = gto + numGroups * (gfrom + offset_fd_p);
            const FP fd_m = finite_diffs[loc_fd_m];
            const FP fd_p = finite_diffs[loc_fd_p];

            // For more information on method, see
            // scipy.interpolate.pchip and
            // http://dx.doi.org/10.1137/1.9780898717952

            const int sign_m = int(0.0 < fd_m) - int(fd_m < 0.0);
            const int sign_p = int(0.0 < fd_p) - int(fd_p < 0.0);

            const FP d = (sign_m * sign_p > 0) ? (fd_m * fd_p) / (f_m * fd_m + f_p * fd_p) : 0.0;
            const UINT64 loc = gto + numGroups * (gfrom + offset);
            derivatives[loc] = d;
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------------------------//
// Sparsify data and print to binary
void Dense_Compton_Data::write_sparse_binary(std::string fileout) {
  Sparse_Compton_Data sd = copy_to_sparse();
  print_sparse(sd);
  write_binary(fileout, sd);
}

//------------------------------------------------------------------------------------------------//
// Sparsify data
Sparse_Compton_Data Dense_Compton_Data::copy_to_sparse() {

  // Size of smallest non-zero value
  //const FP cutoff = 1e5 * std::numeric_limits<FP>::min();
  const FP cutoff = 1e-210;

  // Determine the number of non-zero entries per col (brute-force loop)
  const UINT64 fg_sz = numGroups * numTs;
  std::vector<UINT64> first_groups(fg_sz, numGroups);
  std::vector<UINT64> end_groups(fg_sz, 0U);

  // Ensure diagonal is included
  for (UINT64 iT = 0; iT < numTs; ++iT) {
    for (UINT64 gfrom = 0; gfrom < numGroups; ++gfrom) {
      const UINT64 loc_fg = gfrom + numGroups * iT;
      // first_groups is inclusive, end_groups is exclusive
      first_groups[loc_fg] = gfrom;
      end_groups[loc_fg] = gfrom + 1U;
    }
  }

  for (UINT64 eval = 0; eval < numEvals; ++eval) {
    for (UINT64 iL = 0; iL < numLegMoments; ++iL) {
      for (UINT64 iT = 0; iT < numTs; ++iT) {
        for (UINT64 gfrom = 0; gfrom < numGroups; ++gfrom) {
          for (UINT64 gto = 0; gto < numGroups; ++gto) {
            const UINT64 iT_m = iT > 0 ? iT - 1U : 0;
            const UINT64 iT_p = iT < (numTs - 1U) ? iT + 1U : numTs - 1U;
            const UINT64 loc_m =
                gto +
                numGroups * (gfrom + numGroups * (iT_m + numTs * (iL + numLegMoments * eval)));
            const UINT64 loc =
                gto + numGroups * (gfrom + numGroups * (iT + numTs * (iL + numLegMoments * eval)));
            const UINT64 loc_p =
                gto +
                numGroups * (gfrom + numGroups * (iT_p + numTs * (iL + numLegMoments * eval)));
            const FP val_m = std::fabs(data[loc_m]);
            const FP val = std::fabs(data[loc]);
            const FP val_p = std::fabs(data[loc_p]);
            // If datapoint is nonzero at current or bounding temperatures, include in sparse dataset
            if (val_m > cutoff || val > cutoff || val_p > cutoff) {
              const UINT64 loc_fg = gfrom + numGroups * iT;
              first_groups[loc_fg] = std::min(gto, first_groups[loc_fg]);
              end_groups[loc_fg] = std::max(gto + 1U, end_groups[loc_fg]);
            }
          }
        }
      }
    }
  }

  // Determine sizes and use to compute offsets
  const UINT64 i_sz = numGroups * numTs + 1U;
  std::vector<UINT64> indexes(i_sz, 0U);
  for (UINT64 i = 0; i < fg_sz; ++i) {
    const UINT64 di = end_groups[i] - first_groups[i];
    indexes[i + 1U] = indexes[i] + di;
  }

  // Save all Legendre moments for first (in_lin) eval
  // and save 0th Legendre moment for other (out_lin, nldiff) evals
  const UINT64 numBinaryEvals = (numEvals > 1U) ? 3U : 1U;
  const UINT64 numPoints = (numLegMoments + numBinaryEvals - 1U);
  const UINT64 numPerPoint = indexes[i_sz - 1U];
  const UINT64 numNonZeros = numPerPoint * numPoints;
  std::vector<FP> sparse_data(numNonZeros, 0.0);
  std::vector<FP> sparse_derivatives(numNonZeros, 0.0);

  /*
  std::cout << "DBG numBinaryEvals numPoints numPerPoint numNonZeros "
            << numBinaryEvals << ' ' << numPoints << ' ' << numPerPoint << ' '
            << numNonZeros << '\n';
   */

  const std::array<UINT64, 3> evalsToUse = {0, 1, 4};
  for (UINT64 iuse = 0; iuse < numBinaryEvals; ++iuse) {
    const UINT64 eval = evalsToUse[iuse];
    const UINT64 numLegUse = (eval > 0) ? 1U : numLegMoments;
    for (UINT64 iL = 0; iL < numLegUse; ++iL) {
      const UINT64 point = (eval > 0) ? numLegMoments + iuse - 1U : iL;
      for (UINT64 iT = 0; iT < numTs; ++iT) {
        for (UINT64 gfrom = 0; gfrom < numGroups; ++gfrom) {
          const UINT64 loc_fg = gfrom + numGroups * iT;
          const UINT64 first = first_groups[loc_fg];
          const UINT64 offset = indexes[loc_fg] + point * numPerPoint;
          const UINT64 sz = indexes[loc_fg + 1U] - indexes[loc_fg];
          for (UINT64 dg = 0; dg < sz; ++dg) {
            // d for dense; s for sparse
            const UINT64 gto = dg + first;
            const UINT64 loc_d =
                gto + numGroups * (gfrom + numGroups * (iT + numTs * (iL + numLegMoments * eval)));
            const UINT64 loc_s = dg + offset;
            sparse_data[loc_s] = data[loc_d];
            sparse_derivatives[loc_s] = derivatives[loc_d];
          }
        }
      }
    }
  }

  // Create sparse data structure
  Sparse_Compton_Data sd;
  std::swap(sd.first_groups, first_groups);
  std::swap(sd.indexes, indexes);
  std::swap(sd.data, sparse_data);
  std::swap(sd.derivatives, sparse_derivatives);

  return sd;
}

//------------------------------------------------------------------------------------------------//
// Debug print sparse data
void Dense_Compton_Data::print_sparse(const Sparse_Compton_Data &sd) {
  std::cout << "sparse sizes : " << sd.first_groups.size() << ' ' << sd.indexes.size() << ' '
            << sd.data.size() << ' ' << sd.derivatives.size() << '\n';

  if (false) {
    std::cout << "PRINT CONTENTS (point 0)\n";
    for (UINT64 iT = 0; iT < numTs; ++iT) {
      for (UINT64 gfrom = 0; gfrom < numGroups; ++gfrom) {
        UINT64 loc = gfrom + numGroups * iT;
        UINT64 fg = sd.first_groups[loc];
        UINT64 strt = sd.indexes[loc];
        UINT64 endd = sd.indexes[loc + 1];
        for (UINT64 ii = strt; ii < endd; ++ii) {
          UINT64 gto = (ii - strt) + fg;
          std::cout << iT << ' ' << gfrom << ' ' << gto << ' ' << std::setprecision(2)
                    << sd.data[ii] << ' ' << std::setprecision(2) << sd.derivatives[ii] << '\n';
        }
      }
    }
  }

  UINT64 counter = 0;
  UINT64 line = 8;
  if (false) {
    std::cout << '\n';
    std::cout << "first_groups:\n";
    counter = 0;
    for (const UINT64 fg : sd.first_groups) {
      std::cout << fg << ' ';
      if (++counter % line == 0)
        std::cout << '\n';
    }
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "last_groups (exclusive):\n";
    counter = 0;
    for (UINT64 i = 0; i < sd.first_groups.size(); ++i) {
      UINT64 fg = sd.first_groups[i];
      UINT64 sz = sd.indexes[i + 1] - sd.indexes[i];
      UINT64 lg = fg + sz;
      std::cout << lg << ' ';
      if (++counter % line == 0)
        std::cout << '\n';
    }
    std::cout << '\n';
  }

  if (false) {
    std::cout << '\n';
    std::cout << "indexes:\n";
    counter = 0;
    for (const UINT64 i : sd.indexes) {
      std::cout << i << ' ';
      if (++counter % line == 0)
        std::cout << '\n';
    }
    std::cout << '\n';
  }

  if (false) {
    std::cout << '\n';
    std::cout << "sparse data:\n";
    counter = 0;
    for (const FP d : sd.data) {
      std::cout << std::setprecision(4) << d << ' ';
      if (std::fabs(d) > 1e-210)
        ++counter;
      if (counter % line == 0)
        std::cout << '\n';
    }
    std::cout << '\n';
  }

  if (false) {
    std::cout << '\n';
    std::cout << "sparse derivatives:\n";
    counter = 0;
    for (const FP d : sd.derivatives) {
      std::cout << std::setprecision(4) << d << ' ';
      if (std::fabs(d) > 1e-210)
        ++counter;
      if (counter % line == 0)
        std::cout << '\n';
    }
    std::cout << '\n';
  }
}

//------------------------------------------------------------------------------------------------//
// Write to binary
void Dense_Compton_Data::write_binary(std::string fileout, Sparse_Compton_Data &sd) {
  auto fout = std::ofstream(fileout, std::ios::out | std::ios::binary);

  // binary type
  // (using vector of char to avoid C-style array of char)
  std::vector<char> filetype = {' ', 'c', 's', 'k', ' ', '\0'};
  fout.write(&filetype[0], filetype.size() * sizeof(char));

  UINT64 version_major = 1;
  UINT64 version_minor = 0;
  // ordering: 0 means leg inside; 1 means leg outside
  UINT64 binary_ordering = 1;

  fout.write(reinterpret_cast<char *>(&version_major), sizeof(UINT64));
  fout.write(reinterpret_cast<char *>(&version_minor), sizeof(UINT64));
  fout.write(reinterpret_cast<char *>(&binary_ordering), sizeof(UINT64));

  UINT64 numBinaryEvals = (numEvals > 1U) ? 3U : 1U;

  UINT64 tsz = Ts.size();
  UINT64 egsz = groupBdrs.size();
  UINT64 fgsz = sd.first_groups.size();
  UINT64 isz = sd.indexes.size();
  UINT64 dsz = sd.data.size();

  Check(tsz == numTs);
  Check(egsz == (numGroups + 1U));
  Check(fgsz == (numGroups * numTs));
  Check(isz == (numGroups * numTs + 1U));
  Check(dsz == sd.derivatives.size());

  /*
  std::cout << "DBG binary file v" << version_major << '.' << version_minor
            << " ordering " << binary_ordering << '\n';
   */

  // sizes
  fout.write(reinterpret_cast<char *>(&numTs), sizeof(UINT64));
  fout.write(reinterpret_cast<char *>(&numGroups), sizeof(UINT64));
  fout.write(reinterpret_cast<char *>(&numLegMoments), sizeof(UINT64));
  fout.write(reinterpret_cast<char *>(&numBinaryEvals), sizeof(UINT64));
  fout.write(reinterpret_cast<char *>(&fgsz), sizeof(UINT64));
  fout.write(reinterpret_cast<char *>(&isz), sizeof(UINT64));
  fout.write(reinterpret_cast<char *>(&dsz), sizeof(UINT64));

  // data
  fout.write(reinterpret_cast<char *>(&Ts[0]), tsz * sizeof(FP));
  fout.write(reinterpret_cast<char *>(&groupBdrs[0]), egsz * sizeof(FP));
  fout.write(reinterpret_cast<char *>(&sd.first_groups[0]), fgsz * sizeof(UINT64));
  fout.write(reinterpret_cast<char *>(&sd.indexes[0]), isz * sizeof(UINT64));
  fout.write(reinterpret_cast<char *>(&sd.data[0]), dsz * sizeof(FP));
  fout.write(reinterpret_cast<char *>(&sd.derivatives[0]), dsz * sizeof(FP));

  fout.close();
}

//------------------------------------------------------------------------------------------------//
// print contents of struct
void Dense_Compton_Data::print_contents(int verbosity, int precision) {
  if (verbosity != 0)
    std::cout << "Printing contents at precision " << precision << "...\n";

  // Print sizes
  if (verbosity > 0) {
    std::cout << '\n';
    std::cout << "numEvals " << numEvals << '\n';
    std::cout << "numTs " << numTs << '\n';
    std::cout << "numGroups " << numGroups << '\n';
    std::cout << "numLegMoments " << numLegMoments << '\n';
  }

  // Print grids
  if (verbosity > 1) {
    std::cout << '\n';
    std::cout << "Group boundaries (keV):\n";
    for (UINT64 g = 0; g < numGroups + 1U; ++g) {
      if (g > 0)
        std::cout << ' ';
      std::cout << std::setprecision(precision) << groupBdrs[g];
    }
    std::cout << '\n';

    std::cout << "Temperatures (keV):\n";
    for (UINT64 iT = 0; iT < numTs; ++iT) {
      if (iT > 0)
        std::cout << ' ';
      std::cout << std::setprecision(precision) << Ts[iT];
    }
    std::cout << '\n';
  }

  // Print contents
  if (verbosity > 2) {
    std::cout << '\n';
    std::vector<std::string> evalNames = {"in_lin", "out_lin", "in_nonlin", "out_nonlin", "nldiff"};
    for (UINT64 eval = 0; eval < numEvals; ++eval) {
      std::cout << "Eval: " << evalNames[eval] << '\n';
      for (UINT64 iL = 0; iL < numLegMoments; ++iL) {
        std::cout << "Legendre moment: " << iL << '\n';

        for (UINT64 iT = 0; iT < numTs; ++iT) {
          std::cout << "Temperature (keV): " << Ts[iT] << '\n';

          std::cout << "Data (matrix; cm^2/mole):\n";
          for (UINT64 gto = 0; gto < numGroups; ++gto) {
            if (verbosity <= 3 && gto > 1)
              continue;
            for (UINT64 gfrom = 0; gfrom < numGroups; ++gfrom) {
              if (gfrom > 0)
                std::cout << ' ';
              const UINT64 loc =
                  gto +
                  numGroups * (gfrom + numGroups * (iT + numTs * (iL + numLegMoments * eval)));
              std::cout << std::setprecision(precision) << data[loc];
            }
            std::cout << '\n';
          }

          std::cout << "Derivative in T (matrix; cm^2/mole-keV):\n";
          for (UINT64 gto = 0; gto < numGroups; ++gto) {
            if (verbosity <= 3 && gto > 0)
              continue;
            for (UINT64 gfrom = 0; gfrom < numGroups; ++gfrom) {
              if (gfrom > 0)
                std::cout << ' ';
              const UINT64 loc =
                  gto +
                  numGroups * (gfrom + numGroups * (iT + numTs * (iL + numLegMoments * eval)));
              std::cout << std::setprecision(precision) << derivatives[loc];
            }
            std::cout << '\n';
          }
        }
      }
    }
  }

  if (verbosity != 0) {
    std::cout << '\n';
    std::cout << "...done printing contents\n";
    std::cout << '\n';
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Basic reader of the csk ASCII file format
 */
void read_csk_files(std::string const &basename, int verbosity) {
  // csk data base filename (csk ASCII format required)

#if 1
  std::array<std::string, 2> inouts = {"in", "out"};
  std::array<std::string, 2> lins = {"lin", "nonlin"};
#else
  std::array<std::string, 1> inouts = {"in"};
  std::array<std::string, 1> lins = {"lin"};
  //std::array<std::string, 1> lins = {"nonlin"};
#endif

  // Store Compton data in structure
  Dense_Compton_Data dat;

  // Resize
  {
    UINT64 numfiles = lins.size() * inouts.size();
    std::string const filename = basename + '_' + inouts[0] + '_' + lins[0];
    dat.resize(numfiles, filename);
  }

  // Fill
  UINT64 counter = 0;
  for (std::string lin : lins) {
    for (std::string inout : inouts) {

      ++counter;
      UINT64 eval = counter - 1U;

      std::string const filename = basename + '_' + inout + '_' + lin;
      cout << "Reading file: " << filename << endl;

      bool isnonlin = lin.compare("nonlin") == 0;

      dat.read_from_file(eval, filename, isnonlin);
    }
  }

  // Use filled data to compute derived data
  // (4 evals from raw data plus 1 eval of derived data)
  if (dat.numEvals == 5) {
    dat.compute_nonlinear_difference();
  }
  dat.compute_temperature_derivatives();

  // Save to binary
  std::string fileout = basename + "_b";
  std::cout << "Writing file: " << fileout << '\n';
  dat.write_sparse_binary(fileout);

  // Print
  int precision = 3;
  //precision = 16;
  dat.print_contents(verbosity, precision);

  // Check detailed balance
  if (dat.numEvals == 5) {
    std::cout << "Detailed balance check...\n";
    if (verbosity <= 0)
      std::cout << "T lindiff/nonlindiff\n";

    // Planck MG integral
    std::vector<FP> bg(dat.numGroups, 0.0);

    // indexes for the evaluations
    // I for inscattering, O for outscattering, f for difference
    // L for linear, N for nonlinear
    std::vector<std::string> evalNames = {"in_lin", "out_lin", "in_nonlin", "out_nonlin", "nldiff"};
    const UINT64 e_IL = 0;
    const UINT64 e_OL = 1;
    const UINT64 e_IN = 2;
    const UINT64 e_ON = 3;
    const UINT64 e_fN = 4;
    // 0th Legendre moment
    const UINT64 iL = 0;

    if (verbosity > 0)
      std::cout << '\n';
    for (UINT64 iT = 0; iT < dat.numTs; ++iT) {
      // Get T
      const FP T = dat.Ts[iT];
      if (verbosity > 0)
        std::cout << "Temperature (keV): " << std::setprecision(precision) << T << '\n';

      // Compute bg[T]
      if (verbosity > 1)
        std::cout << "Planck spectrum: ";
      FP bgsum = 0.0;
      for (UINT64 g = 0; g < dat.numGroups; ++g) {
        const FP Elow = dat.groupBdrs[g];
        const FP Ehigh = dat.groupBdrs[g + 1];
        bg[g] = rtt_cdi::integratePlanckSpectrum(Elow, Ehigh, T);
        bgsum += bg[g];

        if (verbosity > 1)
          std::cout << ' ' << std::setprecision(precision) << bg[g];
      }
      if (verbosity > 1)
        std::cout << '\n';
      if (verbosity > 0)
        std::cout << "bgsum (raw): " << std::setprecision(16) << bgsum << '\n';

      // Normalize to bgsum (to match compute_nonlinear_difference)
      for (UINT64 g = 0; g < dat.numGroups; ++g) {
        bg[g] /= bgsum;
      }
      // Compute sums for each eval in equilibrium (I=B)
      std::array<FP, 5> sums = {0.0, 0.0, 0.0, 0.0, 0.0};
      for (UINT64 eval = 0; eval < sums.size(); ++eval) {
        for (UINT64 gfrom = 0; gfrom < dat.numGroups; ++gfrom) {
          FP subsum = 0.0;
          for (UINT64 gto = 0; gto < dat.numGroups; ++gto) {
            // for linear terms, no induced planck[energy_to]
            const FP bgto = (eval >= 2) ? bg[gto] : 1.0;
            const FP bgfrom = bg[gfrom];
            const UINT64 loc =
                gto +
                dat.numGroups *
                    (gfrom + dat.numGroups * (iT + dat.numTs * (iL + dat.numLegMoments * eval)));
            const FP val = dat.data[loc];
            subsum += bgto * val * bgfrom;
          }
          sums[eval] += subsum;
        }
        if (verbosity > 1)
          std::cout << evalNames[eval] << " sum: " << sums[eval] << '\n';
      }

      // Print detailed-balance differences
      const FP scale = 0.75e4 / T;
      const FP lindiff = (sums[e_IL] - sums[e_OL]) * scale;
      const FP nonlindiff_raw = (sums[e_ON] - sums[e_IN]) * scale;
      const FP nonlindiff_use = sums[e_fN] * scale;
      const FP ratio_raw = lindiff / nonlindiff_raw - 1.0;
      const FP ratio_use = lindiff / nonlindiff_use - 1.0;
      if (verbosity > 1)
        std::cout << "lindiff nonlindiff-RAW nonlindiff-USE: " << std::setprecision(6) << lindiff
                  << ' ' << nonlindiff_raw << ' ' << nonlindiff_use << '\n';
      if (verbosity > 0) {
        std::cout << "lindiff / nonlindiff_raw - 1: " << ratio_raw << '\n';
        std::cout << "lindiff / nonlindiff_use - 1: " << ratio_use << '\n';
      } else
        std::cout << T << ' ' << ratio_use << '\n';

      if (verbosity > 0)
        std::cout << '\n';
    }
    std::cout << "...detailed balance check done\n";
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  // Process known command line arguments:
  rtt_dsxx::XGetopt::csmap long_options;
  long_options['h'] = "help";
  long_options['v'] = "version";
  std::map<char, std::string> help_strings;
  help_strings['h'] = "print this message.";
  help_strings['v'] = "print version information and exit.";
  rtt_dsxx::XGetopt program_options(argc, argv, long_options, help_strings);

  std::string const helpstring("\nUsage: cskrw [-hv] "
                               "<csk_base_filename>\n¡Under active development!\n");

  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'v': // --version
      cout << argv[0] << ": version " << rtt_dsxx::release() << endl;
      return 0;

    case 'h': // --help
      cout << argv[0] << ": version " << rtt_dsxx::release() << helpstring << endl;
      return 0;
    }
  }

  // Assume last command line argument is the name of the IPCRESS file.
  std::string const filename = string((argc > 1) ? argv[argc - 1] : "csk");

  try {
    int verbosity = 0;
    verbosity = 1;
    //verbosity = 2;
    //verbosity = 3;
    //verbosity = 4;
    read_csk_files(filename, verbosity);
  } catch (rtt_dsxx::assertion &excpt) {
    cout << "While attempting to read csk file, " << excpt.what() << endl;
    return 1;
  }

  return 0;
}

//------------------------------------------------------------------------------------------------//
// end of cskrw.cc
//------------------------------------------------------------------------------------------------//
