//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/Pseudo_Line_Base.cc
 * \author Kent G. Budge
 * \date   Tue Apr  5 08:42:25 MDT 2011
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Pseudo_Line_Base.hh"
#include "c4/C4_Functions.hh"
#include "cdi/CDI.hh"
#include "ds++/DracoMath.hh"
#include "ds++/Packing_Utils.hh"
#include "ode/quad.hh"
#include "ode/rkqs.hh"
#include <fstream>

namespace rtt_cdi_analytic {
using namespace std;
using namespace rtt_ode;
using namespace rtt_dsxx;
using namespace rtt_cdi;

//------------------------------------------------------------------------------------------------//
#ifdef _MSC_VER
double expm1(double const &x) { return std::exp(x) - 1.0; }
#endif

//------------------------------------------------------------------------------------------------//
void Pseudo_Line_Base::setup_(double emin, double emax) {
  srand(seed_);

  int const number_of_lines = number_of_lines_;
  int const number_of_edges = number_of_edges_;

  // Get global range of energy

  rtt_c4::global_min(emin);
  rtt_c4::global_max(emax);

  if (number_of_lines > 0) {
    center_.resize(number_of_lines);
    for (int i = 0; i < number_of_lines; ++i) {
      center_[i] = (emax - emin) * static_cast<double>(rand()) / RAND_MAX + emin;
    }

    // Sort line centers
    sort(center_.begin(), center_.end());
  }
  // else fuzz model: Instead of lines, we add a random opacity to each opacity bin to simulate very
  // fine, unresolvable line structure.

  unsigned ne = abs(number_of_edges);
  for (unsigned i = 0; i < ne; ++i) {
    if (number_of_edges > 0) {
      // normal behavior is to place edges randomly
      edge_[i] = (emax - emin) * static_cast<double>(rand()) / RAND_MAX + emin;
    } else {
      // placed edges evenly; this makes it easier to choose a group structure that aligns with
      // edges (as would likely be done with a production calculation using a real opacity with
      // strong bound-free components)
      edge_[i] = (emax - emin) * (i + 1) / (ne + 1) + emin;
    }
    double C;
    if (nu0_ < 0) {
      size_t N = continuum_table_.size();
      if (N > 0) {
        C = continuum_table_[static_cast<unsigned int>(edge_[i] * static_cast<double>(N) / emax)];
      } else {
        C = (*continuum_)(vector<double>(1, edge_[i]));
      }
    } else {
      double const nu = edge_[i] / nu0_;
      C = C_ + Bn_ / cube(Bd_ + nu) + R_ * square(square(nu));
    }
    edge_factor_[i] = edge_ratio_ * C;
  }

  // Sort edges
  sort(edge_.begin(), edge_.end());
}

//------------------------------------------------------------------------------------------------//
Pseudo_Line_Base::Pseudo_Line_Base(std::shared_ptr<Expression const> const &continuum,
                                   int number_of_lines, double line_peak, double line_width,
                                   int number_of_edges, double edge_ratio, double Tref, double Tpow,
                                   double emin, double emax, unsigned seed)
    : continuum_(continuum), continuum_table_(std::vector<double>()), emax_(-1.0),
      nu0_(-1), // as fast flag
      C_(-1.0), Bn_(-1.0), Bd_(-1.0), R_(-1.0), seed_(seed), number_of_lines_(number_of_lines),
      line_peak_(line_peak), line_width_(line_width), number_of_edges_(number_of_edges),
      edge_ratio_(edge_ratio), Tref_(Tref), Tpow_(Tpow), center_(std::vector<double>()),
      edge_(abs(number_of_edges)), edge_factor_(abs(number_of_edges)) {
  Require(continuum != std::shared_ptr<Expression>());
  Require(line_peak >= 0.0);
  Require(line_width >= 0.0);
  Require(edge_ratio >= 0.0);
  Require(emin >= 0.0);
  Require(emax > emin);
  // Require parameter (other than emin and emax) to be same on all processors

  setup_(emin, emax);
}

//------------------------------------------------------------------------------------------------//
Pseudo_Line_Base::Pseudo_Line_Base(const string &cont_file, int number_of_lines, double line_peak,
                                   double line_width, int number_of_edges, double edge_ratio,
                                   double Tref, double Tpow, double emin, double emax,
                                   unsigned seed)
    : continuum_(), continuum_table_(std::vector<double>()), emax_(emax), nu0_(-1), // as fast flag
      C_(-1.0), Bn_(-1.0), Bd_(-1.0), R_(-1.0), seed_(seed), number_of_lines_(number_of_lines),
      line_peak_(line_peak), line_width_(line_width), number_of_edges_(number_of_edges),
      edge_ratio_(edge_ratio), Tref_(Tref), Tpow_(Tpow), center_(std::vector<double>()),
      edge_(abs(number_of_edges)), edge_factor_(abs(number_of_edges)) {
  Require(cont_file.size() > 0);
  Require(line_peak >= 0.0);
  Require(line_width >= 0.0);
  Require(edge_ratio >= 0.0);
  Require(emin >= 0.0);
  Require(emax > emin);
  // Require parameter (other than emin and emax) to be same on all processors

  ifstream in(cont_file.c_str());
  if (!in) {
    throw invalid_argument(("could not open " + cont_file).c_str());
  }
  while (in) {
    double x, y, z; // we will ignore x and z and fit table to emax
    in >> x >> y >> z;
    continuum_table_.push_back(y);
  }

  rtt_c4::global_max(emax_);

  setup_(emin, emax);
}

//------------------------------------------------------------------------------------------------//
Pseudo_Line_Base::Pseudo_Line_Base(double nu0, double C, double Bn, double Bd, double R,
                                   int number_of_lines, double line_peak, double line_width,
                                   int number_of_edges, double edge_ratio, double Tref, double Tpow,
                                   double emin, double emax, unsigned seed)
    : continuum_(), continuum_table_(std::vector<double>()), emax_(emax), nu0_(nu0), C_(C), Bn_(Bn),
      Bd_(Bd), R_(R), seed_(seed), number_of_lines_(number_of_lines), line_peak_(line_peak),
      line_width_(line_width), number_of_edges_(number_of_edges), edge_ratio_(edge_ratio),
      Tref_(Tref), Tpow_(Tpow), center_(std::vector<double>()), edge_(abs(number_of_edges)),
      edge_factor_(abs(number_of_edges)) {
  Require(nu0_ > 0.0);
  Require(C_ >= 0.0);
  Require(Bn_ >= 0.0);
  Require(Bd_ >= 0.0);
  Require(R_ >= 0.0);
  Require(line_peak >= 0.0);
  Require(line_width >= 0.0);
  Require(edge_ratio >= 0.0);
  Require(emin >= 0.0);
  Require(emax > emin);
  // Require parameter (other than emin and emax) to be same on all processors

  setup_(emin, emax);
}

//------------------------------------------------------------------------------------------------//
//! Packing function for Pseudo_Line_Base objects.
vector<char> Pseudo_Line_Base::pack() const {
  throw std::range_error("sorry, pack not implemented for Pseudo_Line_Base");
  // Because we haven't implemented packing functionality for Expression trees yet.

#if 0
// caculate the size in bytes
    unsigned const size =
        3 * sizeof(double) + 3 * sizeof(int) + continuum_->packed_size();

    vector<char> pdata(size);

    // make a packer
    rtt_dsxx::Packer packer;

    // set the packer buffer
    packer.set_buffer(size, &pdata[0]);


    // pack the data
    continuum_->pack(packer);
    packer << seed_;
    packer << number_of_lines_;
    packer << line_peak_;
    packer << line_width_;
    packer << number_of_edges_;
    packer << edge_ratio_;

    // Check the size
    Ensure (packer.get_ptr() == &pdata[0] + size);

    return pdata;
#endif
}

//------------------------------------------------------------------------------------------------//
double Pseudo_Line_Base::monoOpacity(double const x, double const T) const {

  int const number_of_lines = number_of_lines_;
  double const width = line_width_;
  double const peak = line_peak_;

  double Result;
  if (nu0_ < 0) {
    size_t N = continuum_table_.size();
    if (N > 0) {
      Result = continuum_table_[static_cast<unsigned int>(x * static_cast<double>(N) / emax_)];
    } else {
      Result = (*continuum_)(vector<double>(1, x));
    }
  } else {
    double const nu = x / nu0_;
    Result = C_ + Bn_ / cube(Bd_ + nu) + R_ * square(square(nu));
  }

  if (number_of_lines >= 0) {
    for (int i = 0; i < number_of_lines; ++i) {
      double const nu0 = center_[i];
      double const d = (x - nu0) / (width * nu0);
      // Result += peak*exp(-d*d);
      Result += peak / (1 + d * d);
    }
  } else {
    // Fuzz model. We had better be precalculating opacities for consistent behavior.
    Result += peak * static_cast<double>(rand()) / RAND_MAX;
  }

  unsigned const number_of_edges = abs(number_of_edges_);

  for (unsigned i = 0; i < number_of_edges; ++i) {
    double const nu0 = edge_[i];
    if (x >= nu0) {
      Result += edge_factor_[i] * cube(nu0 / x);
    }
  }
  // if the power is ~0, then pow(a,0) == 1.0.
  if (!rtt_dsxx::soft_equiv(Tpow_, 0.0, std::numeric_limits<double>::epsilon())) {
    Result *= pow(T / Tref_, Tpow_);
  }
  return Result;
}

} // end namespace rtt_cdi_analytic

//------------------------------------------------------------------------------------------------//
// end of Pseudo_Line_Base.cc
//------------------------------------------------------------------------------------------------//
