//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/Pseudo_Line_Analytic_MultigroupOpacity.cc
 * \author Kent G. Budge
 * \date   Tue Apr  5 08:42:25 MDT 2011
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Pseudo_Line_Analytic_MultigroupOpacity.hh"
#include "cdi/CDI.hh"
#include "ds++/DracoMath.hh"
#include "ds++/Packing_Utils.hh"
#include "ode/quad.hh"
#include "ode/rkqs.hh"

namespace rtt_cdi_analytic {
using namespace std;
using namespace rtt_ode;
using namespace rtt_dsxx;
using namespace rtt_cdi;

using sf_double = Analytic_MultigroupOpacity::sf_double;
using vf_double = Analytic_MultigroupOpacity::vf_double;

//------------------------------------------------------------------------------------------------//
class PLP_Functor {
public:
  using return_type = double;
  PLP_Functor(Pseudo_Line_Base const *ptr, double const T) : ptr_(ptr), T_(T) {}
  double operator()(double x);

private:
  Pseudo_Line_Base const *ptr_;
  double T_;
};

double PLP_Functor::operator()(double x) {
  return ptr_->monoOpacity(x, T_) * Pseudo_Line_Base::BB(T_, x);
}

//------------------------------------------------------------------------------------------------//
class PLPW_Functor {
public:
  using return_type = double;
  explicit PLPW_Functor(double const T) : T_(T) {}
  double operator()(double x);

private:
  double T_;
};

double PLPW_Functor::operator()(double x) { return Pseudo_Line_Base::BB(T_, x); }

//------------------------------------------------------------------------------------------------//
class PLR_Functor {
public:
  using return_type = double;
  PLR_Functor(Pseudo_Line_Base const *ptr, double const T) : ptr_(ptr), T_(T) {}
  double operator()(double x);

private:
  Pseudo_Line_Base const *ptr_;
  double T_;
};

double PLR_Functor::operator()(double x) {
  return Pseudo_Line_Base::DBB(T_, x) / ptr_->monoOpacity(x, T_);
}

//------------------------------------------------------------------------------------------------//
class PLRW_Functor {
public:
  using return_type = double;
  explicit PLRW_Functor(double const T) : T_(T) {}
  double operator()(double x);

private:
  double T_;
};

double PLRW_Functor::operator()(double x) { return Pseudo_Line_Base::DBB(T_, x); }

//------------------------------------------------------------------------------------------------//
Pseudo_Line_Analytic_MultigroupOpacity::Pseudo_Line_Analytic_MultigroupOpacity(
    sf_double const &group_bounds_in, rtt_cdi::Reaction const reaction_in,
    std::shared_ptr<Expression const> const &continuum, unsigned number_of_lines, double line_peak,
    double line_width, unsigned number_of_edges, double edge_ratio, double Tref, double Tpow,
    double emin, double emax, Averaging const averaging, unsigned const qpoints, unsigned seed_in)
    : Analytic_MultigroupOpacity(group_bounds_in, reaction_in),
      Pseudo_Line_Base(continuum, number_of_lines, line_peak, line_width, number_of_edges,
                       edge_ratio, Tref, Tpow, emin, emax, seed_in),
      averaging_(averaging), qpoints_(qpoints) {}

//------------------------------------------------------------------------------------------------//
//! Packing function
//! \bug No unit tests so commenting out until a test can be added.
Analytic_MultigroupOpacity::sf_char Pseudo_Line_Analytic_MultigroupOpacity::pack() const {
  sf_char pdata = Analytic_MultigroupOpacity::pack();
  Insist(false, "Please add a Draco Unit Test before restoring this code.");
  return pdata;
}

//------------------------------------------------------------------------------------------------//
sf_double Pseudo_Line_Analytic_MultigroupOpacity::getOpacity(double T, double /*rho*/) const {
  sf_double const &group_bounds = this->getGroupBoundaries();
  size_t const number_of_groups = group_bounds.size() - 1;
  sf_double Result(number_of_groups, 0.0);

  double line_width = this->line_width();

  switch (averaging_) {
  case NONE: {
    double g1 = group_bounds[0];
    for (size_t g = 0; g < number_of_groups; ++g) {
      double const g0 = g1;
      g1 = group_bounds[g + 1];
      double const nu = 0.5 * (g0 + g1);
      Result[g] = monoOpacity(nu, T);
    }
  } break;

  case ROSSELAND: {
    PLR_Functor rfunctor(this, T);
    PLRW_Functor wfunctor(T);

    double g1 = group_bounds[0];
    for (size_t g = 0; g < number_of_groups; ++g) {
      double const g0 = g1;
      g1 = group_bounds[g + 1];
      double eps = 1e-5;
      double t = 0.0, b = 0.0;
      if (qpoints_ == 0) {
        double x1 = g0;
        while (x1 < g1) {
          double x0 = x1;
          x1 += 2 * line_width;
          if (x1 > g1)
            x1 = g1;

          t += rtt_ode::quad(rfunctor, x0, x1, eps, rkqs<double, Quad_To_ODE<PLR_Functor>>);
          b += rtt_ode::quad(wfunctor, x0, x1, eps, rkqs<double, Quad_To_ODE<PLRW_Functor>>);
        }
      } else {
        for (unsigned ig = 0; ig < qpoints_; ++ig) {
          double const x = (ig + 0.5) * (g1 - g0) / qpoints_ + g0;

          double w = rtt_cdi::integrateRosselandSpectrum(g0, g1, T);

          t += w / monoOpacity(x, T);
          b += w;
        }
      }

      Result[g] = b / t;
    }
  } break;

  case PLANCK: {
    double g1 = group_bounds[0];
    for (size_t g = 0; g < number_of_groups; ++g) {
      PLP_Functor pfunctor(this, T);
      PLPW_Functor wfunctor(T);

      double const g0 = g1;
      g1 = group_bounds[g + 1];
      double eps = 1e-5;
      double t = 0.0, b = 0.0;
      if (qpoints_ == 0) {
        double x1 = g0;
        while (x1 < g1) {
          double x0 = x1;
          x1 += 2 * line_width;
          if (x1 > g1)
            x1 = g1;

          t += rtt_ode::quad(pfunctor, x0, x1, eps, rkqs<double, Quad_To_ODE<PLP_Functor>>);
          b += rtt_ode::quad(wfunctor, x0, x1, eps, rkqs<double, Quad_To_ODE<PLPW_Functor>>);
        }
      } else {
        for (unsigned ig = 0; ig < qpoints_; ++ig) {
          double const x = (ig + 0.5) * (g1 - g0) / qpoints_ + g0;

          double w = rtt_cdi::integratePlanckSpectrum(g0, g1, T);

          t += w * monoOpacity(x, T);
          b += w;
        }
      }

      Result[g] = t / b;
    }
  } break;

  default:
    Insist(false, "bad case");
  }

  return Result;
}

//------------------------------------------------------------------------------------------------//
vf_double Pseudo_Line_Analytic_MultigroupOpacity::getOpacity(sf_double const &T, double rho) const {
  size_t const n = T.size();
  vf_double Result(n);

  for (size_t i = 0; i < n; ++i) {
    double const Ti = T[i];
    Result[i] = getOpacity(Ti, rho);
  }
  return Result;
}

//------------------------------------------------------------------------------------------------//
vf_double Pseudo_Line_Analytic_MultigroupOpacity::getOpacity(double const T,
                                                             sf_double const &rho) const {
  return {rho.size(), getOpacity(T, rho[0])};
}

//------------------------------------------------------------------------------------------------//
std::string Pseudo_Line_Analytic_MultigroupOpacity::getDataDescriptor() const {
  std::string descriptor;

  rtt_cdi::Reaction const rxn = getReactionType();

  if (rxn == rtt_cdi::TOTAL)
    descriptor = "Pseudo Line Multigroup Total";
  else if (rxn == rtt_cdi::ABSORPTION)
    descriptor = "Pseudo Line Multigroup Absorption";
  else if (rxn == rtt_cdi::SCATTERING)
    descriptor = "Pseudo Line Multigroup Scattering";
  else {
    Insist(0, "Invalid Pseudo Line multigroup model opacity!");
  }

  return descriptor;
}

} // end namespace rtt_cdi_analytic

//------------------------------------------------------------------------------------------------//
// end of Pseudo_Line_Analytic_MultigroupOpacity.cc
//------------------------------------------------------------------------------------------------//
