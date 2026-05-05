//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ode/quad.hh
 * \author Kent Budge
 * \brief  Adaptive quadrature of a function over a specified interval.
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef ode_quad_hh
#define ode_quad_hh

#include "Function_Traits.hh"
#include <vector>

namespace rtt_ode {
//================================================================================================//
/*!
 * \class Quad_To_ODE
 * \brief Helper class for quad
 *
 * This class transforms a quadrature into an ODE problem.  This allows adaptive ODE methods to be
 * used for adaptive quadrature.
 *
 * \arg \a Function A function type for which <code> double operator()(double) </code> is defined.
 */
//================================================================================================//

template <typename Function> class Quad_To_ODE {
public:
  explicit Quad_To_ODE(Function &f) : func(f) {}

  void operator()(double x, std::vector<double> const & /*unused*/, std::vector<double> &dydx) {
    dydx[0] = func(x);
  }

private:
  Function &func;
};

//! Adaptive quadrature of a function over a specified interval.
template <typename Function, typename Rule>
typename Function_Traits<Function>::return_type quad(Function func, double x1, double x2,
                                                     double &eps, Rule rule);

//================================================================================================//
// IMPLEMENTATIONS
//================================================================================================//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Adaptive quadrature of a function over a specified interval.
 *
 * \arg \a Function Function type representing the function to be integrated.
 * \arg \a Rule Function type representing the rule for performing the quadrature. This function
 *         type must be compatible with the signature of rtt_ode::rkqs.
 *
 * \param[in] func Function to be integrated
 * \param[in] x1,x2 Interval over which the function is to be integrated.
 * \param[in,out] eps On entry, contains the desired accuracy.  On return, contains the accuracy
 *              actually achieved.
 * \param[in] rule Ordinary differential equation solver to use to integrate the function, such as
 *              `bsstep` or `rkqs`.
 *
 * \return Definite integral of the function over the specified interval.
 *
 * \pre \c x2>x1
 * \pre \c eps>=0
 *
 * \exception range_error Thrown if the step size becomes too small.
 */
template <typename Function, typename Rule>
typename Function_Traits<Function>::return_type quad(Function func, double x1, double x2,
                                                     double &eps, Rule rule) {
  Require(x2 > x1);
  Require(eps >= 0);

  using namespace std;

  Quad_To_ODE<Function> const derivs(func);

  using Field = typename Function_Traits<Function>::return_type;

  double x = x1; // NOLINT(misc-const-correctness)
  double h = x2 - x1;
  vector<Field> y(1, Field(0.0));
  vector<Field> dydx(1);
  vector<double> yscal(1);
  for (;;) {
    dydx[0] = func(x);
    yscal[0] = 1;
    if ((x + h - x2) * (x + h - x1) > 0)
      h = x2 - x;
    double hdid, hnext;
    rule(y, dydx, x, h, eps, yscal, hdid, hnext, derivs);
    if ((x - x2) * (x2 - x1) >= 0) {
      return y[0];
    }
    if (std::fabs(hnext) <= std::numeric_limits<double>::epsilon() * (x2 - x1)) {
      throw std::range_error("step size too small in quad"); // LCOV_EXCL_LINE
    }
    h = hnext;
  }
}

} // end namespace rtt_ode

#endif // ode_quad_hh

//------------------------------------------------------------------------------------------------//
// end of ode/quad.hh
//------------------------------------------------------------------------------------------------//
