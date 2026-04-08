//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   min/mnbrak.hh
 * \brief  Bracket minimum of a function.
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef min_mnbrak_hh
#define min_mnbrak_hh

#include "dsxx/DracoMath.hh"
#include "dsxx/dbc.hh"
#include <algorithm>
#include <cmath>
#include <concepts>
#include <functional>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>

// cSpell:ignore mnbrak, GLIMIT, ulim

namespace rtt_min {

// A C++20 concept for a callable f(x) -> something convertible to double
template <typename F>
concept UnaryDoubleFunction = requires(F f, double x) {
  { std::invoke(f, x) } -> std::convertible_to<double>;
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Bracket minimum of a function starting from ax and bx.
 *
 * \tparam Function A function type that is callable with a single argument and convertible to
 *         double.
 *
 * \param[in,out] ax Lower bound of bracket interval.
 * \param[in,out] bx Estimate of minimum location within bracket. This will **not** be a polished
 *                   minimum; that is the job of rtt_min::brent.
 * \param[in,out] cx Upper bound of bracket interval.
 * \param[out] fa Function value at ax.
 * \param[out] fb Function value at bx.
 * \param[out] fc Function value at cx.
 * \param[in] func Function to be minimized.
 *
 * On return, ax, bx, cx should surround a local minimum: `f(bx) < f(ax) and f(bx) < f(cx)`
 *
 * \note This is based on the classic Numerical Recipes mnbrak routine.
 */
template <UnaryDoubleFunction Function>
void mnbrak(double &ax, double &bx, double &cx, double &fa, double &fb, double &fc, Function func) {
  using namespace std;
  using namespace rtt_dsxx;

  constexpr double GOLD = 1.618034; // default magnification
  constexpr double GLIMIT = 100.0;  // max magnification for parabolic step
  constexpr double TINY = 1.0e-20;

  auto eval = [&](double x) -> double { return static_cast<double>(std::invoke(func, x)); };

  fa = eval(ax);
  fb = eval(bx);

  // Ensure we are moving downhill from ax to bx
  if (fb > fa) {
    std::swap(ax, bx);
    std::swap(fa, fb);
  }

  cx = bx + GOLD * (bx - ax);
  fc = eval(cx);

  // Keep expanding while we are still going downhill
  while (fb > fc) {
    // Parabolic extrapolation
    const double r = (bx - ax) * (fb - fc);
    const double q = (bx - cx) * (fb - fa);
    const double denom = 2.0 * std::copysign(std::max(std::abs(q - r), TINY), q - r);
    double u = bx - ((bx - cx) * q - (bx - ax) * r) / denom;
    const double ulim = bx + GLIMIT * (cx - bx);

    double fu = 0.0;

    // Test various cases for u
    if ((bx - u) * (u - cx) > 0.0) { // u is between bx and cx
      fu = eval(u);
      if (fu < fc) { // Minimum lies between bx and cx
        ax = bx;
        fa = fb;
        bx = u;
        fb = fu;
        return;             // {ax, bx, cx, fa, fb, fc};
      } else if (fu > fb) { // Minimum lies between ax and u
        cx = u;
        fc = fu;
        return; // {ax, bx, cx, fa, fb, fc};
      }

      // Parabolic fit was no help; use default magnification
      u = cx + GOLD * (cx - bx);
      fu = eval(u);
    } else if ((cx - u) * (u - ulim) > 0.0) { // u is between cx and its allowed limit
      fu = eval(u);
      if (fu < fc) {
        bx = cx;
        fb = fc;
        cx = u;
        fc = fu;
        u = cx + GOLD * (cx - bx);
        fu = eval(u);
      }
    } else if ((u - ulim) * (ulim - cx) >= 0.0) { // Limit parabolic u to maximum allowed value
      u = ulim;
      fu = eval(u);
    } else { // Reject parabolic u; use default magnification
      u = cx + GOLD * (cx - bx);
      fu = eval(u);
    }

    // Eliminate oldest point and continue
    ax = bx;
    fa = fb;
    bx = cx;
    fb = fc;
    cx = u;
    fc = fu;
  }
  return;
}

} // end namespace rtt_min

#endif // min_mnbrak_hh

//------------------------------------------------------------------------------------------------//
// end of min/mnbrak.hh
//------------------------------------------------------------------------------------------------//
