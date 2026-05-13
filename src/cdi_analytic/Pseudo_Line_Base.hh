//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/Pseudo_Line_Base.hh
 * \author Kent G. Budge
 * \note   Copyright (C) 2011-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_analytic_Pseudo_Line_Base_hh
#define rtt_cdi_analytic_Pseudo_Line_Base_hh

#include "parser/Expression.hh"
#include <cstdint>
#include <cstdio>

namespace rtt_cdi_analytic {
using rtt_parser::Expression;

#ifdef _MSC_VER
double expm1(double const &x);
#endif

//------------------------------------------------------------------------------------------------//
/*!
 * \class Pseudo_Line_Base
 * \brief Defines a random line spectrum for the opacity.
 *
 * The opacity function is a continuum on which is superimposed a number of lines of the specified
 * peak and width. The line locations are chosen at random.
 *
 * The mass opacity coefficient is assumed independent of temperature or density, which allows
 * precalculation of the opacity structure, an important time saver.
 */
class Pseudo_Line_Base {
public:
  enum Averaging : uint8_t {
    NONE,      //!< evaluate opacity at band center
    ROSSELAND, //!< form a Rosseland (transparency) mean
    PLANCK,    //!< form a Planck (extinction) mean

    END_AVERAGING //!< sentinel value
  };

private:
  // Coefficients
  std::shared_ptr<Expression const> continuum_; //!< continuum opacity [cm^2/g]
  std::vector<double> continuum_table_;
  double emax_;

  double nu0_;
  double C_;
  double Bn_;
  double Bd_;
  double R_;

  unsigned seed_;
  int number_of_lines_;
  double line_peak_;  //!< peak line opacity [cm^2/g]
  double line_width_; //!< line width as fraction of line frequency.
  int number_of_edges_;
  double edge_ratio_;

  double Tref_; //!< reference temperature for temperature dependence
  double Tpow_; //!< temperature dependence exponent

  std::vector<double> center_;      //!< line centers for this realization
  std::vector<double> edge_;        //!< edges for this realization
  std::vector<double> edge_factor_; //!< opacity at threshold

  void setup_(double emin, double emax);

public:
  Pseudo_Line_Base(std::shared_ptr<Expression const> const &cont, int number_of_lines,
                   double line_peak, double line_width, int number_of_edges, double edge_ratio,
                   double Tref, double Tpow, double emin, double emax, unsigned seed);

#if 0
  /*!
   * \brief Second constructor for Pseudo_Line_Base.
   * \bug No unit test, but needed by Capsaicin thermal_data/microphysics_parser.cc ->
   *      cdi_analytic/Pseudo_Line_Analytic_Odfmg_Opacity.cc -> Pseudo_Line_Base
   */
  Pseudo_Line_Base(std::string const &cont_file, int number_of_lines, double line_peak,
                   double line_width, int number_of_edges, double edge_ratio, double Tref,
                   double Tpow, double emin, double emax, unsigned seed);

  /*!
   * \brief Third constructor for Pseudo_Line_Base.
   * \bug No unit test, but needed by Capsaicin thermal_data/microphysics_parser.cc ->
   *      cdi_analytic/Pseudo_Line_Analytic_Odfmg_Opacity.cc -> Pseudo_Line_Base
   */
  Pseudo_Line_Base(double nu0, double C, double Bn, double Bd, double R, int number_of_lines,
                   double line_peak, double line_width, int number_of_edges, double edge_ratio,
                   double Tref, double Tpow, double emin, double emax, unsigned seed);
#endif

  //! Constructor for packed state.
  explicit Pseudo_Line_Base(std::vector<char> const &packed);

  virtual ~Pseudo_Line_Base() = default;
  Pseudo_Line_Base(Pseudo_Line_Base const &rhs) = delete;
  Pseudo_Line_Base(Pseudo_Line_Base &&rhs) noexcept = delete;
  Pseudo_Line_Base &operator=(Pseudo_Line_Base const &rhs) = delete;
  Pseudo_Line_Base &operator=(Pseudo_Line_Base &&rhs) noexcept = delete;

  [[nodiscard]] double line_width() const { return line_width_; }
  [[nodiscard]] double Tref() const { return Tref_; }
  [[nodiscard]] double Tpow() const { return Tpow_; }

  //! Pack up the class for persistence.
  [[nodiscard]] std::vector<char> pack() const;

  //! Compute a monochromatic opacity
  [[nodiscard]] double monoOpacity(double nu, double T) const;

  static double BB(double const T, double const x) {
    double const e = expm1(x / T);
    return (e > 0.0) ? x * x * x / e : x * x * T;
  }

  static double DBB(double const T, double const x) {
    double const e = expm1(x / T);
    if (e > 0) {
      double const de = -exp(x / T) * x / (T * T);
      return -x * x * x * de / (e * e);
    } else {
      return x * x;
    }
  }
};

} // end namespace rtt_cdi_analytic

#endif // rtt_cdi_analytic_Pseudo_Line_Base_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_analytic/Pseudo_Line_Base.hh
//------------------------------------------------------------------------------------------------//
