//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/CDI.hh
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:22:06 2000
 * \brief  CDI class header file.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_CDI_hh
#define rtt_cdi_CDI_hh

#include "CPEloss.hh"
#include "EICoupling.hh"
#include "EoS.hh"
#include "GrayOpacity.hh"
#include "MultigroupOpacity.hh"
#include "ds++/Constexpr_Functions.hh"
#include "ds++/FMA.hh"
#include "ds++/Soft_Equivalence.hh"
#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <memory>

//------------------------------------------------------------------------------------------------//
// UNNAMED NAMESPACE
//------------------------------------------------------------------------------------------------//
// Nested unnamed namespace that holds data and services used by the Planckian integration routines.
// The data in this namespace is accessible by the methods in this file only (internal linkage).

namespace {

// Constants used in the Taylor series expansion of the Planckian:
double constexpr coeff_3 = 1.0 / 3.0;
double constexpr coeff_4 = -1.0 / 8.0;
double constexpr coeff_5 = 1.0 / 60.0;
double constexpr coeff_7 = -1.0 / 5040.0;
double constexpr coeff_9 = 1.0 / 272160.0;
double constexpr coeff_11 = -1.0 / 13305600.0;
double constexpr coeff_13 = 1.0 / 622702080.0;
double constexpr coeff_15 = -6.91 / 196151155200.0;
double constexpr coeff_17 = 1.0 / 1270312243200.0;
double constexpr coeff_19 = -3.617 / 202741834014720.0;
double constexpr coeff_21 = 43.867 / 107290978560589824.0;

double constexpr coeff = 0.1539897338202651; // 15/pi^4
double constexpr NORM_FACTOR = 0.25 * coeff; // 15/(4*pi^4);

//------------------------------------------------------------------------------------------------//
/*!
 * \fn inline double taylor_series_planck(double x)
 * \brief Computes the normalized Planck integral via a 21 term Taylor expansion.
 *
 * The taylor expansion of the planckian integral looks as follows:
 *
 * \code
 * I(x) = c0 ( c3 x^3 + c4 x^4 + c5 x^5 + c7 x^7 + c9 x^9 + c11 x^11 + c13 x^13 + c15 x^15
 *      + c17 x^17 + c19 x^19 + c21 x^21 )
 * \endcode
 *
 * If done naively, this requires 136 multiplications. If you accumulate the powers of \f$ x \f$ as
 * you go, it can be done with 24 multiplications.
 *
 * If you express the polynomial as follows:
 *
 * \code
 * I(x) = c0 x^3 ( c3 + x ( c4 + x (c5 + x^2 ( c7 + x^2 ( c9 + x^2 ( c11 + x^2 ( c13 + x^2
 *      * ( c15 + x^2 ( c17 + x^2 ( c19 + x^2 c21 ) ) ) ) ) ) ) ) ) )
 * \endcode
 *
 * the evaluation can be done with 13 multiplications. Furthermore, we do not need to worry about
 * overflow on large powers of \f$ x \f$, since the largest power we compute is \f$ x^3 \f$.
 *
 * \param[in] x The point at which the Planck integral is evaluated.
 * \return The integral value.
 */
static inline double taylor_series_planck(double x) {
  Require(x >= 0.0);

  double const xsqrd = x * x;

  double taylor = FMA(xsqrd, coeff_21, coeff_19);
  taylor = FMA(taylor, xsqrd, coeff_17);
  taylor = FMA(taylor, xsqrd, coeff_15);
  taylor = FMA(taylor, xsqrd, coeff_13);
  taylor = FMA(taylor, xsqrd, coeff_11);
  taylor = FMA(taylor, xsqrd, coeff_9);
  taylor = FMA(taylor, xsqrd, coeff_7);
  taylor = FMA(taylor, xsqrd, coeff_5);
  taylor = FMA(taylor, x, coeff_4);
  taylor = FMA(taylor, x, coeff_3);
  taylor *= x * xsqrd * coeff;

  Ensure(taylor >= 0.0);

  return taylor;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief return the 10-term Polylogarithmic expansion (minus one) for the Planck integral given
 *        \f$ x \f$ and \f$ e^{-x} \f$ (for efficiency)
 */
static double polylog_series_minus_one_planck(double const x, double const eix) {
  Require(x >= 0.0);
  Require(x < 1.0e154); // value will be squared, make sure it's less than sqrt of max double
  Require(rtt_dsxx::soft_equiv(std::exp(-x), eix));

  double const xsqrd = x * x;

  constexpr std::array<double, 9> i_plus_two_inv = {
      0.5000000000000000, // 1/2
      0.3333333333333333, // 1/3
      0.2500000000000000, // 1/4
      0.2000000000000000, // 1/5
      0.1666666666666667, // 1/6
      0.1428571428571429, // 1/7
      0.1250000000000000, // 1/8
      0.1111111111111111, // 1/9
      0.1000000000000000  // 1/10
  };
  double const *curr_inv = i_plus_two_inv.data();

  // initialize to what would have been the calculation of the i=1 term. This saves a number of "mul
  // by one" ops.
  double eixp = eix;
  double li1 = eix;
  double li2 = eix;
  double li3 = eix;
  double li4 = eix;

  // calculate terms 2..10.  This loop has been unrolled by a factor of 3
  for (size_t i = 2; i < 11; i += 3) {
    double const ip0_inv = *curr_inv++;
    eixp *= eix;
    double eixr_ip0 = eixp * ip0_inv;

    double const ip1_inv = *curr_inv++;
    eixp *= eix;
    double eixr_ip1 = eixp * ip1_inv;

    double const ip2_inv = *curr_inv++;
    eixp *= eix;
    double eixr_ip2 = eixp * ip2_inv;

    double const r10 = eixr_ip0;
    double const r11 = eixr_ip1;
    double const r12 = eixr_ip2;

    double const r20 = (eixr_ip0 *= ip0_inv);
    double const r21 = (eixr_ip1 *= ip1_inv);
    double const r22 = (eixr_ip2 *= ip2_inv);

    li1 += r10 + r11 + r12;

    double const r30 = (eixr_ip0 *= ip0_inv);
    double const r31 = (eixr_ip1 *= ip1_inv);
    double const r32 = (eixr_ip2 *= ip2_inv);

    li2 += r20 + r21 + r22;

    double const r40 = (eixr_ip0 * ip0_inv);
    double const r41 = (eixr_ip1 * ip1_inv);
    double const r42 = (eixr_ip2 * ip2_inv);

    li3 += r30 + r31 + r32;
    li4 += r40 + r41 + r42;
  }

  // calculate the lower polylogarithmic integral
  double const poly = -coeff * (xsqrd * x * li1 + 3 * xsqrd * li2 + 6.0 * (x * li3 + li4));

  Ensure(poly <= 0.0);
  return poly;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute the difference between an integrated Planck and Rosseland curves over \f$ (0,\nu)
 * \f$.
 *
 * This helper function is used by CDI::integrate_planck_rosseland (also in this file).
 *
 * ==> The underlying function x^4/(exp(x) - 1) can be difficult to evaluate with double precision
 *     when x is very small. Instead, when x < 1.e-5, we use the first 2 terms in the expansion
 *     x^4/(exp(x) - 1) ~ x^3(1-x/2), remaining terms are x^5/12 + O(x^7).
 *
 * ==> When x is large, e^x reaches floating point overflow rapidly, thus the function is
 *     modified to the form exp(-x)*x^4/(1-exp(-x)). Accuracy of 1-exp(-x) suffers when x is large.
 *     However, std::expm1 should improve the accuracy of that evaluation.
 *
 * ==> The large x fix might be able to be changed in the future if bug listed in man page is
 *     corrected.
 *
 * \param[in] freq The frequency for the upper limit of the integrand.
 * \param[in] exp_freq exp(-freq)
 * \return The difference between the integrated Planck and Rosseland curves over \f$ (0,\nu) \f$.
 */
static double Planck2Rosseland(double const freq, double const exp_freq) {
  Require(freq >= 0.0);
  Require(rtt_dsxx::soft_equiv(exp_freq, std::exp(-freq)));

  // Case 1: if nu/T is sufficiently large, then the evaluation is 0.0.
  //         this evaluation also prevents overflow when evaluating (nu/T)^4.
  if (freq > std::pow(std::numeric_limits<decltype(freq)>::max(), 1.0 / 4.0))
    return 0.0;

  double const freq_3 = freq * freq * freq;

  // Case 2: if nu/T < 1.0e-5, evaluate via Taylor expansion.
  if (freq < 1.0e-5)
    return NORM_FACTOR * freq_3 * (1.0 - 0.5 * freq);

  // Case 3: All other cases
  return NORM_FACTOR * exp_freq * freq_3 * freq / -std::expm1(-freq);
}

} // end of unnamed namespace

namespace rtt_cdi {

//================================================================================================//
/*!
 * \class CDI
 *
 * \brief This class provides a Common Data Interface (CDI) to Atomic, Nuclear and Equation of State
 *        (EOS) data.
 *
 * \sa The Draco packages cdi_gandolf (opacity data) and cdi_eospac (equation of state data).
 *
 * The client must first instantiate concrete Opacity, Nuclear and EOS classes that are derived from
 * the abstract classes found in the CDI package.  A CDI object is then created using these concrete
 * classes as constructor parameters.  Each CDI object will provide access to data for **one**
 * material.  This material may be a mixture (e.g. water) if that mixture has been defined in the
 * underlying data tables.  However, CDI will not mix data table entries to create a new material.
 * This type of mixing should be done by a separate package or the client code.
 *
 * Since this header file does not include the definitions for GrayOpacity or MultigroupOpacity, the
 * calling routine must include these header files.  If the calling routine does not make use of one
 * of these classes then it's definition file does not need to be included, however this will result
 * in the compile-time warning:
 *
 * <pre>
 *       line 69: warning: delete of pointer to incomplete class
 *          delete p;
 * </pre>
 *
 * The user should not worry about this warning as long he/she is not trying to instantiate the
 * class specified by the error message.
 *
 * Including CDI.hh in the client file will automatically include all of the pertinent CDI headers
 * and definitions (rtt_cdi::GrayOpacity, rtt_cdi::MultigroupOpacity, EoS, rtt_cdi::Model,
 * rtt_cdi::Reaction).
 *
 * CDI also contains static services that allow a client to integrate the (normalized) Planckian and
 * Rosseland (\f$ \partial B/\partial T \f$) integrals.  The overloaded functions that perform these
 * services are:
 *
 * \arg integratePlanckSpectrum(double const lowFreq, double const highFreq, double const T)
 *
 * \arg integratePlanckSpectrum(double const freq, double const T)
 *
 * \arg integrateRosselandSpectrum(double const lowf, double const hif, double T);
 *
 * \arg integrate_Rosseland_Planckian_Spectrum(double const lowf, const double hif, double const T,
 *      double& PL, double& ROSL);
 *
 * \arg integratePlanckSpectrum(int const groupIndex, double const T);
 *
 * \arg integratePlanckSpectrum(double const T);
 *
 * \arg integrateRosselandSpectrum(int const groupIndex,  double const T);
 *
 * \arg integrate_Rosseland_Planckian_Spectrum(int const groupIndex, const double T, double& PL,
 *      double& ROSL);
 *
 * The first four forms can be called (CDI::integratePlanckSpectrum()) anytime.  They simply
 * integrate the normalized Planckian or Rosseland over a frequency range.  The next four forms may
 * only be called after the multigroup frequency boundaries have been set.  These boundaries are set
 * after a call, from any CDI object, to setMultigroupOpacity().  The frequency boundaries are
 * stored statically.  After an initial call by any CDI object of setMultigroupOpacity(), the
 * frequency bounds are checked to make sure they do not change.  Changing the boundaries throws an
 * exception. Thus, clients are allowed to view the group structure through any multigroup data set
 * (cdi.mg()->getGroupBoundaries()) or "globally" by calling CDI::getFrequencyGroupBoundaries().
 * The context of usage dictates which type of call to make; the result is invariant.  See the test
 * (tCDI.cc) for usage examples of the CDI Planckian and Rosseland integration routines.
 *
 * We detail the two types of integration here, instead of in the individual methods:
 *
 * The Planckian functions integrate the normalized Planckian that is defined:
 *
 * \f[
 *    b(x) = \frac{15}{\pi^4} \frac{x^3}{e^x - 1}
 * \f]
 *
 * where
 *
 * \f[
 *    x = \frac{h\nu}{kT}
 * \f]
 *
 * and
 *
 * \f[
 *    B(\nu,T)d\nu = \frac{acT^4}{4\pi} b(x)dx
 * \f]
 *
 * where \f$B(\nu,T)\f$ is the Planckian and is defined
 *
 * \f[
 *    B(\nu,T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{h\nu/kt} - 1}
 * \f]
 *
 * The normalized Planckian, integrated from 0 to \f$\infty\f$, equals one. However, depending upon
 * the maximum and minimum group boundaries, the normalized Planck function may integrate to
 * something less than one.
 *
 * This function performs the following integration:
 *\f[
 *      \int_{x_{low}}^{x_{high}} b(x) dx
 *\f]
 * where \f$x_{low}\f$ is calculated from the input low frequency bound and \f$x_high\f$ is
 * calculated from the input high frequency bound.  This integration uses the method of B. Clark
 * (JCP (70)/2, 1987).  We use a 10-term Polylogarithmic expansion for the normalized Planckian,
 * except in the low-x limit, where we use a 21-term Taylor series expansion.
 *
 * The user is responsible for applying the appropriate constants to the result of this integration.
 * For example, to make the result of this integration equivalent to
 * \f[
 *      \int_{\nu_{low}}^{\nu_{high}} B(\nu,T) d\nu
 * \f]
 * then you must multiply by a factor of \f$\frac{acT^4}{4\pi}\f$ where a is the radiation constant.
 * If you want to evaluate expressions like the following:
 *\f[
 *      \int_{4\pi} \int_{\nu_{low}}^{\nu_{high}} B(\nu,T) d\nu d\Omega
 *\f]
 * then you must multiply by \f$acT^4\f$.
 *
 * In the limit of \f$T \rightarrow 0, b(T) \rightarrow 0\f$, therefore we return a hard zero for a
 * temperature equal to a hard zero.
 *
 * The integral is calculated using a polylogarithmic series approximation, except for low
 * frequencies, where the Taylor series approximation is more accurate.  Each of the Taylor and
 * polylogarithmic approximation has a positive truncation error, so they intersect above the
 * correct solution; therefore, we always use the smaller one for a continuous concatenated
 * function.  When both frequency bounds reside above the Planckian peak (above 2.822 T), we skip
 * the Taylor series calculations and use the polylogarithmic series minus one (the minus one is for
 * roundoff control).
 *
 * This Rosseland functions integrate the normalized Rosseland that is defined:
 * \f[
 *    r(x) = \frac{15}{4\pi^4} \frac{x^4 e^x}{(e^x - 1)^2}
 * \f]
 * where
 * \f[
 *    x = \frac{h\nu}{kT}
 * \f]
 * and
 * \f[
 *    R(\nu, T)d\nu = \frac{4 acT^3}{4\pi}r(x)dx
 * \f]
 * where \f$R(\nu, T)\f$ is the Rosseland and is defined
 * \f[
 *    R(\nu, T) = \frac{\partial B(\nu, T)}{\partial T}
 * \f]
 * \f[
 *    B(\nu, T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{\frac{h\nu}{kT}} - 1}
 * \f]
 * The normalized Rosseland, integrated from 0 to \f$\infty\f$, equals one. However, depending upon
 * the maximum and minimum group boundaries, the normalized Rosseland function may integrate to
 * something less than one.
 *
 * This function performs the following integration:
 * \f[
 *      \int_{x_1}^{x_N} r(x) dx
 * \f]
 * where \f$x_1\f$ is the low frequency bound and \f$x_N\f$ is the high frequency bound of the
 * multigroup data set.  This integration uses the method of B. Clark (JCP (70)/2, 1987).  We use a
 * 10-term Polylogarithmic expansion for the normalized Planckian, except in the low-x limit, where
 * we use a 21-term Taylor series expansion.
 *
 * For the Rosseland we can relate the group interval integration to the Planckian group interval
 * integration, by equation 27 in B. Clark paper.
 * \f[
 *     \int_{x_1}^{x_N} r(x) dx = \int_{x_1}^{x_N} b(x) dx
 *     - \frac{15}{4\pi^4} \frac{x^4}{e^x - 1}
 * \f]
 * Therefore our Rosseland group integration function can simply wrap the Planckian function
 * integratePlanckSpectrum(lowFreq, highFreq, T)
 *
 * The user is responsible for applying the appropriate constants to the result of this integration.
 * For example, to make the result of this integration equivalent to
 * \f[
 *      \int_{\nu_1}^{\nu_N} R(\nu,T) d\nu
 * \f]
 * then you must multiply by a factor of \f$\frac{4 acT^3}{4\pi}\f$ where a is the radiation
 * constant.  If you want to evaluate expressions like the following:
 *\f[
 *      \int_{4\pi} \int_{\nu_1}^{\nu_N} B(\nu,T) d\nu d\Omega
 *\f]
 * then you must multiply by \f$ 4 acT^3 \f$.
 *
 * In the limit of \f$ T \rightarrow 0, r(T) \rightarrow 0\f$, therefore we return a hard zero for a
 * temperature equal to a hard zero.
 */
/*!
 * \example cdi/test/tCDI.cc
 *
 * This test code provides an example of how to use CDI to access an user defined opacity class.  We
 * have created an opacity class called dummyOpacity that is used in the creation of a CDI object.
 * The CDI object is then used to obtain obacity data (via dummyOpacity).
 *
 * The test code also provides a mechanism to test the CDI independent of any "real" data objects.
 */
//================================================================================================//

class CDI {
  // NESTED CLASSES AND TYPEDEFS
  using SP_GrayOpacity = std::shared_ptr<const GrayOpacity>;
  using SP_MultigroupOpacity = std::shared_ptr<const MultigroupOpacity>;
  using SP_EoS = std::shared_ptr<const EoS>;
  using SP_EICoupling = std::shared_ptr<const EICoupling>;
  using SP_CPEloss = std::shared_ptr<const CPEloss>;
  using SF_GrayOpacity = std::vector<SP_GrayOpacity>;
  using VF_GrayOpacity = std::vector<SF_GrayOpacity>;
  using SF_MultigroupOpacity = std::vector<SP_MultigroupOpacity>;
  using VF_MultigroupOpacity = std::vector<SF_MultigroupOpacity>;
  using SF_CPEloss = std::vector<SP_CPEloss>;
  using std_string = std::string;
  // Typedefs for CPT mapping:
  using pt_zaid_pair = const std::pair<int32_t, int32_t>;
  using index_pt_map = std::map<pt_zaid_pair, size_t>;
  using map_it = index_pt_map::const_iterator;

  // DATA

  /*!
   * \brief Array that stores the matrix of possible GrayOpacity types.
   *
   * gray_opacities contains smart pointers that links a CDI object to a GrayOpacity object (any
   * type of gray opacity - Gandolf, EOSPAC, Analytic, etc.).  The smart pointers is entered in the
   * set functions.
   *
   * grayOpacities is indexed [0,num_Models-1][0,num_Reactions-1].  It is accessed by
   * [rtt_cdi::Model][rtt_cdi::Reaction]
   */
  VF_GrayOpacity grayOpacities;

  /*!
   * \brief Array that stores the list of possible MultigroupOpacity types.
   *
   * multigroupOpacities contains a list of smart pointers to MultigroupOpacity objects for
   * different rtt_cdi::Reaction types.  It is indexed [0,num_Models-1][0,num_Reactions-1].  It is
   * accessed by [rtt_cdi::Model][rtt_cdi::Reaction].
   */
  VF_MultigroupOpacity multigroupOpacities;

  //! Array that stores CP Eloss types.
  SF_CPEloss CPElosses;

  //! Map vector index -> particle/target pair
  index_pt_map CPEloss_map;
  /*!
   * \brief Frequency group boundaries for multigroup data.
   *
   * This is a static vector that contains the frequency boundaries for multigroup data sets.  The
   * number of frequency (energy) groups is the size of the vector minus one.
   *
   * This data is stored as static so that the same structure is guaranteed for all multigroup data
   * sets.  Thus, each CDI object will have access to the same energy group structure.
   */
  static std::vector<double> frequencyGroupBoundaries;

  /*!
   * \brief Smart pointer to the equation of state object.
   *
   * spEoS is a smart pointer that links a CDI object to an equation of state object (any type of
   * EoS - EOSPAC, Analytic, etc.).  The pointer is established in the CDI constructor.
   */
  SP_EoS spEoS;

  /*!
   * \brief Smart pointer to the electron-ion coupling object.
   *
   * spEICoupling is a smart pointer that links a CDI object to an electron-ion coupling object (any
   * type of Analytic, etc.).  The pointer is established in the CDI constructor.
   */
  SP_EICoupling spEICoupling;

  //! Material ID.
  std_string matID;

  //! Extend integration to place low and high tails in low and high groups?
  DLL_PUBLIC_cdi static bool extend;

  // IMPLELEMENTATION
  // ================

  //! Integrate the normalized Planckian from 0 to x (hnu/kT).
  inline static double integrate_planck(double const scaled_frequency);

  //! Integrate the normalized Planckian from 0 to x (hnu/kT).
  inline static double integrate_planck(double const scaled_frequency,
                                        double const exp_scaled_freqeuency);

  //! Integrate the normalized Planckian and Rosseland from 0 to x (hnu/kT)
  inline static void integrate_planck_rosseland(double const sclaed_frequency,
                                                double const exp_scaled_frequency, double &planck,
                                                double &rosseland);

public:
  // CONSTRUCTORS
  // ------------

  explicit CDI(std_string id = std::string());
  virtual ~CDI();

  // SETTERS
  // -------

  //! Register a gray opacity (rtt_cdi::GrayOpacity) with CDI.
  void setGrayOpacity(const SP_GrayOpacity &spGOp);

  //! Register a gray opacity (rtt_cdi::CPEloss) with CDI.
  void setCPEloss(const SP_CPEloss &spCPEp);

  //! Register a multigroup opacity (rtt_cdi::MultigroupOpacity) with CDI.
  void setMultigroupOpacity(const SP_MultigroupOpacity &spMGOp);

  //! Register an EOS (rtt_cdi::Eos) with CDI.
  void setEoS(const SP_EoS &in_spEoS);

  //! Register an EICoupling (rtt_cdi::EICoupling) with CDI.
  void setEICoupling(const SP_EICoupling &in_spEICoupling);

  //! Clear all data objects
  void reset();

  //! Set extended group boundaries flag
  static void setExtend() { extend = true; }

  // GETTERS
  // -------

  SP_GrayOpacity gray(rtt_cdi::Model m, rtt_cdi::Reaction r) const;
  SP_MultigroupOpacity mg(rtt_cdi::Model m, rtt_cdi::Reaction r) const;
  SP_CPEloss eloss(rtt_cdi::CPModelAngleCutoff mAC, int32_t proj_zaid, int32_t targ_zaid) const;
  SP_EoS eos() const;
  SP_EICoupling ei_coupling() const;

  //! Collapse Multigroup data to single-interval data with Planck weighting.
  static double collapseMultigroupOpacitiesPlanck(std::vector<double> const &groupBounds,
                                                  // double              const & T,
                                                  std::vector<double> const &opacity,
                                                  std::vector<double> const &planckSpectrum,
                                                  std::vector<double> &emission_group_cdf);

  //! Collapse Multigroup data to single-interval data with Planck weighting
  // (without setting the emission CDF)
  static double collapseMultigroupOpacitiesPlanck(std::vector<double> const &groupBounds,
                                                  // double              const & T,
                                                  std::vector<double> const &opacity,
                                                  std::vector<double> const &planckSpectrum);

  /*!
   * \brief Collapse Multigroup data to single-interval reciprocal data with
   *        Planck weighting.
   */
  static double
  collapseMultigroupReciprocalOpacitiesPlanck(std::vector<double> const &groupBounds,
                                              std::vector<double> const &opacity,
                                              std::vector<double> const &planckSpectrum);

  //! Collapse Multigroup data to single-interval data with Rosseland weighting.
  static double collapseMultigroupOpacitiesRosseland(std::vector<double> const &groupBounds,
                                                     std::vector<double> const &opacity,
                                                     std::vector<double> const &rosselendSpectrum);

  //! Return material ID string.
  const std_string &getMatID() const { return matID; }

  bool isGrayOpacitySet(rtt_cdi::Model, rtt_cdi::Reaction) const;
  bool isMultigroupOpacitySet(rtt_cdi::Model, rtt_cdi::Reaction) const;
  bool isCPElossSet(rtt_cdi::CPModelAngleCutoff mAC, int32_t pz, int32_t tz) const;
  bool isEoSSet() const;
  bool isEICouplingSet() const;

  //! Copies the vector of the stored frequency group boundary vector
  static std::vector<double> getFrequencyGroupBoundaries();

  //! Returns the number of frequency groups in the stored frequency vector.
  static size_t getNumberFrequencyGroups();

  //! Returns the extended group boundaries flag.
  static bool getExtend() { return extend; }

  // INTEGRATORS:
  // ===========

  // Over a frequency range:
  // -----------------------

  //! Integrate the normalized Planckian over a frequency range.
  static double integratePlanckSpectrum(double const low, double const high, double const T);

  //! Integrate the normalized Rosseland over a frequency range.
  static double integrateRosselandSpectrum(double const low, double const high, double const T);

  //! Integrate the Planckian and Rosseland over a frequency range.
  static void integrate_Rosseland_Planckian_Spectrum(double const low, double const high,
                                                     double const T, double &planck,
                                                     double &rosseland);

  // Over a specific group:
  // ---------------------

  //! Integrate the normalized Planckian spectrum over a frequency group.
  static double integratePlanckSpectrum(size_t const groupIndex, double const T);

  //! Integrate the normalized Rosseland spectrum over a frequency group
  static double integrateRosselandSpectrum(size_t const groupIndex, double const T);

  //! Integrate the Planckian and Rosseland over a frequency group.
  static void integrate_Rosseland_Planckian_Spectrum(size_t const groupIndex, double const T,
                                                     double &planck, double &rosseland);

  // Over the entire group spectrum:
  // ------------------------------

  //! Integrate the normalized Planckian spectrum over all frequency groups.
  static double integratePlanckSpectrum(double const T);

  // Over a provided vector of frequency bounds at once:
  // --------------------------------------------------

  //! Integrate the Planckian over all frequency groups
  static void integrate_Planckian_Spectrum(std::vector<double> const &bounds, double const T,
                                           std::vector<double> &planck);

  //! Integrate the Planckian over all frequency groups and return vector
  static std::vector<double> integrate_Planckian_Spectrum(std::vector<double> const &bounds,
                                                          double const T);

  //! Integrate the Rosseland over all frequency groups
  static void integrate_Rosseland_Spectrum(std::vector<double> const &bounds, double const T,
                                           std::vector<double> &rosseland);

  //! Integrate the Planckian and Rosseland over all frequency groups
  static void integrate_Rosseland_Planckian_Spectrum(std::vector<double> const &bounds,
                                                     double const T, std::vector<double> &planck,
                                                     std::vector<double> &rosseland);
};

//------------------------------------------------------------------------------------------------//
// INLINE FUNCTIONS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/**
 * \brief Integrate the normalized Planckian spectrum from 0 to \f$ x (\frac{h\nu}{kT}) \f$.
 *
 * \param scaled_freq upper integration limit, scaled by the temperature.
 *
 * \return integrated normalized Planckian from 0 to x \f$(\frac{h\nu}{kT})\f$
 */
double CDI::integrate_planck(double const scaled_freq) {
  double const exp_scaled_freq = std::exp(-scaled_freq);
  return CDI::integrate_planck(scaled_freq, exp_scaled_freq);
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Integrate the normalized Planckian spectrum from 0 to \f$ x (\frac{h\nu}{kT}) \f$.
 *
 * \param scaled_freq upper integration limit, scaled by the temperature.
 * \param exp_scaled_freq upper integration limit, scaled by an exponential function.
 * \return integrated normalized Planckian from 0 to x \f$(\frac{h\nu}{kT})\f$
 *
 * There are 3 cases to consider:
 * 1. nu/T is very large
 *    If nu/T is large enough, then the integral will be 1.0.
 * 2. nu/T is small
 *    Represent the integral via Taylor series expansion (this will be more efficient than case 3).
 * 3. All other cases. Use the polylog algorithm.
 */
double CDI::integrate_planck(double const scaled_freq, double const exp_scaled_freq) {
  Require(scaled_freq >= 0);

  // Case 1: nu/T very large -> integral == 1.0
  if (scaled_freq > 1.0e100)
    return 1.0;

  // Case 2: nu/T is sufficiently small
  // FWIW the break is at about scaled_freq < 2.06192398071289
  double const taylor = taylor_series_planck(std::min(scaled_freq, 1.0e15));
  // Case 3: all other situations
  double const poly = polylog_series_minus_one_planck(scaled_freq, exp_scaled_freq) + 1.0;

  // Choose between 2&3: For large enough nu/T, the next line will always select the polylog value.
  double const integral = std::min(taylor, poly);

  Ensure(integral >= 0.0);
  Ensure(integral <= 1.0);
  return integral;
}

//------------------------------------------------------------------------------------------------//
/*! \brief Integrate the normalized Planckian and Rosseland spectra from 0 to \f$ x
 *         (\frac{h\nu}{kT}) \f$.
 *
 * \param scaled_freq frequency upper integration limit scaled by temperature
 * \param exp_scaled_freq upper integration limit, scaled by an exponential function.
 * \param planck Variable to return the Planck integral
 * \param rosseland Variable to return the Rosseland integral
 */
void CDI::integrate_planck_rosseland(double const scaled_freq, double const exp_scaled_freq,
                                     double &planck, double &rosseland) {
  Require(scaled_freq >= 0.0);
  Require(rtt_dsxx::soft_equiv(exp_scaled_freq, std::exp(-scaled_freq)));

  // Calculate the Planckian integral
  planck = integrate_planck(scaled_freq, exp_scaled_freq);

  Require(planck >= 0.0);
  Require(planck <= 1.0);

  double const diff_rosseland = Planck2Rosseland(scaled_freq, exp_scaled_freq);

  rosseland = planck - diff_rosseland;

  Ensure(rosseland >= 0.0);
  Ensure(rosseland <= 1.0);

  return;
}

} // end namespace rtt_cdi

#endif // rtt_cdi_CDI_hh

//------------------------------------------------------------------------------------------------//
// end of cdi/CDI.hh
//------------------------------------------------------------------------------------------------//
