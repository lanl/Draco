//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/Analytic_Models.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 21 14:36:15 2001
 * \brief  Analytic_Models implementation file.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Analytic_Models.hh"
#include "ds++/Packing_Utils.hh"
#include "roots/zbrac.hh"
#include "roots/zbrent.hh"

namespace rtt_cdi_analytic {

//================================================================================================//
// EOS_ANALYTIC_MODEL MEMBER DEFINITIONS
//================================================================================================//

/*!
 * \brief Calculate the electron temperature given density and Electron internal energy
 *
 * \f[
 * U_e(T_i) = \int_{T=0}^{T_i}{C_v(\rho,T)dT}
 * \f]
 *
 * Where we assume \f$ U_e(0) \equiv 0 \f$.
 *
 * We have chosen to use absolute electron energy instead of dUe to mimik the behavior of EOSPAC.
 *
 * \todo Consider using GSL root finding with Newton-Raphson for improved efficiency.
 */
double Polynomial_Specific_Heat_Analytic_EoS_Model::calculate_elec_temperature(
    double const /*rho*/, double const Ue, double const Te0) const {

  // Return T=0 given Ue <= 0 every time
  if (Ue > 0.0) {
    // Set up the functor
    find_elec_temperature_functor minimizeFunctor(Ue, a, b, c);

    double const epsilon(std::numeric_limits<double>::epsilon());
    // New temperature should be nearby
    double T_max(std::max(1e-18, 100.0 * Te0)); // no zero max
    double T_min(1e-20);
    double xtol(std::min(epsilon, Te0 * epsilon));
    double ytol(Ue * std::numeric_limits<double>::epsilon());
    unsigned iterations(1000);
    // bracket the root
    rtt_roots::zbrac<find_elec_temperature_functor>(minimizeFunctor, T_min, T_max);

    // Search for the root
    double T_new = rtt_roots::zbrent<find_elec_temperature_functor>(minimizeFunctor, T_min, T_max,
                                                                    iterations, xtol, ytol);

    return T_new;
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Calculate the ion temperature given density and ion internal energy
 *
 * \f[
 * U_ic(T_n) = \int_{T=0}^{T_n}{C_v(\rho,T)dT}
 * \f]
 *
 * Where we assume \f$ U_ic(0) \equiv 0 \f$.
 *
 * We have chosen to use absolute electron energy instead of dUe to mimik the behavior of EOSPAC.
 *
 * \todo Consider using GSL root finding with Newton-Raphson for improved efficiency.
 */
double Polynomial_Specific_Heat_Analytic_EoS_Model::calculate_ion_temperature(
    double const /*rho*/, double const Uic, double const Ti0) const {

  // Return T=0 given Uic <= 0 every time
  if (Uic > 0.0) {
    // Set up the functor
    find_elec_temperature_functor minimizeFunctor(Uic, d, e, f);

    double const epsilon(std::numeric_limits<double>::epsilon());
    // New temperature should be nearby
    double T_max(std::max(1e-18, 100.0 * Ti0)); // no zero max
    double T_min(1e-20);
    double xtol(std::min(epsilon, Ti0 * epsilon));
    double ytol(Uic * std::numeric_limits<double>::epsilon());
    unsigned iterations(1000);
    // bracket the root
    rtt_roots::zbrac<find_elec_temperature_functor>(minimizeFunctor, T_min, T_max);

    // Search for the root
    double T_new = rtt_roots::zbrent<find_elec_temperature_functor>(minimizeFunctor, T_min, T_max,
                                                                    iterations, xtol, ytol);

    return T_new;
  } else {
    return 0.0;
  }
}

//================================================================================================//
// CONSTANT_ANALYTIC_MODEL MEMBER DEFINITIONS
//================================================================================================//

//! Unpacking constructor.
Constant_Analytic_Opacity_Model::Constant_Analytic_Opacity_Model(const sf_char &packed) : sigma(0) {
  // size of stream
  int size(sizeof(int) + sizeof(double));

  Require(packed.size() == static_cast<size_t>(size));

  // make an unpacker
  rtt_dsxx::Unpacker unpacker;

  // set the unpacker
  unpacker.set_buffer(size, &packed[0]);

  // unpack the indicator
  int indicator;
  unpacker >> indicator;
  Insist(indicator == CONSTANT_ANALYTIC_OPACITY_MODEL,
         "Tried to unpack the wrong type in Constant_Analytic_Opacity_Model");

  // unpack the data
  unpacker >> sigma;
  Check(sigma >= 0.0);

  Ensure(unpacker.get_ptr() == unpacker.end());
}

//------------------------------------------------------------------------------------------------//
//! Packing function
Analytic_Opacity_Model::sf_char Constant_Analytic_Opacity_Model::pack() const {
  // get the registered indicator
  int indicator = CONSTANT_ANALYTIC_OPACITY_MODEL;

  // caculate the size in bytes: indicator + 1 double
  int size = sizeof(int) + sizeof(double);

  // make a vector of the appropriate size
  sf_char pdata(size);

  // make a packer
  rtt_dsxx::Packer packer;

  // set the packer buffer
  packer.set_buffer(size, &pdata[0]);

  // pack the indicator
  packer << indicator;

  // pack the data
  packer << sigma;

  // Check the size
  Ensure(packer.get_ptr() == &pdata[0] + size);

  return pdata;
}

//------------------------------------------------------------------------------------------------//
//! Return the model parameters
Analytic_Opacity_Model::sf_double Constant_Analytic_Opacity_Model::get_parameters() const {
  return {sigma};
}

//================================================================================================//
// POLYNOMIAL_ANALYTIC_OPACITY_MODEL DEFINITIONS
//================================================================================================//

//! Unpacking constructor.
Polynomial_Analytic_Opacity_Model::Polynomial_Analytic_Opacity_Model(const sf_char &packed)
    : a(0.0), b(0.0), c(0.0), d(0.0), e(0.0), f(1.0), g(1.0), h(1.0), i(0.0), j(0.0), k(0.0) {
  // size of stream
  size_t size = sizeof(int) + 11 * sizeof(double);

  Require(packed.size() == size);

  // make an unpacker
  rtt_dsxx::Unpacker unpacker;

  // set the unpacker
  unpacker.set_buffer(size, &packed[0]);

  // unpack the indicator
  int indicator;
  unpacker >> indicator;
  Insist(indicator == POLYNOMIAL_ANALYTIC_OPACITY_MODEL,
         "Tried to unpack the wrong type in Polynomial_Analytic_Opacity_Model");

  // unpack the data
  unpacker >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j >> k;

  Ensure(unpacker.get_ptr() == unpacker.end());
}

//------------------------------------------------------------------------------------------------//
//! Packing function
Analytic_Opacity_Model::sf_char Polynomial_Analytic_Opacity_Model::pack() const {
  // get the registered indicator
  int indicator = POLYNOMIAL_ANALYTIC_OPACITY_MODEL;

  // caculate the size in bytes: indicator + 11 * double
  int size = sizeof(int) + 11 * sizeof(double);

  // make a vector of the appropriate size
  sf_char pdata(size);

  // make a packer
  rtt_dsxx::Packer packer;

  // set the packer buffer
  packer.set_buffer(size, &pdata[0]);

  // pack the indicator
  packer << indicator;

  // pack the data
  packer << a;
  packer << b;
  packer << c;
  packer << d;
  packer << e;
  packer << f;
  packer << g;
  packer << h;
  packer << i;
  packer << j;
  packer << k;

  // Check the size
  Ensure(packer.get_ptr() == &pdata[0] + size);

  return pdata;
}

//------------------------------------------------------------------------------------------------//
//! Return the model parameters
Analytic_Opacity_Model::sf_double Polynomial_Analytic_Opacity_Model::get_parameters() const {
  sf_double p(11);
  p[0] = a;
  p[1] = b;
  p[2] = c;
  p[3] = d;
  p[4] = e;
  p[5] = f;
  p[6] = g;
  p[7] = h;
  p[8] = i;
  p[9] = j;
  p[10] = k;

  return p;
}

//================================================================================================//
// POLYNOMIAL_SPECIFIC_HEAT_ANALYTIC_EOS_MODEL DEFINITIONS
//================================================================================================//

//! Unpacking constructor.
Polynomial_Specific_Heat_Analytic_EoS_Model::Polynomial_Specific_Heat_Analytic_EoS_Model(
    const sf_char &packed)
    : a(0.0), b(0.0), c(0.0), d(0.0), e(0.0), f(0.0) {
  // size of stream
  size_t size = sizeof(int) + 6 * sizeof(double);

  Require(packed.size() == size);

  // make an unpacker
  rtt_dsxx::Unpacker unpacker;

  // set the unpacker
  unpacker.set_buffer(size, &packed[0]);

  // unpack the indicator
  int indicator;
  unpacker >> indicator;
  Insist(indicator == POLYNOMIAL_SPECIFIC_HEAT_ANALYTIC_EOS_MODEL,
         "Tried to unpack the wrong type in "
         "Polynomial_Specific_Heat_Analytic_EoS_Model");

  // unpack the data
  unpacker >> a >> b >> c >> d >> e >> f;

  Ensure(unpacker.get_ptr() == unpacker.end());
}

//------------------------------------------------------------------------------------------------//
//! Packing function
Analytic_Opacity_Model::sf_char Polynomial_Specific_Heat_Analytic_EoS_Model::pack() const {
  // get the registered indicator
  int indicator = POLYNOMIAL_SPECIFIC_HEAT_ANALYTIC_EOS_MODEL;

  // caculate the size in bytes: indicator + 6 * double
  int size = sizeof(int) + 6 * sizeof(double);

  // make a vector of the appropriate size
  sf_char pdata(size);

  // make a packer
  rtt_dsxx::Packer packer;

  // set the packer buffer
  packer.set_buffer(size, &pdata[0]);

  // pack the indicator
  packer << indicator;

  // pack the data
  packer << a;
  packer << b;
  packer << c;
  packer << d;
  packer << e;
  packer << f;

  // Check the size
  Ensure(packer.get_ptr() == &pdata[0] + size);

  return pdata;
}

//------------------------------------------------------------------------------------------------//
//! Return the model parameters
Analytic_EoS_Model::sf_double Polynomial_Specific_Heat_Analytic_EoS_Model::get_parameters() const {
  sf_double p(6);
  p[0] = a;
  p[1] = b;
  p[2] = c;
  p[3] = d;
  p[4] = e;
  p[5] = f;

  return p;
}

//================================================================================================//
// CONSTANT_ANALYTIC_EI_COUPLING_MODEL MEMBER DEFINITIONS
//================================================================================================//

//! Unpacking constructor.
Constant_Analytic_EICoupling_Model::Constant_Analytic_EICoupling_Model(const sf_char &packed)
    : ei_coupling(0) {
  // size of stream
  int size(sizeof(int) + sizeof(double));

  Require(packed.size() == static_cast<size_t>(size));

  // make an unpacker
  rtt_dsxx::Unpacker unpacker;

  // set the unpacker
  unpacker.set_buffer(size, &packed[0]);

  // unpack the indicator
  int indicator;
  unpacker >> indicator;
  Insist(indicator == CONSTANT_ANALYTIC_EICOUPLING_MODEL,
         "Tried to unpack the wrong type in Constant_Analytic_Opacity_Model");

  // unpack the data
  unpacker >> ei_coupling;
  Check(ei_coupling >= 0.0);

  Ensure(unpacker.get_ptr() == unpacker.end());
}

//------------------------------------------------------------------------------------------------//
//! Packing function
Analytic_EICoupling_Model::sf_char Constant_Analytic_EICoupling_Model::pack() const {
  // get the registered indicator
  int indicator = CONSTANT_ANALYTIC_EICOUPLING_MODEL;

  // caculate the size in bytes: indicator + 1 double
  int size = sizeof(int) + sizeof(double);

  // make a vector of the appropriate size
  sf_char pdata(size);

  // make a packer
  rtt_dsxx::Packer packer;

  // set the packer buffer
  packer.set_buffer(size, &pdata[0]);

  // pack the indicator
  packer << indicator;

  // pack the data
  packer << ei_coupling;

  // Check the size
  Ensure(packer.get_ptr() == &pdata[0] + size);

  return pdata;
}

//------------------------------------------------------------------------------------------------//
//! Return the model parameters
Analytic_EICoupling_Model::sf_double Constant_Analytic_EICoupling_Model::get_parameters() const {
  return {ei_coupling};
}

//================================================================================================//
// ANALYTIC_KP_ALPHA_ELOSS_MODEL MEMBER DEFINITIONS
//================================================================================================//

/*! \brief Calculate the eloss in units of shk^-1; T given in keV, rho in g/cc, v0 in cm/shk
 *
 * The constants in this formula come directly from the fit in Eq. (2) of Kirkpatrick, R. C. and
 * Wheeler, J. A. (1981).  ``The Physics of DT Ignition In Small Fusion Targets.''  Nuclear Fusion,
 * 21(3):389–401.
 *
 * These constants are ONLY valid for alpha energy loss in DT gas.
 *
 * \param T material temperature in keV
 * \param rho material density in g/cm^3
 * \param v incident particle speed in cm/shk
 * \return eloss (vector of time coefficients) in shk^-1
 *
 */
double Analytic_KP_Alpha_Eloss_Model::calculate_eloss(const double T, const double rho,
                                                      const double v) const {
  Require(T >= 0.0);
  Require(rho >= 0.0);
  Require(v >= 0.0);

  double range = 0.03 * T * (1.0 - 0.24 * std::log(1.0 + T)) *
                 (1.0 + 0.37 * std::log((1.0 + rho) / (1.0 + 0.01 * T * T)));

  // Exponent of energy deposition term, sans minus sign and delta_t:
  double eloss = rho * v / range;

  // This analytic model can return negative elosses in some rho-T regimes; catch these and simply
  // return a large, positive number instead (equivlant to a very small, positive value of "range"
  // in the above formula)
  return (eloss >= 0.0) ? eloss : 1.0e25;
}

} // end namespace rtt_cdi_analytic

//------------------------------------------------------------------------------------------------//
// end of Analytic_Models.cc
//------------------------------------------------------------------------------------------------//
