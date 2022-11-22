//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/Analytic_EICoupling.hh
 * \author Mathew Cleveland
 * \date   March 2019
 * \brief  Analytic_EICoupling class definition.
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_analytic_Analytic_EICoupling_hh
#define rtt_cdi_analytic_Analytic_EICoupling_hh

#include "Analytic_Models.hh"
#include "cdi/EICoupling.hh"
#include <memory>

namespace rtt_cdi_analytic {

//================================================================================================//
/*!
 * \class Analytic_EICoupling
 *
 * \brief Derived rtt_cdi::EICoupling class for analytic electron-ion coupling data.
 *
 * The Analytic_EICoupling class is a derived rtt_cdi::EICoupling class. It provides analytic
 * rtt_cdi_analytic::Analytic_EICoupling_Model base class.  Several pre-built derived classes are
 * electron-ion coupling data.  The specific analytic EICoupling model is derived from the
 * provided in Analytic_Models.hh.
 *
 * Clients of this class can provide any analytic model class as long as it conforms to the
 * rtt_cdi_analytic::Analytic_EICoupling_Model interface.
 *
 * See the member functions for details about the data types and units.
 *
 * \example cdi_analytic/test/tstAnalytic_EICoupling.cc
 *
 * Example usage of Analytic_EICoupling, Analytic_EICoupling_Model, and their incorporation into
 * rtt_cdi::CDI.
 */
//================================================================================================//

class Analytic_EICoupling : public rtt_cdi::EICoupling {
public:
  // Useful typedefs.
  using SP_Analytic_Model = std::shared_ptr<Analytic_EICoupling_Model>;
  using const_SP_Model = std::shared_ptr<const Analytic_EICoupling_Model>;
  using sf_double = std::vector<double>;
  using sf_char = std::vector<char>;

private:
  // Analytic EICoupling model.
  SP_Analytic_Model analytic_model;

public:
  // Constructor.
  explicit Analytic_EICoupling(SP_Analytic_Model model_in);

  // Unpacking constructor.
  explicit Analytic_EICoupling(const sf_char &packed);

  // >>> ACCESSORS
  const_SP_Model get_Analytic_Model() const { return analytic_model; }

  // >>> INTERFACE SPECIFIED BY rtt_cdi::EICoupling

  // Get electron ion coupling.
  double getElectronIonCoupling(const double eTemperature, const double iTemperature,
                                const double density, const double w_e,
                                const double w_i) const override;

  sf_double getElectronIonCoupling(const sf_double &vetemperature, const sf_double &vitemperature,
                                   const sf_double &vdensity, const sf_double &vw_e,
                                   const sf_double &vw_i) const override;

  // Pack the Analytic_EICoupling into a character string.
  sf_char pack() const override;
};

} // end namespace rtt_cdi_analytic

#endif // rtt_cdi_analytic_Analytic_EICoupling_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_analytic/Analytic_EICoupling.hh
//------------------------------------------------------------------------------------------------//
