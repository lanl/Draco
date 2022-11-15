//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/test/DummyEICoupling.hh
 * \author Mathew Cleveland
 * \date   March 2019
 * \brief  DummyEICoupling class header file (derived from ../EICoupling)
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_DummyEICoupling_hh
#define rtt_cdi_DummyEICoupling_hh

#include "cdi/EICoupling.hh"

namespace rtt_cdi_test {

//================================================================================================//
/*!
 * \class DummyEICoupling
 *
 * \brief This is an equation of state class that derives its interface from cdi/EICoupling and is
 *        used for testing purposes only.
 *
 * This electron-ion coupling class always contains the same data (set by the default constructor).
 * The data table has the following properties:
 *
 *    dummy_ei_coupling = etemperature + 10*itemperature + 100*density + 1000*w_e + 10000*w_i
 *
 * \sa cdi/test/tEICoupling.cc
 * \sa cdi/test/tCDI.cc
 */
//================================================================================================//

class DummyEICoupling : public rtt_cdi::EICoupling {
public:
  // -------------------------- //
  // Constructors & Destructors //
  // -------------------------- //

  /*!
   * \brief Constructor for DummyEICoupling object.
   *
   * The constructor assigns fixed values for all of the member data.  Every instance of this object
   * has the same member data.
   */
  DummyEICoupling() = default;

  /*!
   * \brief Default DummyEICoupling() destructor.
   *
   * This is required to correctly release memory when a DummyEICoupling object is destroyed.
   */
  ~DummyEICoupling() override = default;

  //! Disable copy construction
  DummyEICoupling(DummyEICoupling const &rhs) = delete;
  //! Disable move construction
  DummyEICoupling(DummyEICoupling &&rhs) noexcept = delete;
  //! Disable assignment
  DummyEICoupling &operator=(DummyEICoupling const &rhs) = delete;
  //! Disable move-assignment
  DummyEICoupling &operator=(DummyEICoupling &&rhs) noexcept = delete;

  // --------- //
  // Accessors //
  // --------- //

  //! EICoupling accessor.
  double getElectronIonCoupling(const double etemperature, const double itemperature,
                                const double density, const double w_e,
                                const double w_i) const override;

  //! EICoupling accessor.
  std::vector<double> getElectronIonCoupling(const std::vector<double> &etemperature,
                                             const std::vector<double> &itemperature,
                                             const std::vector<double> &density,
                                             const std::vector<double> &w_e,
                                             const std::vector<double> &w_i) const override;

  // Dummy pack function.
  std::vector<char> pack() const override { return {}; }

}; // end of class DummyEICoupling

} // end namespace rtt_cdi_test

#endif // rtt_cdi_DummyEICoupling_hh

//------------------------------------------------------------------------------------------------//
// end of cdi/test/DummyEICoupling.hh
//------------------------------------------------------------------------------------------------//
