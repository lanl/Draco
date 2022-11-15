//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/EICoupling.hh
 * \author Mathew Cleveland
 * \date   March 2019
 * \brief  EICoupling class header file (an abstract class)
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_EICoupling_hh
#define rtt_cdi_EICoupling_hh

#include "ds++/config.h"
#include <vector>

namespace rtt_cdi {

//================================================================================================//
/*!
 * \class EICoupling
 *
 * \brief This is a pure virtual class that defines a standard interface for all derived EICoupling
 *        objects.
 *
 * Any derived EICoupling object must provide as a minumum the functionality outlined in this
 * routine.  This functionality includes access to the data grid and the ability to return
 * interpolated opacity values.
 *
 * \example cdi/test/tDummyEICoupling.cc
 * \sa cdi/test/tCDI.cc
 */
//================================================================================================//
class EICoupling {
  // DATA

  // There is no data for a pure virtual object.  This class provides an interface and does not
  // preserve state.

public:
  // ---------- //
  // Destructor //
  // ---------- //

  /*!
   * \brief Default EoS() destructor.
   *
   * This is required to correctly release memory when any object derived from EICoupling is
   * destroyed.
   */
  virtual ~EICoupling() = default;

  //! default the default constructor
  EICoupling() = default;
  //! Remove the copy constructor
  EICoupling(EICoupling const &rhs) = delete;
  //! Remove the move constructor
  EICoupling(EICoupling &&rhs) noexcept = delete;
  //! Remove assignment operator
  EICoupling &operator=(EICoupling const &rhs) = delete;
  //! Remove the move-assignment operator
  EICoupling &operator=(EICoupling &&rhs) noexcept = delete;

  // --------- //
  // Accessors //
  // --------- //

  /*!
   * \brief EICoupling accessor that returns a single electron-ion coupling coefficient.
   *
   * \param eTemperature The electron temperature value for which an opacity value is being
   *     requested (keV).
   * \param iTemperature The electron temperature value for which an opacity value is being
   *     requested (keV).
   * \param density The density value for which an opacity value is being requested (g/cm^3).
   * \param w_e is the plasma electron frequency (1/s) (as defined by Eq. 3.41 in Brown, Preston,
   *     and Singleton, 'Physics Reports', V410, Issue 4, 2005)
   * \param w_i is the average plasma ion frequency (1/s) (as defined by Eq.  3.61 in Brown,
   *    Preston, and Singleton, 'Physics Reports', V410, Issue 4, 2005)
   * \return A electron-ion coupling coeffiecent (kJ/g/K/s).
   */
  virtual double getElectronIonCoupling(const double eTemperature, const double iTemperature,
                                        const double density, const double w_e,
                                        const double w_i) const = 0;

  /*!
   * \brief EICoupling accessor that returns a vector electron-ion coupling coefficients.
   *
   * \param vetemperature The electron temperature value for which an opacity value is being
   *    requested (keV).
   * \param vitemperature The ion temperature value for which an opacity value is being requested
   *     (keV).
   * \param vdensity The density value for which an opacity value is being requested (g/cm^3).
   * \param vw_e is the plasma electron frequency (1/s) (as defined by Eq. 3.41 in Brown, Preston,
   *     and Singleton, 'Physics Reports', V410, Issue 4, 2005)
   * \param vw_i is the average plasma ion frequency (1/s) (as defined by Eq.  3.61 in Brown,
   *     Preston, and Singleton, 'Physics Reports', V410, Issue 4, 2005)
   * \return A vector of electron-ion coupling coeffiecent (kJ/g/K/s).
   */
  virtual std::vector<double> getElectronIonCoupling(const std::vector<double> &vetemperature,
                                                     const std::vector<double> &vitemperature,
                                                     const std::vector<double> &vdensity,
                                                     const std::vector<double> &vw_e,
                                                     const std::vector<double> &vw_i) const = 0;

  /*!
   * \brief Interface for packing a derived EICoupling object.
   *
   * Note, the user hands the return value from this function to a derived EICoupling constructor.
   * Thus, even though one can pack a EICoupling through a base class pointer, the client must know
   * the derived type when unpacking.
   */
  virtual std::vector<char> pack() const = 0;

}; // end of class EICoupling

} // end namespace rtt_cdi

#endif // rtt_cdi_EICoupling_hh

//------------------------------------------------------------------------------------------------//
// end of cdi/EICoupling.hh
//------------------------------------------------------------------------------------------------//
