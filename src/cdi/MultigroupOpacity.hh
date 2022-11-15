//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/MultigroupOpacity.hh
 * \author Kelly Thompson
 * \date   Mon Jan 8 14:58:55 2001
 * \brief  MultigroupOpacity class header file (an abstract class)
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_MultigroupOpacity_hh
#define rtt_cdi_MultigroupOpacity_hh

#include "OpacityCommon.hh"
#include "ds++/config.h"
#include <string>
#include <vector>

namespace rtt_cdi {

//================================================================================================//
/*!
 * \class MultigroupOpacity
 *
 * \brief This is a pure virtual class that defines a standard interface for all derived
 *        MultigroupOpacity objects.
 *
 * Any derived MultigroupOpacity object must provide as a minimum the functionality outlined in this
 * routine.  This functionality includes access to the data grid and the ability to return
 * interpolated opacity values.
 *
 * \sa cdi/test/tDummyOpacity.cc
 * \sa cdi/test/tCDI.cc
 */
//================================================================================================//

class MultigroupOpacity {
  // DATA

  // There is no data for a pure virtual object.  This class provides an interface and does not
  // preserve state.

public:
  // ---------- //
  // Destructor //
  // ---------- //

  /*!
   * \brief Default Opacity() destructor.
   *
   * This is required to correctly release memory when any object derived from MultigroupOpacity is
   * destroyed.
   */
  virtual ~MultigroupOpacity() = default;
  //! Default constructor
  MultigroupOpacity() = default;
  //! Disable copy construction
  MultigroupOpacity(MultigroupOpacity const &rhs) = delete;
  //! Disable move construction
  MultigroupOpacity(MultigroupOpacity &&rhs) noexcept = delete;
  //! Disable assignment
  MultigroupOpacity &operator=(MultigroupOpacity const &rhs) = delete;
  //! Disable move-assignment
  MultigroupOpacity &operator=(MultigroupOpacity &&rhs) noexcept = delete;

  // --------- //
  // Accessors //
  // --------- //

  /*!
   * \brief Opacity accessor that returns a vector of opacities (a single opacity for each group)
   *        that correspond to the provided temperature and density.
   *
   * \param targetTemperature The temperature value for which these opacity values are being
   *          requested (keV).
   * \param targetDensity The density value for which these opacity values are being requested
   *          (g/cm^3)
   * \return A vector of opacities (a single opacity for each group).
   */
  virtual std::vector<double> getOpacity(double targetTemperature, double targetDensity) const = 0;

  /*!
   * \brief Opacity accessor that returns a vector of multigroup opacity vectors that correspond to
   *        the provided vector of temperatures and a single density value.
   *
   * \param targetTemperature A vector of temperature values for which opacity values are being
   *           requested (keV).
   *
   * \param targetDensity The density value for which an opacity value is being requested (g/cm^3).
   *
   * \return A vector of multigroup opacity vectors (cm^2/g).
   */
  virtual std::vector<std::vector<double>> getOpacity(const std::vector<double> &targetTemperature,
                                                      double targetDensity) const = 0;

  /*!
   * \brief Opacity accessor that returns a vector of multigroup opacity vectors that correspond to
   *        the provided vector of densities and a single temperature value.
   *
   * \param targetTemperature The temperature value for which an opacity value is being requested
   *           (keV).
   * \param targetDensity A vector of density values for which opacity values are being requested
   *           (g/cm^3).
   * \return A vector of multigroup opacity vectors (cm^2/g).
   */
  virtual std::vector<std::vector<double>>
  getOpacity(double targetTemperature, const std::vector<double> &targetDensity) const = 0;

  //! Query whether the data is in tables or functional form.
  virtual bool data_in_tabular_form() const = 0;

  //! Query to determine the reaction model.
  virtual rtt_cdi::Reaction getReactionType() const = 0;

  //! Query to determine the physics model.
  virtual rtt_cdi::Model getModelType() const = 0;

  /*!
   * \brief Returns a string that describes the EnergyPolicy.  Currently this will return either
   *        "mg" or "gray."
   */
  virtual std::string getEnergyPolicyDescriptor() const = 0;

  /*!
   * \brief Returns a "plain English" description of the opacity data that this class
   *        references. (e.g. "Gray Rosseland Scattering".)
   */
  virtual std::string getDataDescriptor() const = 0;

  //! Returns the name of the associated data file (if any).
  virtual std::string getDataFilename() const = 0;

  //! Returns a vector of temperatures that define the cached opacity data table.
  virtual std::vector<double> getTemperatureGrid() const = 0;

  //! Returns a vector of densities that define the cached opacity data table.
  virtual std::vector<double> getDensityGrid() const = 0;

  /*!
   * \brief Returns a vector of energy values (keV) that define the energy boundaries of the cached
   *        multigroup opacity data table.
   */
  virtual std::vector<double> getGroupBoundaries() const = 0;

  //! Returns the size of the temperature grid.
  virtual size_t getNumTemperatures() const = 0;

  //*! Returns the size of the density grid.
  virtual size_t getNumDensities() const = 0;

  //! Returns the number of group boundaries found in the current multigroup data set.
  virtual size_t getNumGroupBoundaries() const = 0;

  //! Returns the number of energy groups ( getNumGroupBoundaries() - 1 ).
  virtual size_t getNumGroups() const = 0;

  /*!
   * \brief Interface for packing a derived MultigroupOpacity object.
   *
   * Note, the user hands the return value from this function to a derived MultigroupOpacity
   * constructor.  Thus, even though one can pack a MultigroupOpacity through a base class pointer,
   * the client must know the derived type when unpacking.
   */
  virtual std::vector<char> pack() const = 0;

  //! Returns the general opacity model type (Analytic or Gandolf), defined in OpacityCommon.hh
  virtual rtt_cdi::OpacityModelType getOpacityModelType() const = 0;

}; // end of class MultigroupOpacity

} // end namespace rtt_cdi

#endif // rtt_cdi_MultigroupOpacity_hh

//------------------------------------------------------------------------------------------------//
// end of cdi/MultigroupOpacity.hh
//------------------------------------------------------------------------------------------------//
