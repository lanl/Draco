//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/test/DummyGrayOpacity.hh
 * \author Kelly Thompson
 * \date   Mon Jan 8 15:29:17 2001
 * \brief  DummyGrayOpacity class header file (derived from ../GrayOpacity)
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_DummyGrayOpacity_hh
#define rtt_cdi_DummyGrayOpacity_hh

#include "cdi/GrayOpacity.hh"
#include "cdi/OpacityCommon.hh"

namespace rtt_cdi_test {

//================================================================================================//
/*!
 * \class DummyGrayOpacity
 *
 * \brief This is an opacity class that derives its interface from cdi/GrayOpacity and is used for
 *        testing purposes only.
 *
 * \sa This opacity class always contains the same data (set by the default constructor).  The data
 *     table has the following properties:
 *
 *     Temperatures = { 1.0, 2.0, 3.0 }
 *     Densities    = { 0.1, 0.2 }
 *
 *     Opacity = temperature + density/1000
 *
 * In addition to providing definitions for the member functions outlined in GrayOpacity this class
 * provides three additional 1-D STL-like accessors for opacity data.
 *
 * \sa cdi/test/tDummyOpacity.cc
 * \sa cdi/test/tCDI.cc
 */
//================================================================================================//

class DummyGrayOpacity : public rtt_cdi::GrayOpacity {
  // DATA - all of these values are set in the constructor.

  // string descriptors
  std::string const dataFilename;           // "none"
  std::string const dataDescriptor;         // "DummyGrayOpacity"
  std::string const energyPolicyDescriptor; // "Gray"

  // data grid size
  size_t const numTemperatures = {3};
  size_t const numDensities = {2};

  // the data grid
  std::vector<double> temperatureGrid; // = { 1.0, 2.0, 3.0 }
  std::vector<double> densityGrid;     // = { 0.1, 0.2 }

  rtt_cdi::Reaction reaction_type;
  rtt_cdi::Model model_type;

public:
  // -------------------------- //
  // Constructors & Destructors //
  // -------------------------- //

  /*!
   * \brief Constructor for DummyGrayOpacity object.
   *
   * The constructor assigns fixed values for all of the member data.  Every instance of this object
   * has the same member data.
   */
  explicit DummyGrayOpacity(rtt_cdi::Reaction reaction = rtt_cdi::TOTAL,
                            rtt_cdi::Model model = rtt_cdi::ANALYTIC);

  // --------- //
  // Accessors //
  // --------- //

  /*!
   * \brief Opacity accessor that returns a single opacity that corresponds to the provided
   *         temperature and density.
   *
   * Opacity = temperature + density/1000
   *
   * \param targetTemperature The temperature value for which an opacity value is being requested
   *           (keV).
   * \param targetDensity The density value for which an opacity value is being requested (g/cm^3).
   * \return A single interpolated opacity (cm^2/g).
   */
  double getOpacity(double targetTemperature, double targetDensity) const override;

  /*!
   * \brief Opacity accessor that returns a vector of opacities that correspond to the provided
   *        vector of temperatures and a single density value.
   *
   *     Opacity[i] = temperature[i] + density/1000
   *
   * \param targetTemperature A vector of temperature values for which opacity values are being
   *           requested (keV).
   * \param targetDensity The density value for which an opacity value is being requested (g/cm^3).
   * \return A vector of opacities (cm^2/g).
   */
  std::vector<double> getOpacity(const std::vector<double> &targetTemperature,
                                 double targetDensity) const override;

  /*!
   * \brief Opacity accessor that returns a vector of opacities that correspond to the provided
   *        vector of densities and a single temperature value.
   *
   *     Opacity[i] = temperature[i] + density/1000
   *
   * \param targetTemperature The temperature value for which an opacity value is being requested
   *           (keV).
   * \param targetDensity A vector of density values for which opacity values are being requested
   *           (g/cm^3).
   * \return A vector of opacities (cm^2/g).
   */
  std::vector<double> getOpacity(double targetTemperature,
                                 const std::vector<double> &targetDensity) const override;

  //! Data is in tables.
  bool data_in_tabular_form() const override { return true; }

  //! Return the reaction type.
  rtt_cdi::Reaction getReactionType() const override { return reaction_type; }

  //! Return the model type.
  rtt_cdi::Model getModelType() const override { return model_type; }

  //! Returns a "plain English" description of the data.
  std::string getDataDescriptor() const override { return dataDescriptor; }

  /*!
   * \brief Returns a "plain English" description of the energy group structure (gray
   *        vs. multigroup).
   */
  std::string getEnergyPolicyDescriptor() const override { return energyPolicyDescriptor; }

  /*!
   * \brief Returns the name of the associated data file.  Since there is no data file associated
   *        with this opacity class the string "none" is returned.
   */
  std::string getDataFilename() const override { return dataFilename; }

  //! Returns a vector of temperatures that define the cached opacity data table.
  std::vector<double> getTemperatureGrid() const override { return temperatureGrid; }

  //*! Returns a vector of densities that define the cached opacity data table.
  std::vector<double> getDensityGrid() const override { return densityGrid; }

  //! Returns the size of the temperature grid.
  size_t getNumTemperatures() const override { return numTemperatures; }

  //! Returns the size of the density grid.
  size_t getNumDensities() const override { return numDensities; }

  //! Dummy pack function.
  std::vector<char> pack() const override { return {}; }

  //! Returns the general opacity model type, defined in OpacityCommon.hh
  rtt_cdi::OpacityModelType getOpacityModelType() const override { return rtt_cdi::DUMMY_TYPE; }
};

} // end namespace rtt_cdi_test

#endif // rtt_cdi_DummyGrayOpacity_hh

//------------------------------------------------------------------------------------------------//
// end of cdi/test/DummyGrayOpacity.hh
//------------------------------------------------------------------------------------------------//
