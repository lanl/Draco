//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ipcress/IpcressGrayOpacity.hh
 * \author Kelly Thompson
 * \date   Mon Jan 22 13:23:37 2001
 * \brief  IpcressGrayOpacity class header file (derived from cdi/GrayOpacity)
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_ipcress_IpcressGrayOpacity_hh
#define rtt_cdi_ipcress_IpcressGrayOpacity_hh

#include "IpcressDataTable.hh" // IpcressDataTable object.
#include "cdi/GrayOpacity.hh"
#include "cdi/OpacityCommon.hh"
#include "ds++/Assert.hh"
#include <string>
#include <vector>

namespace rtt_cdi_ipcress {

// -------------------- //
// Forward declarations //
// -------------------- //

class IpcressFile;

//================================================================================================//
/*!
 * \class IpcressGrayOpacity
 *
 * \brief provides access to gray opacity data located in IPCRESS files.
 *
 * IpcressGrayOpacity allows the client code to retrieve opacity data for a particular material.
 * Each IpcressOpacity object represents a specific type of data defined by five attributes: an
 * IPCRESS File (via a IpcressFile object), a material identifier, an energy model (already selected
 * since this is a Gray Opacity class), a physics model and a reaction type.
 *
 * This is a concrete class derived from cdi/GrayOpacity.  This class allows to client to access the
 * data in IPCRESS files via the Ipcress libraries.
 *
 * This class is designed to be used in conjuction with the CDI.  The client code will create a
 * IpcressGrayOpacity object and use this object as an argument during the CDI instantiation.  The
 * purpose of this class is to provide a mechanism for accessing data in IPCRESS files and works by
 * calling the Ipcress library provided by X-5.  The IpcressGrayOpacity constructor expects four
 * arguments: a hook to IPCRESS data file (spIpcressFile), a material identifier, an opacity model
 * (Rosseland or Plank) and an opacity reaction specifier (total, scattering or absorption).  Once
 * constructed, this object allows the client to access any data found in the IPCRESS file for that
 * one material.  The client code will need to create a separate IpcressGrayOpacity object for each
 * material that it needs information about. Multiple opacity objects can exist per IPCRESS file.
 *
 * This class only provides access to gray opacity data.  If the user needs multigroup opacity
 * IPCRESS data he/she should use the cdi_ipcress/IpcressMultigroupOpacity class.
 *
 * When instantiated, the IpcressGrayOpacity object creates a IpcressDataTable object.  The IPCRESS
 * data is cached in this table object.  When the client requests an opacity value at a specified
 * temperature and density the IpcressGrayOpcity object calls the appropriate IPCRESS library
 * routine, which in turn, interpolates on the data cached in the IpcressDataTable object.
 *
 * \example cdi_ipcress/test/tIpcressOpacity.cc
 *
 * Example of IpcressGrayOpacity usage independent of CDI.  In this example we construct a
 * IpcressGrayOpacity object for the material Aluminum (matID=10001 in our example IPCRESS file).
 * We then use the IpcressGrayOpacity object to compute a Rosseland Gray opacity value for a
 * specified material, temperature and density.  Other forms of the getOpacity() accessor are tested
 * along with accessors that return information about the data set and the cached data table.
 *
 * \example cdi_ipcress/test/tIpcressWithCDI.cc
 *
 * This example tests and demonstrates how to use the cdi_ipcress package as a plug-in for the CDI
 * class.
 */
//================================================================================================//

class IpcressGrayOpacity : public rtt_cdi::GrayOpacity {

  // DATA

  // ----------------------- //
  // Specify unique material //
  // ----------------------- //

  //! The name of the ipcress data file.  This is saved for the packer routines
  std::string mutable ipcressFilename;

  /*!
   * \brief Identification number for one of the materials found in the IPCRESS file pointed to by
   *        spIpcressFile.
   */
  size_t materialID;

  // -------------------- //
  // Available data types //
  // -------------------- //

  // The IPCRESS file only holds specific data for each of its materials.

  //! A list of keys known by the IPCRESS file.
  std::vector<std::string> fieldNames;

  // --------------- //
  // Data specifiers //
  // --------------- //

  /*!
   * \brief The physics model that the current data set is based on. {Rosseland, Plank}.  This
   *        enumeration is defined in cdi/OpacityCommon.hh.
   */
  rtt_cdi::Model opacityModel;

  /*!
   * \brief The type of reaction rates that the current data set represents {Total, Scattering,
   *        Absorption}. This enumeration is defined in cdi/OpacityCommon.hh.
   */
  rtt_cdi::Reaction opacityReaction;

  //! A string that identifies the energy model for this class.
  std::string const energyPolicyDescriptor;

  // -------------------- //
  // Opacity lookup table //
  // -------------------- //

  /*!
   * \brief spIpcressDataTable contains a cached copy of the requested IPCRESS opacity lookup table.
   *
   * There is a one-to-one relationship between IpcressGrayOpacity and IpcressDataTable.
   */
  std::shared_ptr<IpcressDataTable const> spIpcressDataTable;

public:
  // ------------ //
  // Constructors //
  // ------------ //

  /*!
   * \brief This is the default IpcressGrayOpacity constructor.  It requires four arguments plus the
   *        energy policy (this class) to be instantiated.
   *
   * The combiniation of a data file and a material ID uniquely specifies a material.  If we add the
   * Model, Reaction and EnergyPolicy the opacity table is uniquely defined.
   *
   * \param spIpcressFile This smart pointer links an IPCRESS file (via the IpcressFile object) to a
   *     IpcressOpacity object. There may be many IpcressOpacity objects per IpcressFile object but
   *     only one IpcressFile object for each IpcressOpacity object.
   * \param in_materialID An identifier that links the IpcressOpacity object to a single material
   *     found in the specified IPCRESS file.
   * \param in_opacityModel The physics model that the current data set is based on.
   * \param in_opacityReaction The type of reaction rate that the current data set represents.
   */
  IpcressGrayOpacity(std::shared_ptr<IpcressFile const> const &spIpcressFile, size_t in_materialID,
                     rtt_cdi::Model in_opacityModel, rtt_cdi::Reaction in_opacityReaction);

  /*!
   * \brief Unpacking constructor.
   *
   * This constructor unpacks a IpcressGrayOpacity object from a state attained through the pack
   * function.
   *
   * \param packed vector<char> of packed IpcressGrayOpacity state; the packed state is attained by
   *     calling pack()
   */
  explicit IpcressGrayOpacity(std::vector<char> const &packed);

  /*!
   * \brief Default IpcressOpacity() destructor.
   *
   * This is required to correctly release memory when a IpcressGrayOpacity is destroyed.  We define
   * the destructor in the implementation file to avoid including the unnecessary header files.
   */
  ~IpcressGrayOpacity() override = default;

  //! Disable copy construction
  IpcressGrayOpacity(IpcressGrayOpacity const &rhs) = delete;

  //! Disable move construction
  IpcressGrayOpacity(IpcressGrayOpacity &&rhs) noexcept = delete;

  //! Disable assignment
  IpcressGrayOpacity &operator=(IpcressGrayOpacity const &rhs) = delete;

  //! Disable move assignment
  IpcressGrayOpacity &operator=(IpcressGrayOpacity &&rhs) noexcept = delete;

  // --------- //
  // Accessors //
  // --------- //

  /*!
   * \brief Opacity accessor that utilizes STL-like iterators.  This accessor expects a list of
   *        (temperature,density) tuples.  An opacity value will be returned for each tuple.  The
   *        temperature and density iterators are required to be the same length.  The opacity
   *        iterator should also have this same length.
   *
   * \param temperatureFirst The beginning position of a STL container that holds a list of
   *     temperatures.
   * \param temperatureLast The end position of a STL container that holds a list of temperatures.
   * \param densityFirst The beginning position of a STL container that holds a list of densities.
   * \param densityLast container that holds a list of temperatures.  \param opacityFirst The
   *     beginning position of a STL container into which opacity values corresponding to the given
   *     (temperature,density) tuple will be stored.
   * \return A list (of type OpacityIterator) of opacities are returned.  These opacities correspond
   *     to the temperature and density values provied in the two InputIterators.
   */
  template <typename TemperatureIterator, typename DensityIterator, typename OpacityIterator>
  OpacityIterator getOpacity(TemperatureIterator temperatureFirst,
                             TemperatureIterator temperatureLast, DensityIterator densityFirst,
                             DensityIterator densityLast, OpacityIterator opacityFirst) const;

  /*!
   * \brief Opacity accessor that utilizes STL-like iterators.  This accessor expects a list of
   *        temperatures in an STL container.  An opacity value will be returned for each
   *        temperature provided.  The opacity iterator should be the same length as the temperature
   *        STL container.
   *
   * \param temperatureFirst The beginning position of a STL container that holds a list of
   *     temperatures.
   * \param temperatureLast The end position of a STL container that holds a list of temperatures.
   * \param targetDensity The single density value used when computing opacities for each given
   *     temperature.
   * \param opacityFirst The beginning position of a STL container into which opacity values
   *     (corresponding to the provided temperature and density values) will be stored.
   * \return A list (of type OpacityIterator) of opacities are returned.  These opacities correspond
   *     to the temperature and density values provied.
   */
  template <typename TemperatureIterator, typename OpacityIterator>
  OpacityIterator getOpacity(TemperatureIterator temperatureFirst,
                             TemperatureIterator temperatureLast, double targetDensity,
                             OpacityIterator opacityFirst) const;

  /*!
   * \brief Opacity accessor that utilizes STL-like iterators.  This accessor expects a list of
   *     densities in an STL container.  An opacity value will be returned for each density
   *     provided.  The opacity iterator should be the same length as the density STL container.
   *
   * \param targetTemperature The single temperature value used when computing opacities for each
   *     given density.
   * \param densityFirst The beginning position of a STL container that holds a list of densities.
   * \param densityLast The end position of a STL container that holds a list of densities.
   * \param opacityFirst beginning position of a STL container into which opacity values
   *     (corresponding to the provided temperature and density values) will be stored.
   * \return A list (of type OpacityIterator) of opacities are returned.  These
   *     opacities correspond to the temperature and density values provied.
   */
  template <typename DensityIterator, typename OpacityIterator>
  OpacityIterator getOpacity(double targetTemperature, DensityIterator densityFirst,
                             DensityIterator densityLast, OpacityIterator opacityFirst) const;

  /*!
   * \brief Opacity accessor that returns a single opacity that corresponds to the provided
   *     temperature and density.
   *
   * \param targetTemperature The temperature value for which an opacity value is being requested.
   * \param targetDensity The density value for which an opacity value is being requested.
   * \return A single opacity.
   */
  double getOpacity(double targetTemperature, double targetDensity) const override;

  /*!
   * \brief Opacity accessor that returns a vector of opacities that correspond to the provided
   *     vector of temperatures and a single density value.
   *
   * \param targetTemperature A vector of temperature values for which opacity values are being
   * \param targetDensity The density value for which an opacity value is being requested.
   * \return A vector of opacities.
   */
  std::vector<double> getOpacity(std::vector<double> const &targetTemperature,
                                 double targetDensity) const override;

  /*!
   * \brief Opacity accessor that returns a vector of opacities that correspond
   *     to the provided vector of densities and a single temperature value.
   *
   * \param targetTemperature The temperature value for which an opacity value is being requested.
   * \param targetDensity A vector of density values for which opacity values are being requested.
   * \return A vector of opacities.
   */
  std::vector<double> getOpacity(double targetTemperature,
                                 std::vector<double> const &targetDensity) const override;

  //! Query whether the data is in tables or functional form.
  bool data_in_tabular_form() const override { return true; }

  //! Query to determine the reaction model.
  rtt_cdi::Reaction getReactionType() const override { return opacityReaction; }

  //! Query to determine the physics model.
  rtt_cdi::Model getModelType() const override { return opacityModel; }

  /*!
   * \brief Returns a string that describes the templated EnergyPolicy. Currently this will return
   *     either "mg" or "gray."
   */
  std::string getEnergyPolicyDescriptor() const override { return energyPolicyDescriptor; }

  /*!
   * \brief Returns a "plain English" description of the opacity data that this class references.
   *     (e.g. "Gray Rosseland Scattering".)
   *
   * The definition of this function is not included here to prevent the inclusion of the
   * IpcressFile.hh definitions within this header file.
   */
  std::string getDataDescriptor() const override { return spIpcressDataTable->getDataDescriptor(); }

  /*!
   * \brief Returns the name of the associated IPCRESS file.
   *
   * The definition of this function is not included here to prevent the inclusion of the
   * IpcressFile.hh definitions within this header file.
   */
  std::string getDataFilename() const override { return ipcressFilename; }

  /*!
   * \brief Returns a vector of temperatures that define the cached opacity data table.
   *
   * We do not return a const reference because this function must construct this information from
   * more fundamental tables.
   */
  std::vector<double> getTemperatureGrid() const override {
    return spIpcressDataTable->getTemperatures();
  }

  /*!
   * \brief Returns a vector of densities that define the cached opacity data table.
   *
   * We do not return a const reference because this function must construct this information from
   * more fundamental tables.
   */
  std::vector<double> getDensityGrid() const override { return spIpcressDataTable->getDensities(); }

  //! Returns the size of the temperature grid.
  size_t getNumTemperatures() const override { return spIpcressDataTable->getNumTemperatures(); }

  //! Returns the size of the density grid.
  size_t getNumDensities() const override { return spIpcressDataTable->getNumDensities(); }

  /*!
   * \brief Pack a IpcressGrayOpacity object.
   *
   * \return packed state in a vector<char>
   */
  std::vector<char> pack() const override;

  /*!
   * \brief Returns the general opacity model type, defined in OpacityCommon.hh. Since this is a
   *     Ipcress model, return 2 (rtt_cdi::IPCRESS_TYPE)
   */
  rtt_cdi::OpacityModelType getOpacityModelType() const override { return rtt_cdi::IPCRESS_TYPE; }

}; // end of class IpcressGrayOpacity

//------------------------------------------------------------------------------------------------//
// INCLUDE TEMPLATE MEMBER DEFINITIONS FOR AUTOMATIC TEMPLATE INSTANTIATION
//------------------------------------------------------------------------------------------------//

// --------------------------------- //
// STL-like accessors for getOpacity //
// --------------------------------- //

// ------------------------------------------ //
// getOpacity with Tuple of (T,rho) arguments //
// ------------------------------------------ //

template <typename TemperatureIterator, typename DensityIterator, typename OpacityIterator>
OpacityIterator
IpcressGrayOpacity::getOpacity(TemperatureIterator tempIter, TemperatureIterator tempLast,
                               DensityIterator densIter, DensityIterator Remember(densLast),
                               OpacityIterator opIter) const {
  // assert that the two input iterators have compatible sizes.
  Require(std::distance(tempIter, tempLast) == std::distance(densIter, densLast));

  // Loop over all (temperature,density) tuple values.
  for (; tempIter != tempLast; ++tempIter, ++densIter, ++opIter)
    // Call the Ipcress Logorithmic Interpolator for Gray data.
    *opIter = spIpcressDataTable->interpOpac(*tempIter, *densIter);
  return opIter;
}

// ------------------------------------ //
// getOpacity() with container of temps //
// ------------------------------------ //

template <typename TemperatureIterator, typename OpacityIterator>
OpacityIterator IpcressGrayOpacity::getOpacity(TemperatureIterator tempIter,
                                               TemperatureIterator tempLast, double targetDensity,
                                               OpacityIterator opIter) const {
  // loop over all the entries the temperature container and calculate an opacity value for each.
  for (; tempIter != tempLast; ++tempIter, ++opIter)
    // Call the Ipcress Logorithmic Interpolator for Gray data.
    *opIter = spIpcressDataTable->interpOpac(*tempIter, targetDensity);
  return opIter;
}

// ---------------------------------------- //
// getOpacity() with container of densities //
// ---------------------------------------- //

template <typename DensityIterator, typename OpacityIterator>
OpacityIterator IpcressGrayOpacity::getOpacity(double targetTemperature, DensityIterator densIter,
                                               DensityIterator densLast,
                                               OpacityIterator opIter) const {
  // loop over all the entries the density container and calculate an opacity value for each.
  for (; densIter != densLast; ++densIter, ++opIter)
    // Call the Ipcress Logorithmic Interpolator for Gray data.
    *opIter = spIpcressDataTable->interpOpac(targetTemperature, *densIter);
  return opIter;
}

} // end namespace rtt_cdi_ipcress

#endif // rtt_cdi_ipcress_IpcressGrayOpacity_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_ipcress/IpcressGrayOpacity.hh
//------------------------------------------------------------------------------------------------//
