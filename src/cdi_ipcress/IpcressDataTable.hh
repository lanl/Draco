//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ipcress/IpcressDataTable.hh
 * \author Kelly Thompson
 * \date   Wednesday, Nov 16, 2011, 17:07 pm
 * \brief  Header file for IpcressDataTable
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_ipcress_IpcressDataTable_hh
#define rtt_cdi_ipcress_IpcressDataTable_hh

#include "IpcressFile.hh"
#include "cdi/OpacityCommon.hh"
#include <memory>

namespace rtt_cdi_ipcress {

//================================================================================================//
/*!
 * \brief Encapsulates all of the data associated with a specific opacity type
 *        (e.g.: total, plank, multigroup) for a single single material.
 *
 * When the user instantiates a IpcressOpacity object a IpcressDataTable object
 * is also created.  There is a one-to-one correspondence between these two
 * objects.  The IpcressDataTable object will load a single opacity table from
 * the IPCRESS file specified by the associated IpcressFile object.  The table
 * that is loaded is specified by the combination of { opacityModel,
 * opacityReaction and the opacityEnergyDescriptor }.
 *
 * Additional data about keywords and the IPCRESS format is available in
 * - Judd, B., Fontes, C.J., and Zhang, H.L. Gandolf : Interface Routines for
 *   IPCRESS Files," Los Alamos Technical Report LA-UR-01-5543, 2001.
 */
//================================================================================================//
class IpcressDataTable {

  // NESTED CLASSES AND TYPEDEFS

  // Data Descriptors:

  /*!
   * \brief A string that specifies the type of data being stored.  Possible
   *     values are rgray, ragray, rsgray, etc.  This key is provided to the
   *     Ipcress libraries as a data specifier. */
  std::string mutable ipcressDataTypeKey;

  /*!
   * \brief A string that specifies the type of data being stored.  This
   *     variables holds an English version of ipcressDataTypeKey. */
  std::string mutable dataDescriptor;

  /*!
   * \brief A string that specifies the energy model for the data being stored.
   *     Possible values are "mg" or "gray". */
  std::string const opacityEnergyDescriptor;

  /*!
   * \brief An enumerated value defined in IpcressOpacity.hh that specifies the
   *     data model.  Possible values are "Rosseland" or "Plank". */
  rtt_cdi::Model const opacityModel;

  /*!
   * \brief An enumerated valued defined in IpcressOpacity.hh that specifies the
   *     reaction model.  Possible values are "Total", "Absorption" or
   *     "Scattering". */
  rtt_cdi::Reaction const opacityReaction;

  //! A list of keys that are known by the IPCRESS file.
  std::vector<std::string> const &fieldNames;

  /*!
   * \brief The IPCRESS material number assocated with the data contained in
   *     this object. */
  size_t const matID;

  // Data Tables:

  //! The temperature grid for this data set.
  std::vector<double> mutable logTemperatures;
  std::vector<double> mutable temperatures;

  //! The density grid for this data set.
  std::vector<double> mutable logDensities;
  std::vector<double> mutable densities;

  //! The energy group boundary grid for this data set.
  std::vector<double> mutable groupBoundaries;

  //! The opacity data table.
  std::vector<double> mutable logOpacities;

public:
  // CREATORS

  //! Standard IpcressDataTable constructor.
  IpcressDataTable(std::string in_opacityEnergyDescriptor, rtt_cdi::Model in_opacityModel,
                   rtt_cdi::Reaction in_opacityReaction,
                   std::vector<std::string> const &in_fieldNames, size_t in_matID,
                   std::shared_ptr<const IpcressFile> const &spIpcressFile);

  // ACCESSORS

  //! Retrieve the size of the temperature grid.
  size_t getNumTemperatures() const { return temperatures.size(); }

  //! Retrieve the size of the density grid.
  size_t getNumDensities() const { return densities.size(); }

  //! Retrieve the size of the energy boundary grid.
  size_t getNumGroupBoundaries() const { return groupBoundaries.size(); }

  //! Retrieve the logarithmic temperature grid.
  std::vector<double> const &getTemperatures() const { return temperatures; }

  //! Retrieve the logarithmic density grid.
  std::vector<double> const &getDensities() const { return densities; }

  //! Retrieve the energy boundary grid.
  std::vector<double> const &getGroupBoundaries() const { return groupBoundaries; }

  //! Return a "plain English" description of the data table.
  std::string const &getDataDescriptor() const { return dataDescriptor; }

  //! Perform linear interploation of log(opacity) values.
  double interpOpac(double const T, double const rho, size_t const group = 0) const;

private:
  /*!
   * \brief This function sets both "ipcressDataTypeKey" and "dataDescriptor" based on the values
   *     given for opacityEnergyDescriptor, opacityModel and opacityReaction.
   */
  void setIpcressDataTypeKey() const;

  /*!
   * \brief Load the temperature, density, energy boundary and opacity opacity tables from the
   *     IPCRESS file.  Convert all tables (except energy boundaries) to log values.
   */
  void loadDataTable(std::shared_ptr<const IpcressFile> const &spIpcressFile);

  //! Search "keys" for "key".  If found return true, otherwise return false.
  template <typename T> bool key_available(const T &key, const std::vector<T> &keys) const;
};

} // namespace rtt_cdi_ipcress

#endif // rtt_cdi_ipcress_IpcressDataTable_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_ipcress/IpcressDataTable.hh
//------------------------------------------------------------------------------------------------//
