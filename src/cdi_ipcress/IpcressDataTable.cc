//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ipcress/IpcressDataTable.cc
 * \author Kelly Thompson
 * \date   Wednesday, Nov 16, 2011, 17:04 pm
 * \brief  Implementation file for IpcressDataTable objects.
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "IpcressDataTable.hh"
#include "IpcressFile.hh"
#include "cdi/OpacityCommon.hh"
#include "ds++/Assert.hh"
#include <cmath> // we need to define log(double)

// ------------------------- //
// NAMESPACE RTT_CDI_IPCRESS //
// ------------------------- //

namespace rtt_cdi_ipcress {

// helper functions: local to file scope
double unary_log(double x) { return std::log(x); }

//------------------------------------------------------------------------------------------------//
/*!
 * \brief IpcressData Table constructor.
 *
 * The constructor requires that the data state be completely defined.  With
 * this information the DataTypeKey is set, then the data table sizes are loaded
 * and finally the table data is loaded.
 *
 * \param[in] in_opacityEnergyDescriptor This string variable specifies the
 *     energy model { "gray" or "mg" } for the opacity data contained in this
 *     IpcressDataTable object.
 * \param in_opacityModel This enumerated value specifies the physics model
 *     {Rosseland or Planck} for the opacity data contained in this object.
 *     The enumeration is defined in IpcressOpacity.hh
 * \param in_opacityReaction This enumerated value specifies the interaction
 *     model { total, scattering, absorption " for the opacity data contained
 *     in this object.  The enumeration is defined in IpcressOpacity.hh
 * \param in_fieldNames This vector of strings is a list of data keys that the
 *     IPCRESS file knows about.  This list is read from the IPCRESS file when
 *     a IpcressOpacity object is instantiated but before the associated
 *     IpcressDataTable object is created.
 * \param in_matID The material identifier that specifies a particular material
 *     in the IPCRESS file to associate with the IpcressDataTable container.
 * \param spIpcressFile A ds++ SmartPointer to a IpcressFile object.  One
 *     GanolfFile object should exist for each IPCRESS file.  Many
 *     IpcressOpacity (and thus IpcressDataTable) objects may point to the
 *     same IpcressFile object.
 */
IpcressDataTable::IpcressDataTable(std::string in_opacityEnergyDescriptor,
                                   rtt_cdi::Model in_opacityModel,
                                   rtt_cdi::Reaction in_opacityReaction,
                                   std::vector<std::string> const &in_fieldNames, size_t in_matID,
                                   std::shared_ptr<const IpcressFile> const &spIpcressFile)
    : ipcressDataTypeKey(""), dataDescriptor(""),
      opacityEnergyDescriptor(std::move(in_opacityEnergyDescriptor)), opacityModel(in_opacityModel),
      opacityReaction(in_opacityReaction), fieldNames(in_fieldNames), matID(in_matID),
      logTemperatures(), temperatures(), logDensities(), densities(), groupBoundaries(),
      logOpacities() {
  // Obtain the Ipcress keyword for the opacity data type specified by the
  // EnergyPolicy, opacityModel and the opacityReaction.  Valid keywords are:
  // {ramg, rsmg, rtmg, pmg, rgray, ragray, rsgray, pgray} This function also
  // ensures that the requested data type is available in the IPCRESS file.
  setIpcressDataTypeKey();

  // Retrieve the data set and resize the vector containers.
  temperatures = spIpcressFile->getData(matID, "tgrid");
  densities = spIpcressFile->getData(matID, "rgrid");
  groupBoundaries = spIpcressFile->getData(matID, "hnugrid");

  // Retrieve table data (temperatures, densities, group boundaries and
  // opacities.  These are stored as logorithmic values.
  loadDataTable(spIpcressFile);

} // end of IpcressDataTable constructor.

// ----------------- //
// PRIVATE FUNCTIONS //
// ----------------- //

//------------------------------------------------------------------------------------------------//
/*!
 * \brief This function sets both "ipcressDataTypeKey" and "dataDescriptor"
 *     based on the values given for opacityEnergyDescriptor, opacityModel and
 *     opacityReaction.
 */
void IpcressDataTable::setIpcressDataTypeKey() const {
  // Build the Ipcress key for the requested data.  Valid keys are: { ramg,
  // rsmg, rtmg, pmg, rgray, ragray, rsgray, pgray }

  if (opacityEnergyDescriptor == "gray") {
    switch (opacityModel) {
    case (rtt_cdi::ROSSELAND):

      switch (opacityReaction) {
      case (rtt_cdi::TOTAL):
        ipcressDataTypeKey = "rgray";
        dataDescriptor = "Gray Rosseland Total";
        break;
      case (rtt_cdi::ABSORPTION):
        ipcressDataTypeKey = "ragray";
        dataDescriptor = "Gray Rosseland Absorption";
        break;
      default:
        Assert(false);
        break;
      }
      break;

    case (rtt_cdi::PLANCK):

      switch (opacityReaction) {
      case (rtt_cdi::ABSORPTION):
        ipcressDataTypeKey = "pgray";
        dataDescriptor = "Gray Planck Absorption";
        break;
      default:
        Assert(false);
        break;
      }
      break;

    default:
      Assert(false);
      break;
    }
  } else // "mg"
  {
    switch (opacityModel) {
    case (rtt_cdi::ROSSELAND):

      switch (opacityReaction) {
      case (rtt_cdi::TOTAL):
        ipcressDataTypeKey = "rtmg";
        dataDescriptor = "Multigroup Rosseland Total";
        break;
      case (rtt_cdi::ABSORPTION):
        ipcressDataTypeKey = "ramg";
        dataDescriptor = "Multigroup Rosseland Absorption";
        break;
      case (rtt_cdi::SCATTERING):
        ipcressDataTypeKey = "rsmg";
        dataDescriptor = "Multigroup Rosseland Scattering";
        break;
      default:
        Assert(false);
        break;
      }
      break;

    case (rtt_cdi::PLANCK):

      switch (opacityReaction) {
      case (rtt_cdi::ABSORPTION):
        ipcressDataTypeKey = "pmg";
        dataDescriptor = "Multigroup Planck Absorption";
        break;
      default:
        Assert(false);
        break;
      }
      break;

    default:
      Assert(false);
      break;
    }
  }

  // Verify that the requested opacity type is available in the IPCRESS file.
  Insist(key_available(ipcressDataTypeKey, fieldNames),
         "requested opacity type is not available in the IPCRESS file.");
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Load the temperature, density, energy boundary and opacity opacity
 *     tables from the IPCRESS file.  Convert all tables (except energy
 *     boundaries) to log values.
 */
void IpcressDataTable::loadDataTable(std::shared_ptr<const IpcressFile> const &spIpcressFile) {
  // The interpolation routines expect everything to be in log form so we only
  // store the logorithmic temperature, density and opacity data.
  logTemperatures.resize(temperatures.size());
  std::transform(temperatures.begin(), temperatures.end(), logTemperatures.begin(), unary_log);
  logDensities.resize(densities.size());
  std::transform(densities.begin(), densities.end(), logDensities.begin(), unary_log);

  std::vector<double> opacities = spIpcressFile->getData(matID, ipcressDataTypeKey);
  logOpacities.resize(opacities.size());
  std::transform(opacities.begin(), opacities.end(), logOpacities.begin(), unary_log);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief This function returns "true" if "key" is found in the list of "keys".
 *        This is a static member function.
 */
template <typename T>
bool IpcressDataTable::key_available(T const &key, std::vector<T> const &keys) const {
  // Loop over all available keys.  If the requested key matches one in the list
  // return true.  If we reach the end of the list without a match return false.
  for (size_t i = 0; i < keys.size(); ++i)
    if (key == keys[i])
      return true;
  return false;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Calculate and return an interpolated opacity value.
 *
 * \param[in] targetTemperature
 * \param[in] targetDensity
 * \param[in] group Group index
 * \return An interpolated opacity value.
 *
 * \note The opacity array is a 1D array.  group id is the fastest moving index
 *       and temperatures are the slowest moving index.
 */
double IpcressDataTable::interpOpac(double const targetTemperature, double const targetDensity,
                                    size_t const group) const {
  double logT = std::log(targetTemperature);
  double logrho = std::log(targetDensity);

  size_t const numrho = logDensities.size();
  size_t const numT = logTemperatures.size();
  size_t const ng = opacityEnergyDescriptor == std::string("gray") ? 1 : groupBoundaries.size() - 1;

  // Check if we are off the table boundaries.  We don't allow extrapolation,
  // so move the target temperature or density to the table boundary.
  Check(numT > 1);
  Check(numrho > 1);
  if (targetTemperature < temperatures[0])
    logT = std::log(temperatures[0]);
  if (targetTemperature > temperatures[numT - 1])
    logT = std::log(temperatures[numT - 1]);
  if (targetDensity < densities[0])
    logrho = std::log(densities[0]);
  if (targetDensity > densities[numrho - 1])
    logrho = std::log(densities[numrho - 1]);

  /*
   * The grid looks like this:
   *
   *      |   T1     |   T      |   T2
   * -----------------------------------------
   * rho1 |   sig11  |          |   sig13
   * -----------------------------------------
   * rho  |   sig21  |  sig22   |   sig23
   * -----------------------------------------
   * rho2 |   sig31  |          |   sig33
   *
   * - rho1, rho2, T1 and T2 are table values.
   * - sig11, sig13, sig31 and sig33 are table values.
   *
   * Use linear interploation wrt log(rho) to find sig21 and sig23, then use
   * linear interpolation wrt log(T) to find sig22.
   */

  // Find the bracketing table values (T1, T2) and (rho1, rho2) for rho and T.
  size_t irho = logDensities.size() - 1;
  size_t iT = logTemperatures.size() - 1;
  for (size_t i = 0; i < numT - 1; ++i) {
    if (logT >= logTemperatures[i] && logT < logTemperatures[i + 1]) {
      iT = i;
      break;
    }
  }
  for (size_t i = 0; i < numrho - 1; ++i) {
    if (logrho >= logDensities[i] && logrho < logDensities[i + 1]) {
      irho = i;
      break;
    }
  }

  // Perform the linear interpolation.

  // index of cell with lower T and lower rho bound
  size_t i = (iT * numrho + irho) * ng + group;
  size_t k = i + ng * numrho; // index for cell with higher T value

  // If we are on the edge of the opacity table, return the edge values.  So
  // there are 4 cases:
  double logOpacity(0.0);

  // 1. Normal path
  if (irho + 1 < numrho && iT + 1 < numT) {
    double logsig12 = logOpacities[i] + (logrho - logDensities[irho]) /
                                            (logDensities[irho + 1] - logDensities[irho]) *
                                            (logOpacities[i + ng] - logOpacities[i]);

    double logsig32 = logOpacities[k] + (logrho - logDensities[irho]) /
                                            (logDensities[irho + 1] - logDensities[irho]) *
                                            (logOpacities[k + ng] - logOpacities[k]);

    logOpacity = logsig12 + (logT - logTemperatures[iT]) /
                                (logTemperatures[iT + 1] - logTemperatures[iT]) *
                                (logsig32 - logsig12);
  }

  // 2. rho is at high side of table, T is in the table
  else if (irho + 1 >= numrho && iT + 1 < numT) {
    logOpacity = logOpacities[i] + (logT - logTemperatures[iT]) /
                                       (logTemperatures[iT + 1] - logTemperatures[iT]) *
                                       (logOpacities[k] - logOpacities[i]);
  }

  // 3. T is at high side of table, rho is in the table
  else if (irho + 1 < numrho && iT + 1 >= numT) {
    logOpacity = logOpacities[i] + (logrho - logDensities[irho]) /
                                       (logDensities[irho + 1] - logDensities[irho]) *
                                       (logOpacities[i + ng] - logOpacities[i]);
  }

  // 4. Both T and rho are on the high side of the table.
  else if (irho + 1 >= numrho && iT + 1 >= numT) {
    logOpacity = logOpacities[i];
  }

  return std::exp(logOpacity);
}

} // end namespace rtt_cdi_ipcress

//------------------------------------------------------------------------------------------------//
// end of IpcressDataTable.cc
//------------------------------------------------------------------------------------------------//
