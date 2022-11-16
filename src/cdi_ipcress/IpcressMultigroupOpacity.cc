//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ipcress/IpcressMultigroupOpacity.cc
 * \author Kelly Thompson
 * \date   Tue Nov 15 15:51:27 2011
 * \brief  IpcressMultigroupOpacity templated class implementation file.
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "IpcressMultigroupOpacity.hh"
#include "IpcressDataTable.hh"
#include "IpcressFile.hh"
#include "ds++/Assert.hh"
#include "ds++/Packing_Utils.hh"
#include <cmath>
#include <memory>

namespace rtt_cdi_ipcress {

// ------------ //
// Constructors //
// ------------ //

//------------------------------------------------------------------------------------------------//
//! Constructor for IpcressMultigroupOpacity object.
IpcressMultigroupOpacity::IpcressMultigroupOpacity(
    std::shared_ptr<IpcressFile const> const &spIpcressFile, size_t in_materialID,
    rtt_cdi::Model in_opacityModel, rtt_cdi::Reaction in_opacityReaction)
    : ipcressFilename(spIpcressFile->getDataFilename()), materialID(in_materialID), fieldNames(),
      opacityModel(in_opacityModel), opacityReaction(in_opacityReaction),
      energyPolicyDescriptor("mg"), spIpcressDataTable() {
  // Verify that the requested material ID is available in the specified IPCRESS file.
  Insist(spIpcressFile->materialFound(materialID),
         std::string("The requested material ID is not available in the ") +
             std::string("specified Ipcress file."));

  // Retrieve keys available for this material from the IPCRESS file.
  fieldNames = spIpcressFile->listDataFieldNames(materialID);
  Check(fieldNames.size() > 0);

  // Create the data table object and fill it with the table data from the IPCRESS file.
  spIpcressDataTable = std::make_shared<IpcressDataTable>(
      energyPolicyDescriptor, opacityModel, opacityReaction, fieldNames, materialID, spIpcressFile);

} // end of IpcressData constructor

//------------------------------------------------------------------------------------------------//
//! Unpacking constructor for IpcressMultigroupOpacity object.
IpcressMultigroupOpacity::IpcressMultigroupOpacity(std::vector<char> const &packed)
    : ipcressFilename(), materialID(0), fieldNames(), opacityModel(), opacityReaction(),
      energyPolicyDescriptor("mg"), spIpcressDataTable() {
  Require(packed.size() >= 5 * sizeof(int));

  // make an unpacker
  rtt_dsxx::Unpacker unpacker;
  unpacker.set_buffer(packed.size(), &packed[0]);

  // unpack and check the descriptor
  int packed_descriptor_size = 0;
  unpacker >> packed_descriptor_size;
  Check(packed_descriptor_size > 0);

  // make a vector<char> for the packed descriptor
  std::vector<char> packed_descriptor(packed_descriptor_size);

  // unpack it
  std::string descriptor;
  for (size_t i = 0; i < static_cast<size_t>(packed_descriptor_size); i++)
    unpacker >> packed_descriptor[i];
  rtt_dsxx::unpack_data(descriptor, packed_descriptor);

  // make sure it is "gray"
  Insist(descriptor == "mg", "Tried to unpack a non-mg opacity in IpcressMultigroupOpacity.");

  // unpack the size of the packed filename
  int packed_filename_size(0);
  unpacker >> packed_filename_size;

  // make a vector<char> for the packed filename
  std::vector<char> packed_filename(packed_filename_size);

  // unpack it
  for (size_t i = 0; i < static_cast<size_t>(packed_filename_size); i++)
    unpacker >> packed_filename[i];
  rtt_dsxx::unpack_data(ipcressFilename, packed_filename);

  // unpack the material id
  int itmp(0);
  unpacker >> itmp;
  materialID = static_cast<size_t>(itmp);

  // unpack the model and reaction
  int model = 0;
  int reaction = 0;
  unpacker >> model >> reaction;

  opacityModel = static_cast<rtt_cdi::Model>(model);
  opacityReaction = static_cast<rtt_cdi::Reaction>(reaction);

  Ensure(unpacker.get_ptr() == &packed[0] + packed.size());

  // build a new IpcressFile
  std::shared_ptr<IpcressFile> spIpcressFile;
  spIpcressFile = std::make_shared<IpcressFile>(ipcressFilename);
  Check(spIpcressFile);

  // Verify that the requested material ID is available in the specified IPCRESS file.
  Insist(spIpcressFile->materialFound(materialID),
         "Requested material ID is not found in the specified Ipcress file.");

  // Retrieve keys available fo this material from the IPCRESS file.
  fieldNames = spIpcressFile->listDataFieldNames(materialID);
  Check(fieldNames.size() > 0);

  // Create the data table object and fill it with the table data from the IPCRESS file.
  spIpcressDataTable = std::make_shared<IpcressDataTable>(
      energyPolicyDescriptor, opacityModel, opacityReaction, fieldNames, materialID, spIpcressFile);

  Ensure(spIpcressFile);
  Ensure(spIpcressDataTable);
}

// --------- //
// Accessors //
// --------- //

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns a single opacity (or a vector of opacities for the
 *     multigroup EnergyPolicy) that corresponds to the provided temperature and density.
 */
std::vector<double> IpcressMultigroupOpacity::getOpacity(double targetTemperature,
                                                         double targetDensity) const {
  // number of groups in this multigroup set.
  size_t const numGroups = spIpcressDataTable->getNumGroupBoundaries() - 1;

  // temporary opacity vector used by the wrapper.  The returned data will be
  // copied into the opacityIterator.
  std::vector<double> opacity(numGroups, -99.0);

  // logarithmic interpolation:
  for (size_t g = 0; g < numGroups; ++g) {
    opacity[g] = spIpcressDataTable->interpOpac(targetTemperature, targetDensity, g);
    Check(opacity[g] >= 0.0);
  }
  return opacity;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns a vector of opacities (or a vector of
 *     vectors of opacities for the multigroup EnergyPolicy) that correspond to
 *     the provided vector of temperatures and a single density value.
 */
std::vector<std::vector<double>>
IpcressMultigroupOpacity::getOpacity(std::vector<double> const &targetTemperature,
                                     double targetDensity) const {
  std::vector<std::vector<double>> opacity(targetTemperature.size());
  for (size_t i = 0; i < targetTemperature.size(); ++i)
    opacity[i] = getOpacity(targetTemperature[i], targetDensity);
  return opacity;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns a vector of opacities (or a vector of vectors of opacities
 *     for the multigroup EnergyPolicy) that correspond to the provided vector of densities and a
 *     single temperature value.
 */
std::vector<std::vector<double>>
IpcressMultigroupOpacity::getOpacity(double targetTemperature,
                                     const std::vector<double> &targetDensity) const {
  std::vector<std::vector<double>> opacity(targetDensity.size());
  for (size_t i = 0; i < targetDensity.size(); ++i)
    opacity[i] = getOpacity(targetTemperature, targetDensity[i]);
  return opacity;
}

// ------- //
// Packing //
// ------- //

//------------------------------------------------------------------------------------------------//
/*!
 * Pack the IpcressMultigroupOpacity state into a char string represented by a vector<char>. This
 * can be used for persistence, communication, etc. by accessing the char * under the vector
 * (required by implication by the standard) with the syntax &char_string[0]. Note, it is unsafe to
 * use iterators because they are \b not required to be char *.
 */
std::vector<char> IpcressMultigroupOpacity::pack() const {
  using std::string;
  using std::vector;

  // pack up the energy policy descriptor
  vector<char> packed_descriptor;
  rtt_dsxx::pack_data(energyPolicyDescriptor, packed_descriptor);

  // pack up the ipcress file name
  vector<char> packed_filename;
  rtt_dsxx::pack_data(ipcressFilename, packed_filename);

  // determine the total size: 3 ints (reaction, model, material id) + 2 ints for packed_filename
  // size and packed_descriptor size + char in packed_filename and packed_descriptor
  size_t size = 5 * sizeof(int) + packed_filename.size() + packed_descriptor.size();

  // make a container to hold packed data
  vector<char> packed(size);

  // make a packer and set it
  rtt_dsxx::Packer packer;
  packer.set_buffer(size, &packed[0]);

  // pack the descriptor
  packer << static_cast<int>(packed_descriptor.size());
  for (char &c : packed_descriptor)
    packer << c;

  // pack the filename (size and elements)
  packer << static_cast<int>(packed_filename.size());
  for (char &c : packed_filename)
    packer << c;

  // pack the material id
  packer << static_cast<int>(materialID);

  // pack the model and reaction
  packer << static_cast<int>(opacityModel) << static_cast<int>(opacityReaction);

  Ensure(packer.get_ptr() == &packed[0] + size);
  return packed;
}

} // end namespace rtt_cdi_ipcress

//------------------------------------------------------------------------------------------------//
// end of IpcressMultigroupOpacity.cc
//------------------------------------------------------------------------------------------------//
