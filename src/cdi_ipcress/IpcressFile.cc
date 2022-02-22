//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ipcress/IpcressFile.cc
 * \author Kelly Thompson
 * \date   Tue Aug 22 15:15:49 2000
 * \brief  Implementation file for IpcressFile class.
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "IpcressFile.hh"
#include "IpcressFile.t.hh"
#include "ds++/Assert.hh"
#include "ds++/Endian.hh"

namespace rtt_cdi_ipcress {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief The standard IpcressFile constructor.
 *
 * This is the standard IpcressFile constructor. This object is typically instantiated as a smart
 * pointer.
 *
 * \param[in] ipcressDataFilename A string that contains the name of the Ipcress data file in 
 *     IPCRESS format.  The f77 Ipcress vendor library expects a name with 80 characters or less. If
 *     the filename is longer than 80 characters the library will not be able to open the file.
 *
 * 1. Set some defaults (bytes per word, number of fields in the TOC).
 * 2. Try to open the file
 * 3. Load the title keys to verify that this is an ipcress file.
 * 4. Load the TOC.
 */
IpcressFile::IpcressFile(const std::string &ipcressDataFilename)
    : dataFilename(locateIpcressFile(ipcressDataFilename)),
      ipcress_word_size(8),            // bytes per entry in file.
      ipcressFileHandle(), toc(24, 0), // 24 records in the table of contents
      matIDs(), dfo(), ds(), materialData() {

  Require(rtt_dsxx::has_ieee_float_representation());
  Require(rtt_dsxx::fileExists(dataFilename));

  // Attempt to open the ipcress file.  Open at end of file (ate) so that we can
  // know the size of the binary file.
  ipcressFileHandle.open(dataFilename.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  Insist(ipcressFileHandle.is_open(), "IpcressFile: Unable to open ipcress file.");

  // Save the size of the file (bytes) to check against value of toc[21].
  Remember(std::ifstream::pos_type sizeOfFile(ipcressFileHandle.tellg()););

  //
  // Read the title (the first record of binary file).
  //
  // The first record has 2 8-byte words.  The first word should be either the
  // keyword 'nirvana' or possibly 'ipcress'.  The 2nd word is not used.
  size_t byte_offset(0);
  std::vector<std::string> title(2);
  read_strings(byte_offset, title);
  Insist(title[0] == "nirvana " || title[0] == "ipcress ",
         "The specified file is not IPCRESS format.");

  //
  // Read the table records from the ipcress file. See the .hh file for a
  // description of this data.
  //

  // This data block starts with the 3rd word in the file (byte 16).
  byte_offset = title.size() * ipcress_word_size;
  read_v(byte_offset, toc);

  // Checks for consistency

  // Maximum number of records allows on the file
  Remember(size_t const max_num_records = toc[1] - toc[0];);

  // Length of prefix is always 24 records.
  Check(max_num_records == toc[10] - toc[1]);
  Check(toc[2] == 24);
  Check(toc[4] == 2);
  Check(toc[5] == 0);
  Check(static_cast<size_t>(sizeOfFile) == ipcress_word_size * toc[21]);

  //
  // Read the table of memory offsets for data (dfo) and the associated table of
  // offset data sizes (ds).
  //

  // There are TOC[14] data sets in these tables.
  size_t const num_records_in_data_block = toc[14];
  dfo.resize(num_records_in_data_block);
  ds.resize(num_records_in_data_block);

  // The first entry is at location toc[0]
  byte_offset = ipcress_word_size * toc[0];
  read_v(byte_offset, ds);

  byte_offset = ipcress_word_size * toc[1];
  read_v(byte_offset, dfo);

  //
  // read a list of materials from the file:
  // - dfo[0] contains the number of materials.
  // - the memory block {dfo[1] ... dfo[0]+ds[0]} holds a list of material numbers.
  //

  byte_offset = ipcress_word_size * dfo[0];
  std::vector<int> vdata(1);
  read_v(byte_offset, vdata);
  size_t nummat = vdata[0];

  // Consistency check.  ds[0] is the total reserved space in the file for material IDs.
  Check(nummat < ds[0]);

  // Resize the list of material IDs.
  this->matIDs.resize(nummat);
  this->materialData.resize(nummat);

  // Now read the list of material IDs.
  byte_offset = ipcress_word_size * (dfo[0] + 1);
  read_v(byte_offset, this->matIDs);

  // Load the field data for each material and save to materialData[] vector.
  this->loadFieldData();

  // Close the file
  ipcressFileHandle.close();
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Indicate if the requested material id is available in the data file.
 *
 * \param[in] matid unsigned integer that specifies the material to be queried.
 */
bool IpcressFile::materialFound(size_t matid) const {
  return std::find(matIDs.begin(), matIDs.end(), matid) != matIDs.end();
}

//------------------------------------------------------------------------------------------------//
std::string IpcressFile::locateIpcressFile(std::string const &ipcressFile) {
  std::string foundFile;

  // ensure a name is provided
  Insist(ipcressFile.size() > 0,
         std::string("You must provide a filename when constructing an") + " IpcressFile object.");
  Insist(rtt_dsxx::fileExists(ipcressFile), "Could not located requested ipcress file.");

  // if the provided filename looks okay then use it.
  return ipcressFile;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Read 8 character strings from the binary file
 *
 * \param[in]  byte_offset offset into the ipcress file where the data exists.
 * \param[out] vdata       return value
 */
void IpcressFile::read_strings(size_t const byte_offset, std::vector<std::string> &vdata) const {
  Require(ipcressFileHandle.is_open());

  size_t const nitems(vdata.size());

  // temporary space for loading data from file
  std::vector<char> memblock(ipcress_word_size * nitems);

  // Move file pointer to requested location:
  ipcressFileHandle.seekg(byte_offset, std::ios::beg);

  // Read the data
  ipcressFileHandle.read(&memblock[0], ipcress_word_size * nitems);

  // Copy data into vector<string> container
  for (size_t i = 0; i < nitems; ++i)
    vdata[i] = std::string(&memblock[i * ipcress_word_size], ipcress_word_size);

  return;
}

//------------------------------------------------------------------------------------------------//
//! Populate the materialData member data container.
void IpcressFile::loadFieldData() {
  // Attempt to open the ipcress file.
  Insist(ipcressFileHandle.is_open(), "getKeys: Unable to open ipcress file.");

  // number of fields for the material (trid, rgrid, ...)
  size_t const numFields(toc[14]);

  // read the data from the file.  each entry is 24 bytes (3 words).
  size_t byte_offset = toc[10] * ipcress_word_size;
  // container for three 8-byte words.
  std::vector<std::string> thirdofkey(3);

  for (size_t i = 1; i < numFields; ++i) {
    // Note i=0 case is {'mats','fill',fill'}.  We skip this case by starting at
    // i=1.

    // read three 8-byte words.
    read_strings(byte_offset + 9 * i * ipcress_word_size, thirdofkey);
    // store result into 'keys'
    std::string field = thirdofkey[0] + thirdofkey[1] + thirdofkey[2];

    // find associated material index
    size_t id = atoi(field.c_str());
    Check(id > 0);
    size_t matid(0);
    for (size_t j = 0; j < matIDs.size(); ++j)
      if (id == matIDs[j]) {
        matid = j;
        break;
      }

    // read three 8-byte words.
    read_strings(byte_offset + 3 * (3 * i + 1) * ipcress_word_size, thirdofkey);
    // store result into 'keys'
    field = thirdofkey[0] + thirdofkey[1] + thirdofkey[2];

    // Read and save all values at this time.
    std::vector<double> value(ds[i], 0);
    read_v(dfo[i] * ipcress_word_size, value);

    materialData[matid].add_field(field, value);

    // Special treatment for comp
    if (field.substr(0, 4) == std::string("comp")) {
      double z, a;
      if (ds[i] > 0) {
        z = value[0];
        a = value[1];
        if (a <= 0.0) // support pre-7/2005 values
          a = 2.0 * z;
      } else {
        a = 1.0;
        z = 0.5;
      }
      double zoa = z / a;
      materialData[matid].set_zoa(zoa);
    }
  }
  return;
}

} // end namespace rtt_cdi_ipcress

//------------------------------------------------------------------------------------------------//
// end of IpcressFile.cc
//------------------------------------------------------------------------------------------------//
