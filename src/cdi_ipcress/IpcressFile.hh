//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ipcress/IpcressFile.hh
 * \author Kelly Thompson
 * \date   Tue Aug 22 15:15:49 2000
 * \brief  Header file for IpcressFile class
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_ipcress_IpcressFile_hh
#define rtt_cdi_ipcress_IpcressFile_hh

#include "IpcressMaterial.hh"

namespace rtt_cdi_ipcress {

//================================================================================================//
/*!
 * \class IpcressFile
 *
 * \brief This class controls access the physical IPCRESS data file for
 *        IpcressOpacity. Only one IpcressFile object should exist for each data
 *        file. Several IpcressOpacity objects will access the same IpcressFile
 *        object (one per material found in the data file).
 *
 * This class is designed to be used in conjunction with IpcressOpacity. The
 * client code should create a IpcressFile object and that object is passed to
 * the IpcressOpacity constructor to create a link between the opacity object
 * and the IPCRESS data file.
 *
 * Ipcress file layout
 *
 * - all words are 8 byte entries.  [] denotes 1 word.
 * - see doc/IPCRESS_File_Format.pdf
 *
 * \code
 * offset (words)  data
 * --------------  ------------------------------------------------------
 * 0               [title][unused]
 * 2               [toc01][toc02]...[toc24]
 * toc[0]          [ds01][ds02]... (file_map[14] entries)
 * toc[1]          [dfo01][dfo02]... (file_map[14] entries)
 * toc[10]         Table of data fields, TDF
 * dfo[0]          [number of materials on file]
 * dfo[0]+1        [matid01][matid02]...           (dfo[0] entries)
 *
 * Material Data
 * (mat data for matid=10001, field=tgrid...)
 *
 * dtf[1] == {matid, field_name} == {10001,tgrid} (tgrid exists for mat10001)
 * ds[1]           tgrid has this many values.
 * dfo[1]          the tgrid values can be loaded from this address.
 *
 * dtf[2]          Mat 10001 has field rgrid
 * ds[2]           rgrid has ds[2] entries.
 * dfo[2]          rgrid's values can be loaded from this address.
 * \endcode
 *
 * \example cdi_ipcress/test/tIpcressFile.cc
 * Example of IpcressFile use independent of IpcressOpacity or CDI.
 */
//================================================================================================//

class IpcressFile {

  // NESTED CLASSES AND TYPEDEFS

  // DATA

  //! IPCRESS data filename
  std::string const dataFilename;

  //! Each value in the ipcress file uses 8-bytes.
  size_t const ipcress_word_size;

  //! File handle for the open ipcress file
  std::ifstream mutable ipcressFileHandle;

  /*!
   * A map (index) of table values
   * [0]  - disk address of 'number of words of data' array
   * [1]  - disk address of 'array of disk addresses for data' array
   *        (toc[0]+mxrec)
   * [2]  - logical length of prefix (always == 24)
   * [3]  - disk length of prefix
   * [4]  - disk address of prefix (always 2)
   * [5]  - logical length of information block (always 0)
   * [6]  - disk length of information block
   * [7]  - disk address of information block
   * [8]  - logical length of index block
   * [9]  - disk length ofindex block (3*mxkey+2)*mxrec)
   * [10] - disk address of keys ([1]+mxrec)
   * [11] - disk length of data block used.
   * [12] - disk length of data block
   * [13] - no longer used
   * [14] - number of data records in data block
   * [15] - word length of key entry for data record
   * [16] - maximum number of search keys (mxkey)
   * [17] - disk address of last index entry ([0]+[14]-1)
   * [18] - logical length of last data record
   * [19] - disk length of last data record
   * [20] - disk address of last data record
   * [21] - disk length of file
   * [22] - last address on file
   * [23] - logical data space used
   *
   * mxrec = max num records on ipcress file == [1] - [0]
   * mxkey = max num of search keys == [16]
   */
  std::vector<size_t> toc;

  //! A list of material IDs found in the data file.
  std::vector<size_t> matIDs;

  //! An array that hold disk offset to field data (e.g.: tgrid)
  std::vector<size_t> dfo;

  //! This array holds the length of each data set (how many entries in tgrid).
  std::vector<size_t> ds;

  /*!
   * \brief A vector of containers.  Each contains all field data (tgrid,
   *        ramg,...) for one material as loaded from the IPCRESS file. */
  std::vector<IpcressMaterial> materialData;

public:
  // CREATORS

  //! Standard IpcressFile constructor.
  explicit IpcressFile(std::string const &ipcressDataFilename);

  // ACCESSORS

  //! Returns the IPCRESS data filename.
  std::string const &getDataFilename() const { return dataFilename; }

  //! Returns the number of materials found in the data file.
  size_t getNumMaterials() const { return matIDs.size(); }

  //! Returns a list of material identifiers found in the data file.
  std::vector<size_t> const &getMatIDs() const { return matIDs; }

  //! Indicate if the requested material id is available in the data file.
  bool materialFound(size_t matid) const;

  /*!
   * \brief Locate the index into the materialData array for provided material
   *        identifier. */
  size_t getMatIndex(size_t const matid) const {
    size_t pos = std::find(matIDs.begin(), matIDs.end(), matid) - matIDs.begin();
    Ensure(pos < matIDs.size());
    return pos;
  }

  //! Provide a list of loaded data field names for matid.
  std::vector<std::string> listDataFieldNames(size_t const matid) const {
    Require(materialFound(matid));
    size_t matidx = getMatIndex(matid);
    return materialData[matidx].listDataFieldNames();
  }

  //! Provide access to data arrays
  std::vector<double> getData(size_t const matid, std::string const &fieldName) const {
    Require(materialFound(matid));
    size_t matidx = getMatIndex(matid);
    return materialData[matidx].data(fieldName);
  }

  //! Print a summary of the Ipcress file
  void printSummary(std::ostream &out = std::cout) const;
  void printSummary(size_t const matid, std::ostream &out = std::cout) const;

private:
  // IMPLEMENTATION

  //! Attempt to locate the requested ipcress file.
  static std::string locateIpcressFile(std::string const &ipcressFile);

  /*!
   * \brief Load the list of field data names and the associated data arrays for
   *        the requested material from the Ipcress file and save them into the
   *        IpcressMaterial container. */
  void loadFieldData();

  //! Read an array of integers or doubles from the ipcress file.
  template <typename T> void read_v(size_t const offset_bytes, std::vector<T> &vdata) const;

  //! Read strings from the binary file
  void read_strings(size_t const offset_bytes, std::vector<std::string> &vdata) const;
};

} // end namespace rtt_cdi_ipcress

#endif // rtt_cdi_ipcress_IpcressFile_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_ipcress/IpcressFile.hh
//------------------------------------------------------------------------------------------------//
