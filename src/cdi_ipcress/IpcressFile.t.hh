//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_ipcress/IpcressFile.t.hh
 * \author Kelly Thompson
 * \date   Tue Aug 22 15:15:49 2000
 * \brief  Template Implementation file for IpcressFile class.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_ipcress_IpcressFile_t_hh__
#define __cdi_ipcress_IpcressFile_t_hh__

#include "IpcressFile.hh"
#include "ds++/Assert.hh"
#include "ds++/Endian.hh"
#include "ds++/path.hh"

namespace rtt_cdi_ipcress {

//---------------------------------------------------------------------------//
/*! 
 * \brief Read 8 character strings from the binary file
 * 
 * \param byte_offset offset into the ipcress file where the data exists.
 * \param vdata       return value 
 * \return vdata
 */
template <typename T>
void IpcressFile::read_v(size_t const byte_offset,
                         std::vector<T> &vdata) const {
  Require(ipcressFileHandle.is_open());

  size_t const nitems(vdata.size());

  // temporary space for loading data from file
  std::vector<char> memblock(ipcress_word_size * nitems);

  // Move file pointer to requested location:
  ipcressFileHandle.seekg(byte_offset, std::ios::beg);

  // Read the data
  ipcressFileHandle.read(&memblock[0], ipcress_word_size * nitems);

  // Copy data into vector<int> container
  double ddata;
  for (size_t i = 0; i < nitems; ++i) {
    // cast raw cahr data to double and perform a byte swap
    std::memcpy(&ddata, &memblock[i * ipcress_word_size], ipcress_word_size);
    if (!rtt_dsxx::is_big_endian())
      rtt_dsxx::byte_swap(ddata);
    // Save to the vector<int>
    vdata[i] = static_cast<T>(ddata);
  }
  return;
}

} // end namespace rtt_cdi_ipcress

#endif // __cdi_ipcress_IpcressFile_t_hh__

//---------------------------------------------------------------------------//
// end of cdi_ipcress/IpcressFile.t.hh
//---------------------------------------------------------------------------//
