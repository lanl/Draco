//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/ifpstream.hh
 * \author Mathew Cleveland
 * \date   Feb. 2020
 * \brief  ifpstream class to read processor decomposed data in parallel.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
 *
 * This is reader is designed to read parallel decomposed data written by the ofpstream object. It
 * uses similar logic by reading all data with rank==0 and broadcasting requested buffers to the
 * remaining processors.
 */
//------------------------------------------------------------------------------------------------//

#ifndef c4_ifpstream_hh
#define c4_ifpstream_hh

#include <fstream>
#include <sstream>
#include <vector>

#include "c4/config.h"

namespace rtt_c4 {

//================================================================================================//
/*!
 * \class ifpstream
 *
 * \brief Input stream for serializing input to all MPI ranks.
 *
 * This offloads all file IO to a single processor. The master processor reads all data from the
 * file and broadcasts it to the MPI processors.
 *
 * \example c4/test/tstifpstream.cc
 */
//================================================================================================//

class ifpstream : public std::stringstream {
public:
  //! Constructor -- default to standard output mode (ASCII)
  explicit ifpstream(std::string const &filename, ios_base::openmode const mode = ios_base::out);

  //! Fill buffer for every rank
  void fill_buffers(unsigned const buffer_size);

  //! Delete the default constructor
  ifpstream() = delete;

private:
  ios_base::openmode mode_; //!< specify binary or ascii input file type
  std::ifstream in_;        //!< rank 0 file stream
};

} // end namespace rtt_c4

#endif // c4_ifpstream_hh

//------------------------------------------------------------------------------------------------//
// end of c4/ifpstream.hh
//------------------------------------------------------------------------------------------------//
