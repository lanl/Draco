//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/ifpstream.hh
 * \author Mathew Cleveland
 * \brief  ifpstream class to read processor decomposed data in parallel.
 * \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved.
 *
 * This is reader is designed to read parallel decomposed data written by the ofpstream object.
 * All ranks open the file and we use global communication to skip through the file in parallel.
 */
//------------------------------------------------------------------------------------------------//

#ifndef c4_ifpstream_hh
#define c4_ifpstream_hh

#include "c4/C4_Functions.hh"
#include "c4/config.h"
#include <cstddef> // size_t
#include <fstream>
#include <sstream>
#include <vector>

namespace rtt_c4 {

//================================================================================================//
/*!
 * \class ifpstream
 *
 * \brief Input stream for serializing input to all MPI ranks.
 *
 * This allows all ranks to read, but uses parallel communication to ensure that ranks read in order
 * based on their individual buffer sizes. 
 *
 * \example c4/test/tstifpstream.cc
 */
//================================================================================================//

class ifpstream : public std::stringstream { // NOLINT(misc-multiple-inheritance)
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
