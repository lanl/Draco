//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/ifpstream.cc
 * \author Mathew Cleveland
 * \brief  ifpstream class to read processor decomposed data in parallel.
 * \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved.
 *
 * This reader is designed to read parallel decomposed data written by the ofpstream object. 
 * All ranks read from the file. We use parallel communication to read set buffer sizes from the
 * file in the order of the processor id 0->N.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/ifpstream.hh"
#include "dsxx/Assert.hh"
#include "c4/C4_Functions.hh"
#include <cstdint>
#include <ios> // ios_base::openmode
#include <string>

// cSpell: ignore bsi

namespace rtt_c4 {
using namespace std;

//------------------------------------------------------------------------------------------------//
/*! \brief Create an ifpstream to for a parallel file read.
 *
 * Create an ifpstream that reads and brocasts input file streams to all ranks
 *
 * \param[in] filename Name of the file to which is being read
 * \param[in] mode File read mode (ascii/binary)-- defaults to ascii
 */
ifpstream::ifpstream(std::string const &filename, ios_base::openmode const mode) : mode_(mode) {
  in_.open(filename, mode_);
  Insist(in_, filename + " does not exist");
}

//------------------------------------------------------------------------------------------------//
/*! \brief Fill parallel buffers with data from the input file.
 *
 *  Use buffer sizes to sync up the file read for each rank (must be called globally)
 *
 *  \param[in] buffer_size  local process buffer size (already known on each rank)
 */
void ifpstream::fill_buffers(unsigned const buffer_size) {
  using rtt_c4::global_sum;
  using rtt_c4::node;
  using rtt_c4::nodes;
  using rtt_c4::receive;
  using rtt_c4::send;

  // Store the starting position (should be the same on all ranks)
  auto start = in_.tellg();

  // Calculate offset and total read size over all ranks
  const auto bsi64 = static_cast<uint64_t>(buffer_size);
  auto offset = rtt_c4::prefix_sum(bsi64) - bsi64;
  auto total_buffer_size = bsi64;
  rtt_c4::global_sum(total_buffer_size);
  // Move to this ranks read position, read the data, and pass it to the stream
  std::string local_string(buffer_size, '\0');
  in_.seekg(offset, std::ios::cur);
  in_.read(&local_string[0], buffer_size);
  str(local_string);
  // Advance the current read pointer so all ranks match for the next read
  in_.seekg(start + static_cast<std::streamoff>(total_buffer_size), std::ios::beg);
}

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of ifpstream.cc
//------------------------------------------------------------------------------------------------//
