//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/opstream.cc
 * \author Kent G. Budge
 * \date   Mon Jun 25 12:12:31 MDT 2018
 * \brief  Define methods of class opstream
 * \note   Copyright (C) 2018-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "opstream.hh"
#include "C4_Functions.hh"

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
/*! Write all buffered data to console.
 *
 * Causes all buffered data to be written to console in MPI rank order; that is, all data from rank
 * 0 is written first, then all data from rank 1, and so on.
 *
 * /param[in,out] out ostream buffer to write data into. defaults to std::cout.
 */
void opstream::mpibuf::send(std::ostream &myout) {
  unsigned const pid = rtt_c4::node();
  if (pid == 0) {
    buffer_.push_back('\0'); // guarantees that buffer_.size() > 0
    myout << &buffer_[0];
    buffer_.clear();

    unsigned const pids = rtt_c4::nodes();
    for (unsigned i = 1; i < pids; ++i) {
      unsigned N(0);
      receive(&N, 1, i);
      if (N > 0) {
        buffer_.resize(N);
        rtt_c4::receive(&buffer_[0], N, i);
      }
      buffer_.push_back('\0');
      myout << &buffer_[0]; // guarantees that buffer_.size() > 0
    }
  } else {

    Check(buffer_.size() < UINT_MAX);
    auto N = static_cast<unsigned>(buffer_.size());
    rtt_c4::send(&N, 1, 0);
    if (N > 0)
      rtt_c4::send(&buffer_[0], N, 0);
  }
  buffer_.clear();
  rtt_c4::global_barrier();
}

//------------------------------------------------------------------------------------------------//
/*! Add the specified character to the buffer.
 *
 * For simplicity, opstream is currently implemented by treating every character write as an
 * overflow which is intercepted and added to the internal buffer. This is not actually that
 * inefficient for this class, since it means that when the stream using the buffer wants to insert
 * data, it checks the buffer's cursor pointer, always finds that it is null, and calls overflow
 * instead. These are not expensive operations. Should we see any evidence this class is taking
 * significant time, which should not happen for its intended use (synchronizing diagnostic output),
 * we can re-implement to let the stream do explicitly buffered insertions without this change
 * affecting any user code -- this interface is all private.
 *
 * \param[in] c Next character to add to the internal buffer.
 *
 * \return Integer representation of the character just added to the buffer.
 */
/*virtual*/ opstream::mpibuf::int_type opstream::mpibuf::overflow(int_type c) {
  buffer_.push_back(static_cast<char>(c));
  return c;
}

//------------------------------------------------------------------------------------------------//
/*! Shrink the buffer to fit the current data.
 *
 * This is included for completeness, and also to let a user who is really concerned about the last
 * byte of storage shrink the buffer of an opstream that has done some large writes, and which he
 * will be using again later, but which he does not want tying up any memory in the meanwhile.
 */
void opstream::mpibuf::shrink_to_fit() { buffer_.shrink_to_fit(); }

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of opstream.cc
//------------------------------------------------------------------------------------------------//
