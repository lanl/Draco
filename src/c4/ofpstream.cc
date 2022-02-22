//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/ofpstream.cc
 * \author Kent G. Budge
 * \date   Mon Jun 25 11:36:43 MDT 2018
 * \brief  Define methods of class ofpstream
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ofpstream.hh"
#include "C4_Functions.hh"

namespace rtt_c4 {
using namespace std;

//------------------------------------------------------------------------------------------------//
/*! Create an ofpstream directed to a specified file.
 *
 * Create an ofpstream that synchronizes output to the specified file by MPI rank.
 *
 * \param[in] filename Name of the file to which synchronized output is to be
 *               written.
 * \param[in] mode File write mode (ascii/binary)-- defaults to ascii
 */
ofpstream::ofpstream(std::string const &filename, ios_base::openmode const mode)
    : std::ostream(&sb_) {
  sb_.mode_ = mode;
  if (rtt_c4::node() == 0) {
    sb_.out_.open(filename, mode);
  }
}

//------------------------------------------------------------------------------------------------//
/*! Synchronously write all buffered data.
 *
 * Write all buffered data to the specified file by MPI rank. That is, all buffered data for rank 0
 * is written, followed by all buffered data for rank 1, and so on.
 */
void ofpstream::mpibuf::send() {
  unsigned const pid = rtt_c4::node();
  if (pid == 0) {
    if (mode_ == ios_base::binary) {
      if (buffer_.size() > 0) {
        out_.write(&buffer_[0], buffer_.size());
      }
    } else {
      buffer_.push_back('\0'); // guarantees that buffer_.size() > 0
      out_ << &buffer_[0];
    }
    buffer_.clear();
    unsigned const pids = rtt_c4::nodes();
    for (unsigned i = 1; i < pids; ++i) {
      unsigned N(0);
      receive(&N, 1, i);
      buffer_.resize(N);
      if (N > 0) {
        rtt_c4::receive(&buffer_[0], N, i);
      }
      if (mode_ == ios_base::binary) {
        if (N > 0) {
          out_.write(&buffer_[0], buffer_.size());
        }
      } else {
        buffer_.push_back('\0'); // guarantees that buffer_.size() > 0
        out_ << &buffer_[0];
      }
    }

  } else {

    Check(buffer_.size() < UINT_MAX);
    auto N = static_cast<unsigned>(buffer_.size());
    rtt_c4::send(&N, 1, 0);
    if (N > 0)
      rtt_c4::send(&buffer_[0], N, 0);
  }
  buffer_.clear();
  out_.flush();
  rtt_c4::global_barrier();
}

//------------------------------------------------------------------------------------------------//
/*! Add the specified character to the buffer.
 *
 * For simplicity, ofpstream is currently implemented by treating every character write as an
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
/*virtual*/ ofpstream::mpibuf::int_type ofpstream::mpibuf::overflow(int_type c) {
  buffer_.push_back(static_cast<char>(c));
  return c;
}

//------------------------------------------------------------------------------------------------//
/*! Shrink the buffer to fit the current data.
 *
 * This is included for completeness, and also to let a user who is really concerned about the last
 * byte of storage shrink the buffer of an ofpstream that has done some large writes, and which he
 * will be using again later, but which he does not want tying up any memory in the meanwhile.
 */
void ofpstream::mpibuf::shrink_to_fit() { buffer_.shrink_to_fit(); }

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of ofpstream.cc
//------------------------------------------------------------------------------------------------//
