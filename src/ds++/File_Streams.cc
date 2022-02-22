//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/File_Streams.cc
 * \author Rob Lowrie
 * \date   Mon Nov 15 10:03:51 2004
 * \brief  File_Streams implementation file.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "File_Streams.hh"
#include <iomanip>

namespace {

// Define a string to indicate that a file was written in binary mode. This string should be one
// that is unlikely to be used by a client.
static const std::string BINARY_FILE_HEADER = "bInArYfIlE_rtt_dsxx_File_Streams";
} // namespace

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
// File_Output functions.
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param filename The file name to open for writing.  If empty, open() must be used later to open a
 *                 file.
 * \param binary   If true, use binary mode for writing.
 */
File_Output::File_Output(std::string const &filename, bool const binary)
    : d_stream(), d_binary(binary)

{
  if (!filename.empty())
    open(filename, binary);
}

//------------------------------------------------------------------------------------------------//
//! Destructor
File_Output::~File_Output() { close(); }

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Opens a file for writing.
 *
 * \param filename The file name to open for writing.
 * \param binary   If true, use binary mode for writing.
 */
void File_Output::open(const std::string &filename, const bool binary) {
  Require(!filename.empty());

  if (d_stream.is_open())
    close();

  d_last_was_char = false;
  d_binary = binary;

  if (d_binary) {
    d_stream.open(filename.c_str(), std::ios::binary);
    d_stream << BINARY_FILE_HEADER;
  } else {
    d_stream.open(filename.c_str());
  }

  Ensure(d_stream);
}

//------------------------------------------------------------------------------------------------//
//! Closes the stream.
void File_Output::close() {
  if (d_stream.is_open()) {
    if (d_last_was_char)
      d_stream << std::endl;
    d_stream.close();
  }
}

//------------------------------------------------------------------------------------------------//
//! Stream output for type char.
File_Output &File_Output::operator<<(const char c) {
  Require(d_stream.is_open());

  if (d_binary) {
    d_stream.write(&c, 1);
  } else // ascii mode
  {
    // For char, we don't add a newline, in case its part of a character string.
    d_last_was_char = true;
    d_stream << c;
  }

  Ensure(d_stream.good());

  return *this;
}

//------------------------------------------------------------------------------------------------//
// File_Input functions.
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param filename The file name to open for reading.  If empty, open() must be used later to open a
 *                 file.
 */
File_Input::File_Input(std::string const &filename) : d_stream(), d_line(std::string()) {
  if (!filename.empty())
    open(filename);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Destructor
 */
File_Input::~File_Input() { close(); }

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Opens a file for reading.
 *
 * \param filename The file name to open for reading.
 */
void File_Input::open(const std::string &filename) {
  Require(!filename.empty());

  using std::string;

  d_char_line = -1;

  if (d_stream.is_open())
    close();

  // Start by opening the file in binary mode.

  d_stream.open(filename.c_str(), std::ios::binary);
  d_binary = true;

  // Check if the binary header is present.

  for (char s : BINARY_FILE_HEADER) {
    char c('x');
    d_stream >> c;
    if (c != s || (!d_stream.good())) {
      d_binary = false;
      break;
    }
  }

  // If the file is not binary, re-open in ascii mode.
  if (!d_binary) {
    d_stream.close();
    d_stream.open(filename.c_str());
  }

  Ensure(d_stream);
}

//------------------------------------------------------------------------------------------------//
//! Closes the stream.
void File_Input::close() {
  if (d_stream.is_open()) {
    d_stream.close();
  }
}

//------------------------------------------------------------------------------------------------//
//! Stream input for type char.
File_Input &File_Input::operator>>(char &c) {
  Require(d_stream.is_open());

  if (d_binary) {
    d_stream.read(&c, 1);
  } else // ascii mode
  {
    if (d_char_line < 0) {
      std::getline(d_stream, d_line);
      Check(!d_line.empty());
      d_char_line = 0;
    }

    Check(static_cast<size_t>(d_char_line) < d_line.size());
    c = d_line[d_char_line];
    ++d_char_line;
  }

  Ensure(d_stream.good());

  return *this;
}

} // namespace rtt_dsxx

//------------------------------------------------------------------------------------------------//
// end of ds++/File_Streams.cc
//------------------------------------------------------------------------------------------------//
