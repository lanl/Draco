//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   viz/Ensight_Stream.cc
 * \author Rob Lowrie
 * \brief  Ensight_Stream implementation file.
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Ensight_Stream.hh"
#include "dsxx/Assert.hh"
#include "dsxx/Packing_Utils.hh"
#include "c4/C4_Functions.hh"
#include "c4/ofpstream.hh"
#include <climits>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <ios>
#include <memory>
#include <string>

namespace rtt_viz {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief The endl manipulator.
 *
 * Note that this is a function within the rtt_viz namespace, NOT a member function of
 * Ensight_Stream.
 */
Ensight_Stream &endl(Ensight_Stream &s) {
  Require(s.d_stream);

  if (!s.d_binary)
    *s.d_stream << '\n';

  Require(s.d_stream->good());

  return s;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * This constructor opens the stream, if \a file_name is non-empty.  See open() for more
 * information.
 *
 * \param file_name  Name of output file.
 * \param binary     If true, output binary.  Otherwise, output ascii.
 * \param geom_file  If true, then a geometry file will be dumped.
 * \param decomposed If true, input is domain decomposed. Otherwise domain replicated.
 */
Ensight_Stream::Ensight_Stream(const std::string &file_name, const bool binary,
                               const bool geom_file, const bool decomposed)
    : d_decomposed_stream(), d_serial_stream(), d_binary(binary) {
  if (!file_name.empty())
    open(file_name, d_binary, geom_file, decomposed);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Destructor
 *
 * Automatically closes stream, if open.
 */
Ensight_Stream::~Ensight_Stream() { close(); }

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Opens the stream.
 *
 * \a geom_file is used only so that the "C Binary" header may be dumped when \a binary is true.  If
 * the geometry file is binary, Ensight assumes that all data files are also binary.  This class
 * does NOT check whether \a binary is consistent across all geometry and data files.
 *
 * \param file_name  Name of output file.
 * \param binary     If true, output binary.  Otherwise, output ascii.
 * \param geom_file  If true, then a geometry file will be dumped.
 * \param decomposed If true, input is domain decomposed. Otherwise domain replicated.
 */
void Ensight_Stream::open(const std::string &file_name, const bool binary, const bool geom_file,
                          const bool decomposed) {
  Require(!file_name.empty());

  d_binary = binary;

  // Open the stream.
  if (decomposed) {
    if (binary)
      d_decomposed_stream = std::make_unique<rtt_c4::ofpstream>(file_name, std::ios::binary);
    else
      d_decomposed_stream = std::make_unique<rtt_c4::ofpstream>(file_name);
    // set to a generic ostream
    d_stream = &*d_decomposed_stream;
  } else {
    Insist(rtt_c4::node() == 0, "Ensight_Stream, called by nonzero rank without "
                                "the domain decomposed flag");
    if (binary)
      d_serial_stream = std::make_unique<std::ofstream>(file_name, std::ios::binary);
    else
      d_serial_stream = std::make_unique<std::ofstream>(file_name);
    // set to a generic ostream
    d_stream = &*d_serial_stream;
  }
  Check(d_stream);

  // Set up the file.

  if (binary) {
    if (geom_file)
      *this << "C Binary";
  } else {
    // set precision for ascii mode
    d_stream->precision(5);
    d_stream->setf(std::ios::scientific, std::ios::floatfield);
  }

  Ensure(d_stream->good());
}

//------------------------------------------------------------------------------------------------//
//! Closes the stream.
void Ensight_Stream::close() {
  flush();
  if (d_decomposed_stream)
    d_decomposed_stream.reset();
  if (d_serial_stream) {
    d_serial_stream->close();
  }
  d_stream = nullptr;
}

//------------------------------------------------------------------------------------------------//
//! Output for ints.
Ensight_Stream &Ensight_Stream::operator<<(const int32_t i) {
  Require(d_stream);

  if (d_binary)
    binary_write(i);
  else
    *d_stream << std::setw(10) << i;

  Ensure(d_stream->good());

  return *this;
}

//------------------------------------------------------------------------------------------------//
/*
 * \brief Output for unsigned
 *
 * This is a convience function.  It simply casts to int.  Ensight does not support output of
 * unsigned ints.
 */
Ensight_Stream &Ensight_Stream::operator<<(const unsigned i) {
  Check(i < INT_MAX);
  int const j = static_cast<int>(i);
  *this << j;
  return *this;
}

//------------------------------------------------------------------------------------------------//
/*
 * \brief Output for int64_t.
 *
 * This is a convience function.  It simply casts to int.  Ensight does not support output of
 * unsigned ints.
 *
 * \bug Not tested so commented out.
 */
// Ensight_Stream &Ensight_Stream::operator<<(const int64_t i) {
//   Check(i < INT_MAX && i > -1 * INT_MAX);
//   int const j = static_cast<int>(i);
//   *this << j;
//   return *this;
// }

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Output for uint64_t.
 *
 * This is a convience function.  It simply casts to int.  Ensight does not support output of
 * unsigned ints.
 *
 * \bug Not tested so commented out.
 */
// Ensight_Stream &Ensight_Stream::operator<<(const uint64_t i) {
//   Check(i < INT_MAX);
//   int const j = static_cast<int>(i);
//   *this << j;
//   return *this;
// }

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Output for doubles.
 *
 * Note that Ensight only supports "float" for binary mode.
 */
Ensight_Stream &Ensight_Stream::operator<<(const double d) {
  Require(d_stream);
  if (d_binary)
    binary_write(static_cast<float>(d));
  else
    *d_stream << std::setw(12) << d;
  Ensure(d_stream->good());
  return *this;
}

//------------------------------------------------------------------------------------------------//
//! Output for strings.
Ensight_Stream &Ensight_Stream::operator<<(const std::string &s) {
  Require(d_stream);

  if (d_binary) {
    // Ensight demands all character strings be 80 chars.  Make it so.
    std::string sc(s);
    sc.resize(80);
    d_stream->write(sc.c_str(), 80);
  } else
    *d_stream << s;

  Ensure(d_stream->good());

  return *this;
}

//------------------------------------------------------------------------------------------------//
//! Output for function pointers.
Ensight_Stream &Ensight_Stream::operator<<(FP f) {
  Require(d_stream);

  Require(f);

  f(*this);

  Ensure(d_stream->good());

  return *this;
}

//------------------------------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Does binary write of \a v.
 *
 * The type \a T must support sizeof(T).
 *
 * The template implementation is defined here because only functions within this translation unit
 * should be calling this function.
 */
template <typename T> void Ensight_Stream::binary_write(const T v) {
  Require(d_stream);

  char *vc = new char[sizeof(T)];

  rtt_dsxx::Packer p;
  p.set_buffer(sizeof(T), vc);
  p.pack(v);

  d_stream->write(vc, sizeof(T));
  delete[] vc;

  Ensure(d_stream->good());
}

} // namespace rtt_viz

//------------------------------------------------------------------------------------------------//
// end of Ensight_Stream.cc
//------------------------------------------------------------------------------------------------//
