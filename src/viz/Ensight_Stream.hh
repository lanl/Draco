//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   viz/Ensight_Stream.hh
 * \author Rob Lowrie
 * \date   Fri Nov 12 22:28:37 2004
 * \brief  Header for Ensight_Stream.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_viz_Ensight_Stream_hh
#define rtt_viz_Ensight_Stream_hh

#include "c4/ofpstream.hh"
#include "ds++/config.h"
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

namespace rtt_viz {

// Foward declarations

class Ensight_Stream;

//! A specific "endl" manipulator for Ensight_Stream.
Ensight_Stream &endl(Ensight_Stream &s);

//================================================================================================//
/*!
 * \class Ensight_Stream
 * \brief Output file stream for Ensight files.
 *
 * This class handles output to an Ensight file.  It takes care of binary or ascii mode, and the
 * proper data formatting for each mode.  The data formatting follows the Ensight Gold data format.
 * For binary mode, note that Ensight supports the following data types:
 *    - 80 character strings
 *    - float
 *    - int
 * So for example, before output, a double will be cast to a float, and a size_t will be cast to an
 * int.  Note that double floating point accuracy is not preserved by using ascii format, because
 * Ensight requires output as e12.5.
 */
//================================================================================================//

class Ensight_Stream {
private:
  // TYPEDEFS

  // FP is a function pointer.  This is usd for stream manipulators, such as rtt_viz::endl.
  using FP = Ensight_Stream &(*)(Ensight_Stream &);

  // DATA

  //! An optional parallel stream for domain decomposed geometry
  std::unique_ptr<rtt_c4::ofpstream> d_decomposed_stream;

  //! An optional serial stream for replicated domain geometry
  std::unique_ptr<std::ofstream> d_serial_stream;

  //! The standard stream for writing
  std::ostream *d_stream{};

  //! If true, in binary mode.  Otherwise, ascii mode.
  bool d_binary;

public:
  // CREATORS

  //! Constructor.
  explicit Ensight_Stream(const std::string &file_name = "", const bool binary = false,
                          const bool geom_file = false, const bool domain_decomposed = false);

  //! Destructor.
  ~Ensight_Stream();

  //! Disable copy construction
  Ensight_Stream(Ensight_Stream const &rhs) = delete;

  //! Disable move construction
  Ensight_Stream(Ensight_Stream &&rhs) noexcept = delete;

  //! Disable assignment
  Ensight_Stream &operator=(Ensight_Stream const &rhs) = delete;

  //! Disable move assignment
  Ensight_Stream &operator=(Ensight_Stream &&rhs) noexcept = delete;

  // MANIPULATORS

  //! Opens the stream.
  void open(const std::string &file_name, const bool binary = false, const bool geom_file = false,
            const bool domain_decomposed = false);

  //! Closes the stream.
  void close();

  //! Write parallel buffers
  void flush() {
    if (d_decomposed_stream)
      d_decomposed_stream->send();
    if (d_serial_stream)
      d_serial_stream->flush();
    return;
  }

  //! Expose is_open().
  bool is_open() { return bool(d_stream); }

  // The supported output stream functions.

  Ensight_Stream &operator<<(const int32_t i);
  Ensight_Stream &operator<<(const unsigned i);
  // Ensight_Stream &operator<<(const int64_t i);
  // Ensight_Stream &operator<<(const uint64_t i);
  Ensight_Stream &operator<<(const double d);
  Ensight_Stream &operator<<(const std::string &s);
  Ensight_Stream &operator<<(FP f);

  friend Ensight_Stream &endl(Ensight_Stream &s);

private:
  // Does binary write of v.
  template <typename T> void binary_write(const T v);
};

} // end namespace rtt_viz

#endif // rtt_viz_Ensight_Stream_hh

//------------------------------------------------------------------------------------------------//
// end of viz/Ensight_Stream.hh
//------------------------------------------------------------------------------------------------//
