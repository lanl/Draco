//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/File_Streams.hh
 * \author Rob Lowrie
 * \date   Fri Nov 19 12:42:18 2004
 * \brief  Header for File_Output and File_Input.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_File_Streams_hh
#define rtt_dsxx_File_Streams_hh

#include "Assert.hh"
#include <array>
#include <cstring>
#include <fstream>
#include <sstream>

namespace rtt_dsxx {

//================================================================================================//
/*!
 * \class File_Output
 * \brief A file output stream.
 *
 * This class wraps \c std::ofstream and can write either ascii or binary files.  File_Output is
 * interchangeable with \c rtt_dsxx::Packer() and File_Input is interchangeable with \c
 * rtt_dsxx::Unpacker().  By templating an I/O class on stream type, the same code can be used to
 * both pack data and save data to files.
 *
 * Users of this class can also write the same code to write both ascii or binary files.  For
 * example,
 * \code
 *    void my_write(bool binary)
 *    {
 *        File_Output f("file.out", binary);
 *        double x = 1.5;
 *        int i = 3;
 *        char c = 'x';
 *        f << x << i << c;
 *    }
 * \endcode
 *  writes data to file.out as either binary or ascii, depending on the value of \a binary.  This
 *  file may be read using File_Input, as
 * \code
 *    void my_read()
 *    {
 *        File_Input f("file.out"); // automatically determines binary or not.
 *        double x;
 *        int i;
 *        char c;
 *        f >> x >> i >> c;
 *    }
 * \endcode
 *
 * Using the stream syntax (\c operator<< and \c operator>>), data written with File_Output will be
 * read properly by File_Input.  To guarantee proper reads, the ascii files are written in a certain
 * format: - Except for type char, each value is placed on its own line.  - Values of type char are
 * placed on the same line, if the previous value written was type char.  This way, character
 * strings may be written.
 *
 * Manipulators such as \c std::endl are not supported.  File_Output is not intended for "pretty
 * printing;" instead, the intent is to be able to save data to a file to be read in later by
 * File_Input.
 *
 * In binary mode, each type \a T must return the proper size from \c sizeof(T).  This restricts \c
 * File_Output and \c File_Input to types such as "Plain Old Data" (POD: int, double, char, float,
 * ...).  More complicated objects (such as \c std::string) are not supported; they must be broken
 * into their POD attributes.
 *
 * Note that binary files are generally \b not cross-platform compatible.
 *
 * \sa File_Streams.cc for detailed descriptions.
 */
//================================================================================================//

class File_Output {
private:
  // DATA

  // The stream to which data is written.
  std::ofstream d_stream;

  // If true, last datatype written was a char.  Used only in non-binary mode.
  bool d_last_was_char{false};

  // If true, in binary mode.
  bool d_binary;

public:
  // Constructor.
  explicit File_Output(const std::string &filename = "", const bool binary = false);

  // Disable other creators
  File_Output(const File_Output &) = delete;
  File_Output(const File_Output &&) = delete;
  File_Output &operator=(const File_Output &) = delete;
  File_Output &operator=(const File_Output &&) = delete;

  // Destructor.
  ~File_Output();

  // Opens filename.
  void open(const std::string &filename, const bool binary = false);

  // Closes the stream.
  void close();

  // General stream output.
  template <class T> inline File_Output &operator<<(const T i);

  // Overloaded output for type char.
  File_Output &operator<<(const char c);
};

//================================================================================================//
/*!
 * \class File_Input
 * \brief A file input stream.
 *
 * This class wraps std::ifstream and can read either ascii or binary files.  Users of this class
 * can write the same code to read both ascii or binary files.  However, binary file support
 * requires that the file be written with the File_Output stream (this restriction is so that
 * File_Input may use file header info to determine whether the file was written binary or ascii).
 * See File_Output for more information.
 *
 * \sa File_Streams.cc for detailed descriptions.
 */
//================================================================================================//

class File_Input {
private:
  // DATA

  // The stream from which data is read.
  std::ifstream d_stream;

  // The last line read from d_stream.  Used only in non-binary mode.
  std::string d_line;

  // Location within d_line to read type char.  Used only in non-binary mode.
  int d_char_line{-1};

  // If true, in binary mode.
  bool d_binary{false};

public:
  // Constructor.
  explicit File_Input(const std::string &filename = "");

  // Disable other creators (rule of 5)
  File_Input() = delete;
  File_Input(const File_Input &) = delete;
  File_Input(const File_Input &&) = delete;
  File_Input &operator=(const File_Input &) = delete;
  File_Input &operator=(const File_Input &&) = delete;

  // Destructor.
  ~File_Input();

  // Opens filename.
  void open(const std::string &filename);

  // Closes the stream.
  void close();

  // General stream input.
  template <class T> inline File_Input &operator>>(T &i);

  // Overloaded input for type char.
  File_Input &operator>>(char &c);
};

//------------------------------------------------------------------------------------------------//
//  INLINE FUNCTIONS.
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief General stream output.
 *
 * \param i The data to be written.
 * \tparam T The datatype.
 */
template <class T> File_Output &File_Output::operator<<(const T i) {
  Require(d_stream.is_open());

  if (d_binary) {
    std::array<char, sizeof(T)> buffer{};
    std::memcpy(buffer.data(), const_cast<T *>(&i), sizeof(T));
    d_stream.write(buffer.data(), sizeof(T));
  } else // ascii mode
  {
    if (d_last_was_char)
      d_stream << '\n';

    d_last_was_char = false;
    d_stream << i << '\n';
  }

  Ensure(d_stream.good());

  return *this;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief General stream input.
 *
 * \param i The data to be read.
 * \tparam T The datatype.
 */
template <class T> File_Input &File_Input::operator>>(T &i) {
  Require(d_stream.is_open());

  if (d_binary) {
    std::array<char, sizeof(T)> buffer{};
    d_stream.read(buffer.data(), sizeof(T));
    std::memcpy(&i, buffer.data(), sizeof(T));
  } else // ascii mode
  {
    std::getline(d_stream, d_line);
    Check(!d_line.empty());
    std::istringstream s(d_line);

    s >> i;
    d_char_line = -1;
  }

  Ensure(d_stream.good());

  return *this;
}

} // end namespace rtt_dsxx

#endif // rtt_dsxx_File_Streams_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/File_Streams.hh
//------------------------------------------------------------------------------------------------//
