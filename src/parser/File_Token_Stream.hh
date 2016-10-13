//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/File_Token_Stream.hh
 * \author Kent G. Budge
 * \brief  Definition of class File_Token_Stream.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef CCS4_File_Token_Stream_HH
#define CCS4_File_Token_Stream_HH

#include "Text_Token_Stream.hh"
#include <fstream>

namespace rtt_parser {
using std::string;
using std::set;
using std::ifstream;

//-------------------------------------------------------------------------//
/*!
 * \brief File-based token stream
 *
 * File_Token_Stream represents a text token stream that derives its text
 * stream from a file in the file system.  It reports errors to the standard
 * console error stream \c cerr.
 */

class DLL_PUBLIC_parser File_Token_Stream : public Text_Token_Stream {
public:
  // CREATORS

  //! Construct an uninitialized File_Token_Stream.
  File_Token_Stream();

  //! Construct a File_Token_Stream from a file.
  explicit File_Token_Stream(string const &filename);

  //! Construct a File_Token_Stream from a file.
  File_Token_Stream(string const &filename, set<char> const &whitespace,
                    bool no_nonbreaking_ws = false);

  // MANIPULATORS

  //! Attach the File_Token_Stream to a file.
  void open(string const &filename);

  virtual void rewind();

  virtual void report(Token const &token, string const &message);

  virtual void report(string const &message);

protected:
  // IMPLEMENTATION

  virtual string location_() const;

  virtual void fill_character_buffer_();
  virtual bool error_() const;
  virtual bool end_() const;

private:
  // DATA

  string filename_; //!< File from which to take token text.
  ifstream infile_; //!< Stream from which to take token text.
};

} // rtt_parser

#endif // CCS4_File_Token_Stream_HH

//---------------------------------------------------------------------------//
// end of File_Token_Stream.hh
//---------------------------------------------------------------------------//
