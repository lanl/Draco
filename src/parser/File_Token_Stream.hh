//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/File_Token_Stream.hh
 * \author Kent G. Budge
 * \brief  Definition of class File_Token_Stream.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_File_Token_Stream_HH
#define rtt_File_Token_Stream_HH

#include "Text_Token_Stream.hh"
#include <fstream>

namespace rtt_parser {
using std::ifstream;
using std::set;
using std::string;

//------------------------------------------------------------------------------------------------//
/*!
 * \brief File-based token stream
 *
 * File_Token_Stream represents a text token stream that derives its text stream
 * from a file in the file system.  It reports errors to the standard console
 * error stream \c cerr.
 */
class File_Token_Stream : public Text_Token_Stream {
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

  void rewind() override;

  void report(Token const &token, string const &message) override;

  void report(string const &message) override;

  void comment(std::string const &message) override;

protected:
  // IMPLEMENTATION

  string location_() const override;

  void fill_character_buffer_() override;
  bool error_() const override;
  bool end_() const override;

  void push_include(std::string &file_name) override;
  void pop_include() override;

private:
  struct letter {
    // IMPLEMENTATION

    //! Constructor
    explicit letter(string const &file_name);

    bool check_class_invariants() const {
      return true; // none at present
    }

    //! Open the input stream.
    void open_();

    //! Rewind the stream.
    void rewind();

    // DATA
    string filename_; //!< File from which to take token text.
    ifstream infile_; //!< Stream from which to take token text.
  };

  // IMPLEMENTATION

  // DATA

  std::stack<std::shared_ptr<letter>> letters_;
  std::shared_ptr<letter> letter_;
};

} // namespace rtt_parser

#endif // rtt_File_Token_Stream_HH

//------------------------------------------------------------------------------------------------//
// end of File_Token_Stream.hh
//------------------------------------------------------------------------------------------------//
