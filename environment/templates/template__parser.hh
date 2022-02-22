//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   <pkg>/<class>__parser.hh
 * \author <user>
 * \brief  Define parse table for <class>
 * \note   Copyright (C) 2015-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

// clang-format off

#ifndef <spkg>_<class>__parser_hh
#define <spkg>_<class>__parser_hh

#include "<class>.hh"
#include "parser/Class_Parse_Table.hh"

namespace rtt_parser {
using<namespace>::<class>;

//================================================================================================//
template <> class Class_Parse_Table << class >> {
public:
  // NESTED CLASSES AND TYPEDEFS

  typedef<class> Return_Class;

  // CREATORS

  //! Default constructors.
  Class_Parse_Table();

  // MANIPULATORS

  // ACCESSORS

  // SERVICES

  bool allow_exit() const { return true; }

  void check_completeness(Token_Stream & tokens);

  SP << class >> create_object();

  // STATICS

  static Parse_Table const &parse_table() { return parse_table_; }

protected:
  // DATA

  // STATIC

private:
  // NESTED CLASSES AND TYPEDEFS

  // IMPLEMENTATION

  // DATA

  // STATIC

  static Class_Parse_Table *current_;
  static Parse_Table parse_table_;
};

} // end namespace rtt_parser

#endif // <spkg>_<class>__parser_hh

//------------------------------------------------------------------------------------------------//
// end of <pkg>/<class>__parser.hh
//------------------------------------------------------------------------------------------------//
