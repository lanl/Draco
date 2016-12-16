//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Abstract_Class_Parser.i.hh
 * \author Kent Budge
 * \date   Thu Jul 17 14:08:42 2008
 * \brief  Member definitions of class Abstract_Class_Parser
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef utils_Abstract_Class_Parser_i_hh
#define utils_Abstract_Class_Parser_i_hh

//===========================================================================//
/*!
 * \class c_string_vector
 * \brief Table of raw strings (C strings) that does not leak.
 *
 * This class exists to hold a table of pointers to raw strings, allocated by
 * strdup, that will be properly cleaned up (using free) at program termination
 * so as not to show a memory leak in memory debugging tools.
 */
class DLL_PUBLIC_parser c_string_vector {
public:
  ~c_string_vector();
  c_string_vector(void) : data(0) { /* empty */
  }
  vector<char *> data;
};

extern c_string_vector abstract_class_parser_keys;

//===========================================================================//
/*
 * The following rather lengthy and clumsy declaration declares storage for
 * the parse functions.
 *
 * Remember: typedef SP<Abstract_Class> Parse_Function(Token_Stream &);
 */
template <typename Class, Parse_Table &get_parse_table(),
          SP<Class> &get_parsed_object()>
vector<typename Abstract_Class_Parser<Class, get_parse_table,
                                      get_parsed_object>::Parse_Function *>
    Abstract_Class_Parser<Class, get_parse_table, get_parsed_object>::map_;

//---------------------------------------------------------------------------//
/*!
 * This function allows a host code to register children of the abstract class
 * with the parser. This helps support extensions by local developers.
 *
 * \param keyword Keyword associated with the child class
 *
 * \param parse_function Parse function that reads a specification from a
 * Token_Stream and returns a corresponding object of the child class.
 */
template <typename Class, Parse_Table &get_parse_table(),
          SP<Class> &get_parsed_object()>
void Abstract_Class_Parser<Class, get_parse_table, get_parsed_object>::
    register_child(string const &keyword, Parse_Function *parse_function) {
  using namespace rtt_parser;

  char *cptr = new char[keyword.size() + 1];
  std::strcpy(cptr, keyword.c_str());
  abstract_class_parser_keys.data.push_back(cptr);

  int const N = static_cast<int>(map_.size());

  map_.push_back(parse_function);

  Keyword key = {cptr, parse_child_, N, ""};

  get_parse_table().add(&key, 1);

  Ensure(check_static_class_invariants());
}

//---------------------------------------------------------------------------//
/*!
 * This is the generic parse function associated with all child keywords. It
 * makes use of the Parse_Function associated with each child keyword.
 */
template <typename Class, Parse_Table &get_parse_table(),
          SP<Class> &get_parsed_object()>
void Abstract_Class_Parser<Class, get_parse_table, get_parsed_object>::
    parse_child_(Token_Stream &tokens, int const child) {
  Check(static_cast<unsigned>(child) < map_.size());

  if (get_parsed_object() != SP<Class>()) {
    tokens.report_semantic_error("specification already exists");
  }

  get_parsed_object() = map_[child](tokens);

  Ensure(check_static_class_invariants());
}

//---------------------------------------------------------------------------//
template <typename Class, Parse_Table &get_parse_table(),
          SP<Class> &get_parsed_object()>
bool Abstract_Class_Parser<Class, get_parse_table,
                           get_parsed_object>::check_static_class_invariants() {
  return true; // no significant invariant for now
}

#endif // utils_Abstract_Class_Parser_i_hh

//---------------------------------------------------------------------------//
// end of utils/Abstract_Class_Parser.i.hh
//---------------------------------------------------------------------------//
