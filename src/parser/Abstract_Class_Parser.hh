//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/Abstract_Class_Parser.hh
 * \author Kent Budge
 * \brief  Define class Abstract_Class_Parser
 * \note   Copyright (C) 2010-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef parser_Abstract_Class_Parser_hh
#define parser_Abstract_Class_Parser_hh

#include "parser/Contextual_Parse_Functor.hh"
#include "parser/Parse_Table.hh"
#include <functional>
#include <iostream>

namespace rtt_parser {

//================================================================================================//
/*!
 * \class Abstract_Class_Parser
 * \brief Template for parser that produces a class object.
 *
 * This template is meant to be specialized for parse tables that select one of a set of child
 * classes of a single abstract class. It simplifies and regularizes the task of allowing additional
 * child classes to be added to the table by a local developer working on his own version of one of
 * the Capsaicin drivers.
 *
 * - \a Abstract_Class The abstract class whose children are to be parsed.
 * - \a get_parse_table A function that returns a reference to the parse table for the abstract
 *   class.
 * - \a get_parsed_object A function that returns a reference to a storage location for a pointer
 *   to the abstract class.
 *
 * The key to this class is the register_child function, which is called for each child class prior
 * to attempting any parsing. It specifies a keyword for selecting each child class and a function
 * that does the actual parsing of the class specification. This assumes an input grammar of the
 * form
 *
 * \code
 * abstract class keyword
 *   child class keyword
 *     (child class specification)
 *   end
 * end
 * \endcode
 *
 * Note that Abstract_Class_Parser does not actually do any parsing itself. It is simply a
 * repository for keyword-parser combinations that is typically used by the Class_Parser for the
 * abstract class.
 *
 * See test/tstAbstract_Class_Parser.cc for an example of its use.
 *
 * This template has proven useful but does not provide a fully satisfactory solution to the problem
 * of abstract class keywords other than those specifying a child class.
 */
//================================================================================================//
template <typename Abstract_Class, Parse_Table &get_parse_table(),
          std::shared_ptr<Abstract_Class> &get_parsed_object(),
          typename Parse_Function = std::function<std::shared_ptr<Abstract_Class>(Token_Stream &)>>
class Abstract_Class_Parser {
public:
  //! Register children of the abstract class
  static void register_child(std::string const &keyword, Parse_Function parse_function);

  //! Register children of the abstract class
  static void register_child(std::string const &keyword,
                             std::shared_ptr<Abstract_Class> parse_function(Token_Stream &));

  //! Check the class invariants
  static bool check_static_class_invariants();

private:
  // IMPLEMENTATION

  //! Parse the child type
  static void parse_child_(Token_Stream &, int);

  // DATA

  //! Map of child keywords to child creation functions
  static std::vector<Parse_Function> map_;
};

// Implementation
#include "Abstract_Class_Parser.i.hh"

} // end namespace rtt_parser

#endif // parser_Abstract_Class_Parser_hh

//------------------------------------------------------------------------------------------------//
// end of parser/Abstract_Class_Parser.hh
//------------------------------------------------------------------------------------------------//
