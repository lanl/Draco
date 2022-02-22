//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   <pkg>/<class>__parser.cc
 * \author <user>
 * \brief  Define parse table for <class>
 * \note   Copyright (C) 2015-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "<class>__parser.hh"
#include "parser/utilities.hh"

namespace rtt_parser {
using namespace rtt_dsxx;
using namespace<namespace>;

Class_Parse_Table << class >> *Class_Parse_Table << class >> ::current_;
Parse_Table Class_Parse_Table << class >> ::parse_table_;

//------------------------------------------------------------------------------------------------//
Class_Parse_Table << class >> ::Class_Parse_Table() {
  static bool first_time = true;
  if (first_time) {
    const Keyword keywords[] = {
        // {"sample", parse_sample, 0, ""},
    };

    const unsigned number_of_keywords = sizeof(keywords) / sizeof(Keyword);

    parse_table_.add(keywords, number_of_keywords);

    first_time = false;
  }

  //    sample = sample_sentinel_value;

  current_ = this;
}

//------------------------------------------------------------------------------------------------//
void Class_Parse_Table << class >> ::check_completeness(Token_Stream &tokens) {}

//------------------------------------------------------------------------------------------------//
SP << class >> Class_Parse_Table << class >> ::create_object() {
  return SP << class >> (new<class>());
}

//------------------------------------------------------------------------------------------------//
template <> SP << class >> parse_class << class >> (Token_Stream & tokens) {
  return parse_class_from_table<Class_Parse_Table << class> >> (tokens);
}

} // end namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// end of <pkg>/<class>__parser.hh
//------------------------------------------------------------------------------------------------//
