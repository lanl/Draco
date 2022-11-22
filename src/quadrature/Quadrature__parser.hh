//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Quadrature__parser.hh
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  Parser for various quadrature classes.
 * \note   Copyright (C) 2015-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Quadrature.hh"
#include "parser/Class_Parse_Table.hh"

namespace rtt_parser {
using rtt_quadrature::Quadrature;

//================================================================================================//
/*!
 * \class Quadrature_Parse_Table
 * \brief Parse table for Quadrature objects
 */
//================================================================================================//

template <> class Class_Parse_Table<Quadrature> {
public:
  // TYPEDEFS

  using Return_Class = Quadrature;

  // MANAGEMENT

  Class_Parse_Table();

  virtual ~Class_Parse_Table() { child_.reset(); }
  Class_Parse_Table(Class_Parse_Table const &rhs) = delete;
  Class_Parse_Table(Class_Parse_Table &&rhs) noexcept = delete;
  Class_Parse_Table &operator=(Class_Parse_Table const &rhs) = delete;
  Class_Parse_Table &operator=(Class_Parse_Table &&rhs) noexcept = delete;

  // SERVICES

  Parse_Table const &parse_table() const { return parse_table_; }

  bool allow_exit() const { return true; }

  void check_completeness(Token_Stream &tokens);

  std::shared_ptr<Quadrature> create_object();

  // STATICS

  static void register_quadrature(string const &keyword,
                                  std::shared_ptr<Quadrature> parse_function(Token_Stream &));

private:
  // STATICS

  static Parse_Table &get_parse_table() { return parse_table_; }
  static std::shared_ptr<Quadrature> &get_parsed_object();
  static Class_Parse_Table *current_;
  static Parse_Table parse_table_;
  static std::shared_ptr<Quadrature> child_;
};

} // namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// end of quadrature/Quadrature__parser.hh
//------------------------------------------------------------------------------------------------//
