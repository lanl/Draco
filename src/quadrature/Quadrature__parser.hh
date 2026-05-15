//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Quadrature__parser.hh
 * \author Kelly Thompson
 * \brief  Parser for various quadrature classes.
 * \note   Copyright (C) 2015-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "parser/Abstract_Class_Parser.hh"
#include "parser/Class_Parse_Table.hh"
#include "quadrature/Quadrature.hh"

namespace rtt_parser {

//================================================================================================//
/*!
 * \class Quadrature_Parse_Table
 * \brief Parse table for Quadrature objects
 */
//================================================================================================//

template <> class Class_Parse_Table<rtt_quadrature::Quadrature> {
public:
  // TYPEDEFS

  using Return_Class = rtt_quadrature::Quadrature;

  // MANAGEMENT

  Class_Parse_Table();

  virtual ~Class_Parse_Table() { child_.reset(); }
  Class_Parse_Table(Class_Parse_Table const &rhs) = delete;
  Class_Parse_Table(Class_Parse_Table &&rhs) noexcept = delete;
  Class_Parse_Table &operator=(Class_Parse_Table const &rhs) = delete;
  Class_Parse_Table &operator=(Class_Parse_Table &&rhs) noexcept = delete;

  // SERVICES

  [[nodiscard]] Parse_Table const &parse_table() const { return parse_table_; }

  [[nodiscard]] bool allow_exit() const { return true; }

  void check_completeness(rtt_parser::Token_Stream &tokens);

  std::shared_ptr<rtt_quadrature::Quadrature> create_object();

  // STATICS

  /*!
   * \brief This function allows local developers to add new transport models to Capsaicin, by
   *        adding their names to the set of transport models recognized by the \c Transport_Model
   *        parser.
   *
   * \param keyword Name of the new transport model. Must be unique.
   * \param parse_function Function to call to parse a specification for the new transport model.
   *        Must return an instance of the transport model class.
   */
  static void register_quadrature(
      std::string const &keyword,
      std::shared_ptr<rtt_quadrature::Quadrature> parse_function(rtt_parser::Token_Stream &)) {
    Abstract_Class_Parser<rtt_quadrature::Quadrature, get_parse_table,
                          get_parsed_object>::register_child(keyword, parse_function);
  }

private:
  // STATICS

  static Parse_Table &get_parse_table() { return parse_table_; }
  static std::shared_ptr<rtt_quadrature::Quadrature> &get_parsed_object();
  static Class_Parse_Table *current_;
  static Parse_Table parse_table_;
  static std::shared_ptr<rtt_quadrature::Quadrature> child_;
};

} // namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// end of quadrature/Quadrature__parser.hh
//------------------------------------------------------------------------------------------------//
