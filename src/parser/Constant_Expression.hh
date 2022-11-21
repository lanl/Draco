//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/Constant_Expression.hh
 * \author Kent Budge
 * \date   Wed Jul 26 07:53:32 2006
 * \brief  Definition of class Constant_Expression
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef parser_Constant_Expression_hh
#define parser_Constant_Expression_hh

#include "utilities.hh"

namespace rtt_parser {
//------------------------------------------------------------------------------------------------//
void write_c(Unit const &units, ostream &out);

//================================================================================================//
/*!
 * \class Constant_Expression
 * \brief Represents a constant.
 *
 * This is the only concrete type of Expression that is currently available to clients. The need has
 * never arisen to make any other concrete type of Expression directly available.
 */
//================================================================================================//
class Constant_Expression : public Expression {
public:
  // CREATORS

  /*!
   * \brief Create a dimensionless constant.
   *
   * \param[in] number_of_variables Number of indepedent variables in the expression.
   * \param[in] value Value of the constant
   */
  Constant_Expression(unsigned const number_of_variables, double const value)
      : Expression(number_of_variables, value * dimensionless) {}

  /*!
   * \brief Create a dimensioned constant.
   *
   * \param[in] number_of_variables Number of indepedent variables in the expression.
   * \param[in] value Dimensions and value of the constant. The value is stored in the \c conv
   *                  member of this argument.
   */
  Constant_Expression(unsigned const number_of_variables, Unit const &value)
      : Expression(number_of_variables, value) {}

  // ACCESSORS

  bool is_constant() const override { return true; }

private:
  // IMPLEMENTATION

  double evaluate_(double const *const /*x*/) const override { return units().conv; }

  void write_(Precedence const /*precedence*/, vector<string> const & /*vars*/,
              ostream &out) const override {
    if (is_compatible(units(), dimensionless)) {
      out << units().conv;
    } else {
      write_c(units(), out);
    }
  }

  /*virtual*/ bool is_constant_(unsigned /*unused*/) const override { return true; }
};

} // end namespace rtt_parser

#endif // parser_Constant_Expression_hh

//------------------------------------------------------------------------------------------------//
// end of parser/Constant_Expression.hh
//------------------------------------------------------------------------------------------------//
