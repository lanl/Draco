//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Field_Traits.hh
 * \author Kent Budge
 * \brief  Define the Field_Traits class template.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef dsxx_Field_Traits_hh
#define dsxx_Field_Traits_hh

#include <complex>

namespace rtt_dsxx {

//================================================================================================//
/*!
 * \class Field_Traits
 * \brief Standardized description of field types
 *
 * Field_Traits is a traits class for field types.  Field types are types such as double or
 * std::complex<double> that represent a mathematical field. The Field_Traits class is useful in
 * template programming for capturing the properties of such field types.
 */
//================================================================================================//

template <typename Field> class Field_Traits {
public:
  //! Field types can be "labeled." For example, a value-plus-derivatives class has a field value
  //! that is labeled with its derivatives. The following typedef specifies the unlabeled type, by
  //! default the field type itself.
  using unlabeled_type = Field;

  //! Return the unique zero element of the field. By default, this is convertible from double 0.
  static Field zero() { return 0.0; }
  //! Return the unique unit element of the field. By default, this is convertible from double 1.
  static Field one() { return 1.0; }
};

//------------------------------------------------------------------------------------------------//
/*!
 * Sepcialization: const double
 *
 * This specialization eliminates a gcc warning about discarding the const qualifier on the return
 * type of the two member functions.
 */
template <> class Field_Traits<const double> {
public:
  //! Field types can be "labeled." For example, a value-plus-derivatives class has a field value
  //! that is labeled with its derivatives. The following typedef specifies the unlabeled type, by
  //! default the field type itself.
  using unlabeled_type = double;

  //! Return the unique zero element of the field. By default, this is convertible from double 0.
  static double zero() { return 0.0; }
  //! Return the unique unit element of the field. By default, this is convertible from double 1.
  static double one() { return 1.0; }
};

//------------------------------------------------------------------------------------------------//
// value.hh
//------------------------------------------------------------------------------------------------//
/*! Strip a field type of any labeling.
 *
 * Implicit conversion of a labeled field type (such as a value-plus-derivatives class) to an
 * underlying unlabeled field type can be dangerous. However, unlike conversion constructors,
 * conversion operators cannot be flagged explicit. Our alternative is to templatize a conversion
 * function from labeled field type to unlabeled field type.
 *
 * The default implementation assumes there is no labeling to strip.
 *
 * \arg \a Field A field type
 */
template <class Field> inline typename Field_Traits<Field>::unlabeled_type &value(Field &x) {
  return x;
}

//------------------------------------------------------------------------------------------------//
//! A version of the value function template for const arguments.
template <class Field>
inline typename Field_Traits<Field>::unlabeled_type const &value(Field const &x) {
  return x;
}

} // end namespace rtt_dsxx

#endif // traits_Field_Traits_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/Field_Traits.hh
//------------------------------------------------------------------------------------------------//
