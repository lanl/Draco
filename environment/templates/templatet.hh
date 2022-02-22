//-----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <pkg>/<class>.hh
 * \author <user>
 * \date   <date>
 * \brief  <start>
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

// clang-format off

#ifndef <spkg>_<class>_hh
#define <spkg>_<class>_hh

namespace<namespace> {

//================================================================================================//
/*!
 * \class <class>
 * \brief
 *
 * Long description or discussion goes here.  Information about Doxygen commands
 * can be found at http://www.doxygen.org.
 *
 * \sa <class>.cc for detailed descriptions.
 *
 * Code Sample:
 * \code
 *     cout << "Hello, world." << endl;
 * \endcode
 */
/*!
 * \example <pkg>/test/tst<class>.cc
 *
 * Test of <class>.
 */
//================================================================================================//

  template <typename T> class <class> {
  public:
    // NESTED CLASSES AND TYPEDEFS

    // CREATORS

    //! Default constructors.
    <class>();

    //! Copy constructor (the long doxygen description is in the .cc file).
    <class>(const<class><T> &rhs);

    //! Destructor.
    ~<class>();

    // MANIPULATORS

    //! Assignment operator for <class>.
    <class> &operator=(const<class><T> &rhs);

    // ACCESSORS

  private:
    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    // DATA
  };

} // end namespace <namespace>

#endif // <spkg>_<class>_hh

//------------------------------------------------------------------------------------------------//
// end of <pkg>/<class>.hh
//------------------------------------------------------------------------------------------------//
