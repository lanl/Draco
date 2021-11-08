//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_eospac/EospacException.hh
 * \author Kelly Thompson
 * \date   Fri Apr  6 13:59:06 2001
 * \brief  Header file for the cdi_eospac exception handler class.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_eospac_EospacException_hh
#define rtt_cdi_eospac_EospacException_hh

#include <stdexcept>
#include <string>

namespace rtt_cdi_eospac {

//================================================================================================//
/*!
 * \class EospacException
 *
 * \brief This class handles exceptions thrown by Eospac when it calls the EOSPAC library functions.
 *
 * This class provides an Eospac exception data object derived from std::exception.  When an
 * exception is thrown in the cdi_eospac package an EospacException object is created and may be
 * trapped by the calling routine by using a try/catch block.
 *
 * \sa cdi_eospac/test/tEospac.cc
 *
 * This example demonstrates how the try/catch blocks may be used to trap errors thrown by the
 * cdi_eospac package.  They also demonstrate how the calling program can extract information about
 * the exception.
 */
//================================================================================================//
class EospacException : public std::logic_error {

public:
  // CREATORS

  /*!
   * \brief The standard EospacException constructor.
   *
   * When an error is thrown, the data member "message" is set and then std::exception is
   * instantiated.
   *
   * The noexcept at the end of these declarations are exception specifiers. I have not limited the
   * types of exceptions that may be caught by this object (hence the empty parentheses). See
   * Lippman and Lajoie, "C++ Primer" p. 1042ff for more information.
   *
   * \param[in] msg A simple description of the exception. In some cases this may simply be the name
   *        of the Gandolf function that failed.  In other cases it might contain a detailed
   *        description of the error.  This string is used in the construction of std::exception.
   */
  explicit EospacException(std::string const &msg) noexcept : std::logic_error(msg) { /* empty */
  }

  /*!
   * \brief The standard EospacException destructor.
   *
   * \sa what()
   *
   * The destructor frees memory that may have been allocated by the accessor what().  If what was
   * never called then "message" was never allocated and no deallocation is required.
   *
   * Do not allow destructor to throw.
   */
  // ~EospacException() noexcept override = default;
};

//------------------------------------------------------------------------------------------------//
class EospacUnknownDataType : public EospacException {
public:
  explicit EospacUnknownDataType(std::string const &msg) noexcept
      : EospacException(msg) { /* empty */
  }
  // ~EospacUnknownDataType() noexcept override = default;
};

} // end namespace rtt_cdi_eospac

#endif // rtt_cdi_eospac_EospacException_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_eospac/EospacException.hh
//------------------------------------------------------------------------------------------------//
