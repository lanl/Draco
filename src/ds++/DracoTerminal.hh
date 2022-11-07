//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/DracoTerminal.hh
 * \author Kelly Thompson
 * \date   Sat Oct 05 2019
 * \brief  Wrapper for a Terminal class that provides colored output.
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved.
 *
 * \sa https://github.com/certik/terminal/blob/master/terminal.h
 *
 * \todo Consider an enum class for colors that derives from $LS_COLORS on Linux.  This would allow
 *       color selection based on users's terminal colors (e.g.: light vs dark scheme). */
//------------------------------------------------------------------------------------------------//

#ifndef dsxx_dracoterminal_hh
#define dsxx_dracoterminal_hh

#include "terminal/terminal.h"
#include "ds++/config.h"
#include <array>

namespace Term {

//================================================================================================//
/*!
 * \class DracoTerminal
 * \brief Global scope singleton object to ensure terminal setup/teardown is done correctly.
 *
 * This will self construct on first call to Term::ccolor and it will remain in scope until the
 * program exits.  This design is based on the Meyer's Singleton Pattern, but the creation of the
 * instance MUST occur in the .cc file for proper exposure from MSVC DLLs. Some useful discussion
 * about this failure mode is at [this github project>(https://github.com/AlexWorx/ALib-Singletons).
 * The Meyer's singleton pattern has the advantage of avoiding heap allocation and possible memory
 * leaks.
 */
class DracoTerminal {

  // >> DATA <<

  /*! Construct a terminal object on creation.  This is will do some terminal initialization that is
   *  required for some older software (like Windows cmd.exe). */
  Term::Terminal term;

  //! Private constructor and destructor so that no objects can be created or destroyed.
  DracoTerminal() = default;
  ~DracoTerminal() = default;

  // >> DATA <<

  //! Set to false to disable all color.
  bool useColor = true;

public:
  //! Disable copy/move construction and assignment
  DracoTerminal(const DracoTerminal &) = delete;
  DracoTerminal(const DracoTerminal &&) = delete;
  DracoTerminal &operator=(const DracoTerminal &) = delete;
  DracoTerminal &operator=(const DracoTerminal &&) = delete;

  // >> DATA <<

  //! Colors based on condition (error, warn, etc).
  DLL_PUBLIC_dsxx static const std::array<uint32_t, 2> error;
  DLL_PUBLIC_dsxx static const std::array<uint32_t, 2> warning;
  DLL_PUBLIC_dsxx static const std::array<uint32_t, 2> note;
  DLL_PUBLIC_dsxx static const std::array<uint32_t, 1> quote;
  DLL_PUBLIC_dsxx static const std::array<uint32_t, 1> pass;
  DLL_PUBLIC_dsxx static const std::array<uint32_t, 2> fail;
  DLL_PUBLIC_dsxx static const std::array<uint32_t, 2> reset;

  // >> MANIPULATORS <<

  /*! \brief Calling this public function should always return true. It creates the singleton object
   *         on the first call. */
  static DracoTerminal &getInstance();

  //! Enable/Disable color text
  inline void enable_color_text() { useColor = true; }
  inline void disable_color_text() { useColor = false; }

  // >> ACCESSORS <<
  inline bool use_color() { return useColor; }
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Replacement for terminal/terminal.h's color function that will ensure the global
 *         singleton terminal object has been constructed prior to the use of color strings.
 *
 * \tparam T is an arbitrary type that can be converted to an integer that represents ansi-color
 *           codes.
 * \param[in] value represents a set of integer values that can be converted into ansi-color
 *               strings.  For example, "\033[01m" indicates bold text.  Use the codes provided in
 *               the DracoTerminal to make these more readable.
 */
template <typename T> std::string ccolor(T const value) {
  // Access the singleton wrapper for Term::Terminal
  Term::DracoTerminal &term_inst = Term::DracoTerminal::getInstance();
  std::string retVal;
  if (term_inst.use_color()) {
    retVal += "\033[" + std::to_string(static_cast<int>(value)) + "m";
  }
  return retVal;
}

//! specialization of ccolor to be used with Term::DracoTerminal::error, etc.
std::string ccolor(std::array<uint32_t, 1> const value);
//! specialization of ccolor to be used with Term::DracoTerminal::error, etc.
std::string ccolor(std::array<uint32_t, 2> const value);

} // namespace Term

#endif // dsxx_dracoterminal_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/DracoTerminal.hh
//------------------------------------------------------------------------------------------------//
