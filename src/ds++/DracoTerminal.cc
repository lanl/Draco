//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/DracoTerminal.cc
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

#include "DracoTerminal.hh"

//------------------------------------------------------------------------------------------------//
/*! \brief Initialize other static const data (color map)
 *
 * \todo Eval ENV{GCC_COLORS} to override the defaults listed below? This probably would need to be
 *       done by cmake and saved to ds++/config.h.
 *
 * This should provide a mechanism to distinguish between light and dark terminal backgrounds.
 *
 * \sa https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Message-Formatting-Options.html
 */
const std::array<uint32_t, 2> Term::DracoTerminal::error = {01, 31};
const std::array<uint32_t, 2> Term::DracoTerminal::warning = {01, 35};
const std::array<uint32_t, 2> Term::DracoTerminal::note = {01, 36};
const std::array<uint32_t, 1> Term::DracoTerminal::quote = {01};
const std::array<uint32_t, 1> Term::DracoTerminal::pass = {32};
const std::array<uint32_t, 2> Term::DracoTerminal::fail = {01, 31};
const std::array<uint32_t, 2> Term::DracoTerminal::reset = {00, 39};

//------------------------------------------------------------------------------------------------//
//! specialization of ccolor for T=array<uint32_t> to be us

//! Instantiate the singleton as a static object that exists for the lifetime of the program.
Term::DracoTerminal &Term::DracoTerminal::getInstance() {
  static Term::DracoTerminal instance;
  return instance;
}

//! specialization of ccolor to be used with Term::DracoTerminal::error, etc.
std::string Term::ccolor(std::array<uint32_t, 1> const value) {
  // Access the singleton wrapper for Term::Terminal
  Term::DracoTerminal &term_inst = Term::DracoTerminal::getInstance();
  if (term_inst.use_color())
    return "\033[" + std::to_string(static_cast<int>(value[0])) + "m";
  return "";
}
//! specialization of ccolor to be used with Term::DracoTerminal::error, etc.
std::string Term::ccolor(std::array<uint32_t, 2> const value) {
  // Access the singleton wrapper for Term::Terminal
  Term::DracoTerminal &term_inst = Term::DracoTerminal::getInstance();
  std::string retVal;
  if (term_inst.use_color())
    for (uint32_t const &it : value)
      retVal += "\033[" + std::to_string(static_cast<int>(it)) + "m";
  return retVal;
}

//------------------------------------------------------------------------------------------------//
// end of ds++/DracoTerminal.cc
//------------------------------------------------------------------------------------------------//
