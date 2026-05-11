//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/runtime_check.cc
 * \author Kent Grimmett Budge
 * \brief  Member definitions of class runtime_check
 * \note   Copyright (C) 2018-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "diagnostics/runtime_check.hh"
#include "dsxx/Assert.hh"
#include "dsxx/DracoStrings.hh"
#include "c4/C4_Functions.hh"
#include <stdexcept>

namespace rtt_diagnostics {

//-------------------------------------------------------------------------//
void runtime_check(bool condition, char const *message) noexcept(false) {
  Require(message != nullptr);

  unsigned sum = !condition;
  rtt_c4::global_sum(sum);
  if (sum != 0) // some processors failed the condition
    throw std::runtime_error("runtime error detected on " + rtt_dsxx::to_string(sum) +
                             " processor(s): " + message);
}

} // end namespace rtt_diagnostics

//------------------------------------------------------------------------------------------------//
// end of diagnostics/runtime_check.cc
//------------------------------------------------------------------------------------------------//
