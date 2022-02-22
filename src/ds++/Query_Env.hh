//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Query_Env.hh
 * \author Tim Kelley
 * \date   Fri Jun 7 08:06:53 2019
 * \brief  Functions for working with your environment
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef Query_Env_hh
#define Query_Env_hh

#include <cstdlib> // getenv, setenv
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility> // pair

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Get a value from the environment
 *
 * \tparam T the type of the value, must be default constructible
 *
 * \param key the key to which you would like the corresponding value
 * \param default_value a default value to return if key undefined.
 * \return: {whether key was defined, corresponding value}
 *
 * \note Calls POSIX getenv() function, which is not required to be re-entrant.  If the key was set
 *       in the environment, get_env_val() returns the value converted to type T, that was set in
 *       the environment. If the key was not set in the environment, it returns the default_value
 *       that the caller provides. The first argument of the pair describes whether the key was
 *       defined.
 *
 * \example ds++/test/tstQuery_Env.cc
 */
template <typename T>
std::pair<bool, T> get_env_val(std::string const &key, T default_value = T{}) {
  static_assert(std::is_default_constructible<T>::value,
                "get_env_val only usable with default constructible types");
  T val{default_value};
  char *s_val = getenv(key.c_str());
  bool is_defined(false);
  if (s_val) {
    std::stringstream val_str(s_val);
    val_str >> val;
    is_defined = true;
  }
  return std::make_pair(is_defined, val);

} // get_env_val

} // namespace rtt_dsxx

#endif // Query_Env

//------------------------------------------------------------------------------------------------//
// end ds++/Query_Env.hh
//------------------------------------------------------------------------------------------------//
