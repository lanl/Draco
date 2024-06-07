//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/Release.hh
 * \author Thomas Evans
 * \date   Thu Jul 15 09:31:44 1999
 * \brief  Header file for dsxx library release function.
 * \note   Copyright (C) 2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_ds_Release_hh
#define rtt_ds_Release_hh

#include "dsxx/config.h"
#include <algorithm>
#include <functional>
#include <map>
#include <string>

namespace rtt_dsxx {

// Typedefs
using mmdevs = std::multimap<int, std::string, std::greater<>>;
using fomdev = std::pair<int, std::string>;

//! Query package for the release number.
const std::string release();
//! Return a list of Draco authors
const std::string author_list(bool const use_doxygen_formatting = false);
//! Return a list of Draco authors
const std::string copyright();

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Format list of authors to do correct line breaks and ensures total line length is less
 *        than a specified maximum.
 *
 * \arg[in] maxlen Maximum line length
 * \arg[in] line_name Category title
 * \arg[in] list of developers
 * \return A formatted message.
 */
std::string print_devs(size_t const maxlinelen, std::string const &line_name, mmdevs const &devs);

} // namespace rtt_dsxx

//! This version can be called by Fortran and wraps the C++ version.
extern "C" void ec_release(char *release_string, size_t maxlen);

#endif // rtt_ds_Release_hh

//------------------------------------------------------------------------------------------------//
// end of dsxx/Release.hh
//------------------------------------------------------------------------------------------------//
