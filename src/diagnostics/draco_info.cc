//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/draco_info.cc
 * \author Kelly Thompson
 * \date   Wednesday, Nov 07, 2012, 18:49 pm
 * \brief  Small executable that prints the version and copyright strings.
 * \note   Copyright (C) 2012-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "draco_info.hh"
#include "c4/config.h"
#include "diagnostics/config.h"
#include "ds++/DracoStrings.hh"
#include "ds++/Release.hh"
#include "ds++/UnitTest.hh"
#include <algorithm> // tolower
#include <iostream>
#include <iterator>
#include <sstream>

namespace rtt_diagnostics {

//------------------------------------------------------------------------------------------------//
//! Constructor
DracoInfo::DracoInfo()
    : release(rtt_dsxx::release()), copyright(rtt_dsxx::copyright()),
      contact("For information, send e-mail to draco@lanl.gov."),
      build_type(rtt_dsxx::string_toupper(CBT)), library_type("static"), system_type("Unknown"),
      site_name("Unknown"), mpirun_cmd(""), diagnostics_level("disabled"), cxx(CMAKE_CXX_COMPILER),
      cxx_flags(CMAKE_CXX_FLAGS), cc(CMAKE_C_COMPILER), cc_flags(CMAKE_C_FLAGS), fc("none"),
      fc_flags("none"), cuda_compiler("none"), cuda_flags("none") {
#ifdef DRACO_SHARED_LIBS
  library_type = "Shared";
#endif
#ifdef CMAKE_SYSTEM_NAME
  system_type = CMAKE_SYSTEM_NAME_STRING;
#endif
#ifdef SITENAME
  site_name = SITENAME;
#endif
#ifdef HAVE_CUDA
  cuda = true;
#endif
#ifdef C4_MPI
  mpi = true;
  mpirun_cmd = std::string(MPIEXEC_EXECUTABLE) + std::string(" ") +
               std::string(MPIEXEC_NUMPROC_FLAG) + std::string(" <N> ");
#ifdef MPIEXEC_PREFLAGS
  mpirun_cmd += std::string(MPIEXEC_PREFLAGS);
#endif
#endif
#ifdef OPENMP_FOUND
  openmp = true;
#endif
#ifdef DRACO_DIAGNOSTICS
  std::ostringstream msg;
  msg << DRACO_DIAGNOSTICS;
  diagnostics_level = msg.str();
#endif
#ifdef DRACO_TIMING
  diagnostics_timing = true;
#endif
  if (build_type == std::string("RELEASE")) {
    cxx_flags += CMAKE_CXX_FLAGS_RELEASE;
    cc_flags += CMAKE_C_FLAGS_RELEASE;
  } else if (build_type == std::string("RELWITHDEBINFO")) {
    cxx_flags += CMAKE_CXX_FLAGS_RELWITHDEBINFO;
    cc_flags += CMAKE_C_FLAGS_RELWITHDEBINFO;
  } else if (build_type == std::string("DEBUG")) {
    cxx_flags += CMAKE_CXX_FLAGS_DEBUG;
    cc_flags += CMAKE_C_FLAGS_DEBUG;
  }
#ifdef CMAKE_Fortran_COMPILER
  fc = CMAKE_Fortran_COMPILER;
  fc_flags = CMAKE_Fortran_FLAGS;
  if (build_type == std::string("RELEASE"))
    fc_flags += CMAKE_Fortran_FLAGS_RELEASE;
  else if (build_type == std::string("DEBUG"))
    fc_flags += CMAKE_Fortran_FLAGS_DEBUG;
#endif
#ifdef CMAKE_CUDA_COMPILER
  cuda_compiler = CMAKE_CUDA_COMPILER;
  cuda_flags = CMAKE_CUDA_FLAGS;
  if (build_type == std::string("RELEASE"))
    cuda_flags += CMAKE_CUDA_FLAGS_RELEASE;
  else if (build_type == std::string("DEBUG"))
    cuda_flags += CMAKE_CUDA_FLAGS_DEBUG;
#endif
}

//------------------------------------------------------------------------------------------------//
void print_text_with_word_wrap(std::string const &longstring, size_t const indent_column,
                               size_t const max_width, std::ostringstream &msg,
                               std::string const &delimiters = " ") {
  // Preserve leading slash
  if (longstring[0] == '/')
    msg << "/";
  std::vector<std::string> const tokens = rtt_dsxx::tokenize(longstring, delimiters);
  std::string const delimiter(delimiters.substr(0, 1));
  size_t i(indent_column);
  for (auto item : tokens) {
    if (i + rtt_dsxx::remove_color(item).length() + 1 > max_width) {
      msg << "\n" << std::string(indent_column, ' ');
      i = indent_column;
    }
    msg << item;
    if (item != tokens.back())
      msg << delimiter;
    i += rtt_dsxx::remove_color(item).length() + 1;
  }
}

//------------------------------------------------------------------------------------------------//
std::string DracoInfo::fullReport() const {
  using std::cout;
  using std::endl;

  std::ostringstream infoMessage;

  // Create a list of features for DbC
  std::vector<std::string> dbc_info;
  dbc_info.emplace_back("Insist");
#ifdef REQUIRE_ON
  dbc_info.emplace_back("Require");
#endif
#ifdef CHECK_ON
  dbc_info.emplace_back("Check");
#endif
#ifdef ENSURE_ON
  dbc_info.emplace_back("Ensure");
#endif
#if DBC & 8
  dbc_info.emplace_back("no-throw version");
#endif
#if DBC & 16
  dbc_info.emplace_back("check-deferred version");
#endif

  // Print version and copyright information to the screen:
  infoMessage << briefReport();

  // Build Information
  //------------------
  size_t constexpr max_width(100);
  size_t const hanging_indent(std::string("    CXX Compiler      : ").length());

  infoMessage << "Build information:"
              << "\n    Build type        : " << build_type
              << "\n    Library type      : " << library_type
              << "\n    System type       : " << system_type
              << "\n    Site name         : " << site_name
              << "\n    CUDA support      : " << (cuda ? "enabled" : "disabled")
              << "\n    MPI support       : " << (mpi ? "enabled" : "disabled (c4 scalar mode)");

  if (mpi) {
    infoMessage << "\n      mpirun cmd      : ";
    print_text_with_word_wrap(mpirun_cmd, hanging_indent, max_width, infoMessage, "/");
  }

  infoMessage << "\n    OpenMP support    : " << (openmp ? "enabled" : "disabled")
              << "\n    Design-by-Contract: " << DBC << ", features = ";
  std::copy(dbc_info.begin(), dbc_info.end() - 1,
            std::ostream_iterator<std::string>(infoMessage, ", "));
  infoMessage << dbc_info.back();
  infoMessage << "\n    Diagnostics       : " << diagnostics_level
              << "\n    Diagnostics Timing: " << (diagnostics_timing ? "enabled" : "disabled");

  // Compilers and Flags
  infoMessage << "\n    CXX Compiler      : ";
  print_text_with_word_wrap(cxx, hanging_indent, max_width, infoMessage, "/");
  infoMessage << "\n    CXX_FLAGS         : ";
  print_text_with_word_wrap(cxx_flags, hanging_indent, max_width, infoMessage);
  infoMessage << "\n    C Compiler        : ";
  print_text_with_word_wrap(cc, hanging_indent, max_width, infoMessage, "/");
  infoMessage << "\n    C_FLAGS           : ";
  print_text_with_word_wrap(cc_flags, hanging_indent, max_width, infoMessage);
#ifdef CMAKE_Fortran_COMPILER
  infoMessage << "\n    Fortran Compiler  : ";
  print_text_with_word_wrap(fc, hanging_indent, max_width, infoMessage, "/");
  infoMessage << "\n    Fortran_FLAGS     : ";
  print_text_with_word_wrap(fc_flags, hanging_indent, max_width, infoMessage);
#endif
  if (cuda) {
    infoMessage << "\n    Cuda Compiler     : ";
    print_text_with_word_wrap(cuda_compiler, hanging_indent, max_width, infoMessage, "/");
    infoMessage << "\n    CUDA_FLAGS        : ";
    print_text_with_word_wrap(cuda_flags, hanging_indent, max_width, infoMessage);
  }

  infoMessage << "\n" << endl;

  return infoMessage.str();
}

//------------------------------------------------------------------------------------------------//
std::string DracoInfo::briefReport() const {
  std::ostringstream infoMessage;

  // Print version and copyright information to the screen:
  infoMessage << "\n";
  print_text_with_word_wrap(release, 5, 100, infoMessage, " ");
  infoMessage << "\n\n" << copyright << "\n" << contact << "\n" << std::endl;
  return infoMessage.str();
}

//------------------------------------------------------------------------------------------------//
//! extract the single-line version info from release and return it
std::string DracoInfo::versionReport() const {
  std::ostringstream infoMessage;
  print_text_with_word_wrap(release, 5, 100, infoMessage, ";");
  infoMessage << "\n" << std::endl;
  return infoMessage.str();
}

} // end namespace rtt_diagnostics

//------------------------------------------------------------------------------------------------//
// end of draco_info.cc
//------------------------------------------------------------------------------------------------//
