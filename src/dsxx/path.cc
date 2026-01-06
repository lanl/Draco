//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/path.cc
 * \brief  Encapsulate path information (path separator, etc.)
 * \note   Copyright (C) 2024-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "path.hh"
#include <filesystem>
#include <string>

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/*! \brief Does the file exist?
 *
 * 1. Read the file attributes using stat.
 * 2. If we were able to get the file attributes so the file obviously exists.
 * 3. If we were not able to get the file attributes.  This may mean that we don't have permission
 *    to access the folder which contains this file. If you need to do that level of checking,
 *    lookup the return values of stat which will give you more details on why stat failed.
 */
[[deprecated]] bool fileExists(std::string const &strFilename) {
  return std::filesystem::exists(strFilename);
}

//------------------------------------------------------------------------------------------------//
/*! \brief Does the 'path' represent a directory?
*
* * https://en.cppreference.com/w/cpp/filesystem/is_directory
*/
[[deprecated]] bool isDirectory(std::string const &path) {
  return std::filesystem::is_directory(path);
}

//------------------------------------------------------------------------------------------------//
//! Recursively remove a directory.
void draco_remove_dir(std::string const &path) {
  draco_walk_directory_tree(path, wdtOpRemove());
  return;
}

//------------------------------------------------------------------------------------------------//
//! Recursively print a directory tree.
void draco_dir_print(std::string const &path) {
  draco_walk_directory_tree(path, wdtOpPrint());
  return;
}

} // end namespace rtt_dsxx

//------------------------------------------------------------------------------------------------//
// end of path.hh
//------------------------------------------------------------------------------------------------//
