//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/path.cc
 * \brief  Encapsulate path information (path separator, etc.)
 * \note   Copyright (C) 2011-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "path.hh"

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/*! \brief Does the file exist?
 *
 * 1. Read the file attributes using stat.
 * 2. If we were able to get the file attributes so the file obviously exists.
 * 3. If we were not able to get the file attributes.  This may mean that we don't have permission
 *    to access the folder which contains this file. If you need to do that level of checking,
 *    lookup the return values of stat which will give you more details on why stat failed.
 *
 * \todo add attribute [[deprecated]] once C++17 is minimum supported language standard.
 */
bool fileExists(std::string const &strFilename) {
#if defined(__cplusplus) & __cplusplus > 201402L // C++17 or newer
  return std::filesystem::exists(strFilename);
#else // older than C++17
  return draco_getstat(strFilename).errorCode() == 0;
#endif
}

//------------------------------------------------------------------------------------------------//
/*! \brief Does the 'path' represent a directory?
*
* * https://en.cppreference.com/w/cpp/filesystem/is_directory
*
* \todo add attribute [[deprecated]] once C++17 is minimum supported language standard.
*/
bool isDirectory(std::string const &path) {
#if defined(__cplusplus) & __cplusplus > 201402L // C++17 or newer
  return std::filesystem::is_directory(path);
#else // older than C++17
  // If the path does not exist, then it cannot be a directory.
  if (!fileExists(path))
    return false;

  bool retVal(false);
  draco_getstat fileStatus(path.c_str());
  Check(fileStatus.valid());
  retVal = fileStatus.isdir();
  return retVal;
#endif
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
