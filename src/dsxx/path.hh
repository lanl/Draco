//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/path.hh
 * \brief  Encapsulate path information (path separator, etc.)
 * \note   Copyright (C) 2024-2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_path_hh
#define rtt_dsxx_path_hh

#include "Assert.hh"
#include "SystemCall.hh"
#include <filesystem>
#include <iostream>

namespace rtt_dsxx {

enum FilenameComponent {
  FC_PATH, //!< Extract path portion of fully qualified filename
  FC_ABSOLUTE,
  FC_NAME, //!< Extract filename portion (w/o path).
  FC_EXT,
  FC_NAME_WE,
  FC_REALPATH,
  FC_NATIVE, //!< Convert directory separator to native format.
  FC_LASTVALUE
};

//------------------------------------------------------------------------------------------------//
//! Get a specific component of a full filename.
std::string getFilenameComponent(std::string const &fqName, FilenameComponent fc);

//------------------------------------------------------------------------------------------------//
//! Does the file exist?
bool fileExists(std::string const &filename);
bool isDirectory(std::string const &path);

//------------------------------------------------------------------------------------------------//
//! Functor for printing all items in a directory tree
class wdtOpPrint {
public:
  void operator()(std::string const &dirpath) const { std::cout << dirpath << '\n'; }
};

//------------------------------------------------------------------------------------------------//
//! Functor for removing one item in a directory tree
class wdtOpRemove {
public:
  void operator()(std::string const &dirpath) const {
    std::cout << "Removing \"" << dirpath << "\"" << '\n';
    draco_remove(dirpath);
  }
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Walk a directory tree structure, perform myOperator() action on each entry.
 * \arg dirname String representing the top node of the directory to be parsed.
 * \arg myOperator Functor that defines action to be taken on each entry in the directory. Recommend
 *      using wdtOpPrint or wdtOpRemove
 *
 * \sa draco_remove_dir Helper function to recursively delete a directory and all its contents.
 * \sa draco_dir_print Helper function that will print a directory and all its contents.
 *
 * Sample implementation for Win32 (uses Win API which I don't want to do)
 * http://forums.codeguru.com/showthread.php?239271-Windows-SDK-File-System-How-to-delete-a-directory-and-subdirectories
 * http://stackoverflow.com/questions/1468774/recursive-directory-deletion-with-win32
 * http://msdn.microsoft.com/en-us/library/aa365488%28VS.85%29.aspx
 *
 * Sample implementation for Unix
 * http://www.linuxquestions.org/questions/programming-9/deleting-a-directory-using-c-in-linux-248696/
 */
template <typename T>
void draco_walk_directory_tree(std::string const &dirname, T const &myOperator) {

  namespace fs = std::filesystem;

  // If file does not exist, report and continue.
  if (!fs::exists(dirname)) {
    std::cout << "File/directory \"" << dirname << "\"does not exist.  Continuing...\n";
    return;
  }

  if (fs::is_directory(dirname)) {
    // Loop over all entries in the current directory.
    for (auto const &dir_entry : fs::directory_iterator(dirname)) {
      if (fs::is_directory(dir_entry.path()))
        draco_walk_directory_tree(dir_entry.path().generic_string(), myOperator);
      else
        myOperator(dir_entry.path().generic_string());
    }
  }
  // Perform action on the top level entry
  myOperator(dirname);

  return;
}

//------------------------------------------------------------------------------------------------//
//! Recursively remove a directory.
void draco_remove_dir(std::string const &path);
//! Recursively print a directory tree.
void draco_dir_print(std::string const &path);

} // end namespace rtt_dsxx

#endif

//------------------------------------------------------------------------------------------------//
// end of path.hh
//------------------------------------------------------------------------------------------------//
