//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/path.hh
 * \brief  Encapsulate path information (path separator, etc.)
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 *
 * \bug Consider replacing path.cc and path.hh with C++17 filesystem */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_path_hh
#define rtt_dsxx_path_hh

#include "Assert.hh"
#include "SystemCall.hh"
#include <iostream>
#if defined UNIX || defined MINGW
#include <dirent.h>   // struct DIR
#include <sys/stat.h> // struct stat; S_ISDIR
#endif

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
  void operator()(std::string const &dirpath) const { std::cout << dirpath << std::endl; }
};

//------------------------------------------------------------------------------------------------//
//! Functor for removing one item in a directory tree
class wdtOpRemove {
public:
  void operator()(std::string const &dirpath) const {
    std::cout << "Removing \"" << dirpath << "\"" << std::endl;
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
 *
 * Consider using Boost.FileSystem
 * \c boost::filesystem::remove_all(path);
 *
 * \code
 #include "boost/filesystem.hpp"
 #include <iostream>
 using namespace boost::filesystem;
 int main()
 {
 path current_dir("."); //
 for (recursive_directory_iterator iter(current_dir), end; iter != end; ++iter)
 std::cout << iter->path() << "\n";
 return 0;
 }
 * \endcode
 */
template <typename T>
void draco_walk_directory_tree(std::string const &dirname, T const &myOperator) {
  // If file does not exist, report and continue.
  if (!fileExists(dirname)) {
    std::cout << "File/directory \"" << dirname << "\"does not exist.  Continuing..." << std::endl;
    return;
  }

#ifdef _MSC_VER
  /*! \note If path contains the location of a directory, it cannot contain a trailing backslash. If
   * it does, -1 will be returned and errno will be set to ENOENT. */
  std::string d_name;
  if (dirname[dirname.size() - 1] == rtt_dsxx::WinDirSep ||
      dirname[dirname.size() - 1] == rtt_dsxx::UnixDirSep)
    d_name = dirname.substr(0, dirname.size() - 1);
  else
    d_name = dirname;

  // If this is not a directory, no recursion is needed:
  if (isDirectory(dirname)) {
    // Handle to the file/directory
    WIN32_FIND_DATA FileInformation;

    // Pattern to match all items in the current directory.
    std::string strPattern = d_name + "\\*.*";

    // Handle to directory
    HANDLE hFile = ::FindFirstFile(strPattern.c_str(), &FileInformation);

    // sanity check
    Insist(hFile != INVALID_HANDLE_VALUE, "Invalid file handle.");

    // Loop over all files in the current directory.
    do {
      // Do not process '.' or '..'
      if (FileInformation.cFileName[0] == '.')
        continue;

      std::string itemPath(d_name + "\\" + FileInformation.cFileName);

      // if the entry is a directory, recursively delete it, otherwise, delete the file:
      if (draco_getstat(itemPath).isdir())
        draco_walk_directory_tree(itemPath, myOperator);
      else
        myOperator(itemPath);

    } while (::FindNextFile(hFile, &FileInformation) == TRUE);

    // Close handle
    ::FindClose(hFile);
  }

  // Perform action on the top level entry

  myOperator(d_name);

#else
  // If this is not a directory, no recursion is needed:
  if (isDirectory(dirname)) {
    DIR *dir; // Handle to directory
    struct dirent *entry;
    // struct stat statbuf;

    dir = opendir(dirname.c_str());
    Insist(dir != nullptr, "Error opendir()");

    // Loop over all entries in the current directory.
    while ((entry = readdir(dir)) != nullptr) {
      std::string d_name(entry->d_name); // NOLINT [hicpp-no-array-decay]

      // Don't include "." or ".." entries.
      if (d_name[0] == '.')
        continue;

      std::string itemPath;
      std::string const uds(1, UnixDirSep);
      if (dirname[dirname.length() - 1] == UnixDirSep)
        itemPath.append(dirname).append(d_name);
      else
        itemPath.append(dirname).append(uds).append(d_name);

      // if the entry is a directory, recursively delete it, otherwise, delete the file
      if (draco_getstat(itemPath).isdir())
        draco_walk_directory_tree(itemPath, myOperator);
      else
        myOperator(itemPath);
    }
    closedir(dir);
  }

  // Perform action on the top level entry
  myOperator(dirname);

#endif

  return;
}

//------------------------------------------------------------------------------------------------//
//! Recursively remove a directory.
void draco_remove_dir(std::string const &path);
//! Recursively print a directory tree.
void draco_dir_print(std::string const &path);

} // end namespace rtt_dsxx

#endif // rtt_dsxx_path_hh

//------------------------------------------------------------------------------------------------//
// end of path.hh
//------------------------------------------------------------------------------------------------//
