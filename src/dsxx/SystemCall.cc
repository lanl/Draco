//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file  dsxx/SystemCall.cc
 * \brief Implementation for the Draco wrapper for system calls. This routine attempts to hide
 *        differences between Unix/Windows system calls.
 * \note  Copyright (C) 2024-2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "SystemCall.hh"
#include "Assert.hh"
#include <array>
#include <filesystem>
#include <string>

#if defined UNIX || defined MINGW
#include <cstring>    // strncpy()
#include <sys/stat.h> // stat
#include <unistd.h>   // gethostname
#endif
#ifdef _MSC_VER
#include <process.h> // _getpid
#endif

namespace rtt_dsxx {

//================================================================================================//
// FREE FUNCTIONS
//================================================================================================//

//------------------------------------------------------------------------------------------------//
/*! \brief Wrapper for system dependent hostname call.
 *
 * * Windows: \c HOST_NAME_MAX set to \c MAX_COMPUTERNAME_LENGTH in config.h
 * * Unix/Linux: \c HOST_NAME_MAX loaded from \<climits\>
 * * Mac OSX: \c HOST_NAME_MAX set to \c _POSIX_HOST_NAME_MAX in config.h
 *
 * \sa sec_HOST_NAME_MAX
 */
std::string draco_gethostname() {
// Windows: gethostname from <winsock2.h>
#ifdef _MSC_VER
  std::array<char, HOST_NAME_MAX> hostname;
  hostname.fill('x');
  int err = gethostname(&hostname[0], sizeof(hostname));
  if (err) {
    return "gethostname() failed!";
  }
  return hostname.data();

#else

// Linux: gethostname from <unistd.h>
#ifdef HAVE_GETHOSTNAME
  std::array<char, HOST_NAME_MAX> hostname{}; // NOLINT(misc-include-cleaner)
  hostname.fill('y');
  int const err = gethostname(&hostname[0], HOST_NAME_MAX); // NOLINT(misc-include-cleaner)
  if (err)
    strncpy(&hostname[0], "gethostname() failed", HOST_NAME_MAX); // NOLINT(misc-include-cleaner)
  return hostname.data();

// Catamount systems do not have gethostname().
#else
  return std::string("Host (unknown)");
#endif
#endif
} // draco_gethostname

//------------------------------------------------------------------------------------------------//
/*! \brief Wrapper for system dependent pid call..
 *
 * Catamount systems do not have getpid().  This function will return -1.
 */
int draco_getpid() {
#ifdef _MSC_VER
  int i = _getpid();
  return i;
#else

#ifdef HAVE_GETHOSTNAME
  return getpid();
#else
  // Catamount systems do not have getpid().  This function will return -1.
  return -1;
#endif
#endif
} // draco_pid

//------------------------------------------------------------------------------------------------//
/*! \brief Wrapper for system dependent getcwd call.  Identify the current working directory.
 *
 *  This should always return a trailing directory separator. */
std::string draco_getcwd() {
  return std::filesystem::current_path().generic_string() + rtt_dsxx::dirSep;
}

//------------------------------------------------------------------------------------------------//
/*! \brief Wrapper for system dependent stat call.
 *
 * Windows
 *    http://msdn.microsoft.com/en-us/library/14h5k7ff%28v=vs.80%29.aspx
 * Linux
 *    http://en.wikipedia.org/wiki/Stat_%28system_call%29
 *
 * \todo Once C++17 is the minimum language standard for Draco, mark this function as
 *       [[deprecated]].  Use C++17's filesystem functions instead of this cumbersome structure.
 */
#ifdef _MSC_VER
draco_getstat::draco_getstat(std::string const &fqName) {
  filefound = true;
  /*! \note If path contains the location of a directory, it cannot contain a trailing backslash. If
   * it does, -1 will be returned and errno will be set to ENOENT. */
  std::string clean_fqName;
  if (fqName[fqName.size() - 1] == rtt_dsxx::WinDirSep ||
      fqName[fqName.size() - 1] == rtt_dsxx::UnixDirSep)
    clean_fqName = fqName.substr(0, fqName.size() - 1);
  else
    clean_fqName = fqName;

  stat_return_code = _stat(clean_fqName.c_str(), &buf);
  if (stat_return_code != 0) {
    switch (errno) {
    case ENOENT: // file not found
      filefound = false;
      break;
    case EINVAL:
      Insist(stat_return_code == 0, "invalid parameter given to _stat.");
      break;
    default:
      /* should never get here */
      Insist(stat_return_code == 0, "_stat returned an error.");
    }
  }

  if (filefound) {
    // Handle to directory
    HANDLE hFile = ::FindFirstFile(clean_fqName.c_str(), &FileInformation);
    // sanity check
    Insist(hFile != INVALID_HANDLE_VALUE, "Invalid file handle.");
    // Close handle
    ::FindClose(hFile);
  }
}
#else
draco_getstat::draco_getstat(std::string const &fqName) {
  stat_return_code = stat(fqName.c_str(), &buf);
}
#endif

//------------------------------------------------------------------------------------------------//
//! Is this a regular file?
bool draco_getstat::isreg() {
#ifdef _MSC_VER
  bool b = FileInformation.dwFileAttributes & FILE_ATTRIBUTE_NORMAL;
  return filefound && b;
#else
  bool b = S_ISREG(buf.st_mode); // NOLINT [hicpp-signed-bitwise]
  return b;
#endif
}

//------------------------------------------------------------------------------------------------//
//! Is this a directory?
bool draco_getstat::isdir() {
#ifdef _MSC_VER
  bool b = FileInformation.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY;
  return filefound && b;
#else
  bool b = S_ISDIR(buf.st_mode); // NOLINT [hicpp-signed-bitwise]
  return b;
#endif
}

//------------------------------------------------------------------------------------------------//
//! Is a Unix permission bit set?
#ifdef _MSC_VER
bool draco_getstat::has_permission_bit(int /*mask*/) {
  Insist(isreg(), "Can only check permission bit for regular files.");
  Insist(false, "draco_getstat::has_permission_bit() not implemented for WIN32");
  return false;
}
#else
bool draco_getstat::has_permission_bit(int mask) {
  Insist(isreg(), "Can only check permission bit for regular files.");
  // check execute bit (buf.st_mode & 0111)
  return (buf.st_mode & mask); // NOLINT [hicpp-signed-bitwise]
}
#endif

//------------------------------------------------------------------------------------------------//
//! Wrapper for system dependent realpath call.
std::string draco_getrealpath(std::string const &path) {
  return std::filesystem::canonical(path).generic_string() + rtt_dsxx::dirSep;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Make a directory
 *
 * \sa  boost::filesystem::create_directories("/tmp/a/b/c");
 *
 * \todo Do we need a permissions argument?  For now, just create directories with 0700 permission.
 */
void draco_mkdir(std::string const &path) {
  namespace fs = std::filesystem;
  fs::create_directories(path);
  fs::permissions(path, fs::perms::owner_all, fs::perm_options::replace);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Delete a single file or directory (not recursive)
 *
 * \sa wdtOpRemove
 *
 * Sample implementation for Win32 (uses Win API which I don't want to do)
 * http://forums.codeguru.com/showthread.php?239271-Windows-SDK-File-System-How-to-delete-a-directory-and-subdirectories
 * http://stackoverflow.com/questions/1468774/recursive-directory-deletion-with-win32
 * http://msdn.microsoft.com/en-us/library/aa365488%28VS.85%29.aspx
 *
 * Sample implementation for Unix
 * http://www.linuxquestions.org/questions/programming-9/deleting-a-directory-using-c-in-linux-248696/
 */
void draco_remove(std::string const &dirpath) {

  std::filesystem::remove_all(dirpath);

  // If the file still exists this check will fail.
  Ensure(!draco_getstat(dirpath).valid());
  return;
}

} // end of namespace rtt_dsxx

//------------------------------------------------------------------------------------------------//
// end of SystemCall.cc
//------------------------------------------------------------------------------------------------//
