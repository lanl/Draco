//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file  ds++/SystemCall.cc
 * \brief Implementation for the Draco wrapper for system calls. This routine attempts to hide
 *        differences between Unix/Windows system calls.
 * \note  Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "SystemCall.hh"
#include "Assert.hh"

#include <cerrno>  // errno
#include <cstdio>  // remove()
#include <cstdlib> // _fullpath
#include <cstring> // strncpy()
#if defined UNIX || defined MINGW
#include <sys/param.h> // MAXPATHLEN
#include <unistd.h>    // gethostname
#endif
#ifdef _MSC_VER
#include <direct.h>  // _getcwd
#include <process.h> // _getpid
#include <sstream>
#endif
#include <array>
#include <iostream>

namespace rtt_dsxx {

//================================================================================================//
// FREE FUNCTIONS
//================================================================================================//

//------------------------------------------------------------------------------------------------//
/*! \brief Wrapper for system dependent hostname call.
 *
 * Windows:
 *     \c HOST_NAME_MAX set to \c MAX_COMPUTERNAME_LENGTH in config.h
 *
 * Catamount systems:
 *     \c HOST_NAME_MAX hard coded by CMake in config.h
 *
 * Unix/Linux:
 *     \c HOST_NAME_MAX loaded from \<climit\>
 *
 * Mac OSX:
 *     \c HOST_NAME_MAX set to \c _POSIX_HOST_NAME_MAX in config.h
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
  std::array<char, HOST_NAME_MAX> hostname{};
  hostname.fill('y');
  int err = gethostname(&hostname[0], HOST_NAME_MAX);
  if (err)
    strncpy(&hostname[0], "gethostname() failed", HOST_NAME_MAX);
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
/*! \brief Wrapper for system dependent getcwd call.
 *
 *  This should always return a trailing directory separator.
 */
std::string draco_getcwd() {
// Identify the current working directory.
#ifdef _MSC_VER
  char *buffer;
  Insist((buffer = _getcwd(nullptr, 0)) != nullptr,
         std::string("getcwd failed: " + std::string(strerror(errno))));
  std::string cwd(buffer, buffer + strnlen(buffer, MAXPATHLEN));
  free(buffer);
#else
  std::array<char, MAXPATHLEN> curr_path{};
  curr_path.fill('z');
  Insist(getcwd(&curr_path[0], MAXPATHLEN) != nullptr,
         "getcwd failed: " + std::string(strerror(errno)));
  std::string cwd(curr_path.data());
#endif

  // Ensure that the last character is a directory separator.
  if (cwd.at(cwd.length() - 1) != rtt_dsxx::dirSep)
    cwd = cwd + rtt_dsxx::dirSep;

  return cwd;
}

//------------------------------------------------------------------------------------------------//
/*! \brief Wrapper for system dependent stat call.
 *
 * Windows
 *    http://msdn.microsoft.com/en-us/library/14h5k7ff%28v=vs.80%29.aspx
 * Linux
 *    http://en.wikipedia.org/wiki/Stat_%28system_call%29
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
  Insist(false, "draco_getstat::hsa_permission_bit() not implemented for WIN32");
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
  std::array<char, MAXPATHLEN> buffer{}; // _MAX_PATH
  buffer.fill('a');
#ifdef _MSC_VER
  // http://msdn.microsoft.com/en-us/library/506720ff%28v=vs.100%29.aspx
  Insist(_fullpath(&buffer[0], path.c_str(), MAXPATHLEN) != nullptr, "Invalid path.");
  std::string retVal(buffer.data());
#else
  Insist((realpath(path.c_str(), &buffer[0])) != nullptr, "Invalid path.");
  std::string retVal(buffer.data());
  if (draco_getstat(retVal).isdir()) {
    // realpath trims the trailing slash, append now.
    retVal += std::string(&dirSep, 1);
  }
#endif
  return retVal;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Make a directory
 *
 * \sa  boost::filesystem::create_directories("/tmp/a/b/c");
 *
 * \todo Do we need a permissions argument?
 */
void draco_mkdir(std::string const &path) {
#ifdef _MSC_VER
  draco_getstat dirInfo(path);
  if (!dirInfo.isdir()) {
    /*! \note If path contains the location of a directory, it cannot contain a trailing
     * backslash. If it does, -1 will be returned and errno will be set to ENOENT. */
    std::string clean_fqName;
    if (path[path.size() - 1] == rtt_dsxx::WinDirSep ||
        path[path.size() - 1] == rtt_dsxx::UnixDirSep)
      clean_fqName = path.substr(0, path.size() - 1);
    else
      clean_fqName = path;

    // make the directory
    int errorCode = _mkdir(clean_fqName.c_str());
    if (errorCode == -1) {
      switch (errno) {
      case EEXIST:
        Insist(errno != EEXIST, "ERROR: Unable to create directory," + clean_fqName +
                                    ", because it already exists.");
        break;
      case ENOENT:
        Insist(errno != ENOENT, "ERROR: Unable to create directory," + clean_fqName +
                                    ", because the path is not found.");
        break;
      default:
        /* should never get here */
        Insist(errno == 0, "ERROR: Unable to create directory, " + clean_fqName);
      }
    }
  }
#else
  mkdir(path.c_str(), 0700);
// S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
#endif
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
 *
 * Consider using Boost.FileSystem
 */
void draco_remove(std::string const &dirpath) {
  // remove() works for all unix items but only for files (not directories) for windows.
#ifdef _MSC_VER
  if (draco_getstat(dirpath).isdir()) {
    // Clear any special directory attributes.
    bool ok = ::SetFileAttributes(dirpath.c_str(), FILE_ATTRIBUTE_NORMAL);
    if (!ok) {
      int myerr = ::GetLastError();
      std::ostringstream msg;
      msg << "ERROR: File attribute not normal. myerr = " << myerr << ", file = " << dirpath;
      Insist(ok, msg.str());
    }

    // Delete directory
    ok = ::RemoveDirectory(dirpath.c_str());
    if (!ok) {
      int myerr = ::GetLastError();
      std::ostringstream msg;
      msg << "ERROR: Error deteting file, myerr = " << myerr << ", file = " << dirpath;
      Insist(ok, msg.str());
    }
  } else {
    remove(dirpath.c_str());
  }
#else
  remove(dirpath.c_str());
#endif
  // If the file still exists this check will fail.
  Ensure(!draco_getstat(dirpath).valid());
  return;
}

} // end of namespace rtt_dsxx

//------------------------------------------------------------------------------------------------//
// end of SystemCall.cc
//------------------------------------------------------------------------------------------------//
