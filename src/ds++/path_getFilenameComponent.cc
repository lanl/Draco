//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/path_getFilenameComponent.cc
 * \brief  Encapsulate path information (path separator, etc.)
 * \note   Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "path.hh"
#include <algorithm> // std::replace()
#include <sstream>

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Helper function to extract path information from a filename.
 *
 * \param[in] fqName A fully qualified filename (/path/to/the/unit/test)
 * \param fc Enum type FilenameComponent that specifies the action.
 * \return filename only, or path to file only.
 *
 * This function expects a fully qualified name of a unit test (e.g.: argv[0]).  It strips off the
 * path and returns the name of the unit test.
 *
 * Options:
 *    FC_PATH        Return path portion only for fqName
 *    FC_ABSOLUTE    not implemented.
 *    FC_NAME        Return filename w/o path for fqName
 *    FC_EXT         not implemented.
 *    FC_NAME_WE     not implemented.
 *    FC_REALPATH    resolve all symlinks
 *    FC_NATIVE      convert path to use native slashes for the current filesystem
 *    FC_LASTVALUE
 */
std::string getFilenameComponent(std::string const &fqName, FilenameComponent fc) {
  using std::string;
  string retVal;
  string::size_type idx;
  std::string fullName(fqName);

  switch (fc) {
  case FC_PATH: {
    // if fqName is a directory and ends with "/", trim the trailing dirSep.
    if (fqName.rfind(rtt_dsxx::UnixDirSep) == fqName.length() - 1)
      fullName = fqName.substr(0, fqName.length() - 1);
    if (fqName.rfind(rtt_dsxx::WinDirSep) == fqName.length() - 1)
      fullName = fqName.substr(0, fqName.length() - 1);

    idx = fullName.rfind(rtt_dsxx::UnixDirSep);
    if (idx == string::npos) {
      // Didn't find directory separator, as 2nd chance look for Windows directory separator.
      idx = fullName.rfind(rtt_dsxx::WinDirSep);
    }
    // If we still cannot find a path separator, return "./"
    if (idx == string::npos)
      retVal = string(".") + rtt_dsxx::dirSep;
    else
      retVal = fullName.substr(0, idx + 1);
    break;
  }

  case FC_NAME: {
    // if fqName is a directory and ends with "/", trim the trailing dirSep.
    if (fqName.rfind(rtt_dsxx::UnixDirSep) == fqName.length() - 1)
      fullName = fqName.substr(0, fqName.length() - 1);
    if (fqName.rfind(rtt_dsxx::WinDirSep) == fqName.length() - 1)
      fullName = fqName.substr(0, fqName.length() - 1);

    idx = fullName.rfind(UnixDirSep);
    if (idx == string::npos) {
      // Didn't find directory separator, as 2nd chance look for Windows directory separator.
      idx = fullName.rfind(WinDirSep);
    }
    // If we still cannot find a path separator, return the whole string as the test name.
    if (idx == string::npos)
      retVal = fullName;
    else
      retVal = fullName.substr(idx + 1);
    break;
  }

  case FC_REALPATH: {
    std::string path(getFilenameComponent(fqName, FC_PATH));
    std::string name(getFilenameComponent(fqName, FC_NAME));
    if (!draco_getstat(path).valid()) {
      retVal = std::string(); // On error, return empty string.
    } else {
      retVal = draco_getrealpath(path) + name;
    }
    break;
  }

  case FC_ABSOLUTE: {
    // Same as FC_REALPATH?
    std::string path(getFilenameComponent(fqName, FC_PATH));
    std::string name(getFilenameComponent(fqName, FC_NAME));
    if (!draco_getstat(path).valid()) {
      retVal = std::string(); // On error, return empty string.
    } else {
      retVal = draco_getrealpath(path) + name;
    }
    break;
  }

  case FC_EXT: {
    Insist(false, "case for FC_EXT not implemented.");
    break;
  }

  case FC_NAME_WE: {
    Insist(false, "case for FC_NAME_WE not implemented.");
    break;
  }

  case FC_NATIVE: {
    // This is always done before returning (see implementation found after the case statement)
    retVal = fullName;
    break;
  }

  default: {
    std::ostringstream msg;
    msg << "Unknown mode for rtt_dsxx::setName(). fc = " << fc;
    Insist(false, msg.str());
  }
  }

  // Always convert paths to use native format.
#ifdef _MSC_VER
  std::replace(retVal.begin(), retVal.end(), UnixDirSep, dirSep);
#else
  std::replace(retVal.begin(), retVal.end(), WinDirSep, dirSep);
#endif
  return retVal;
}

} // end namespace rtt_dsxx

//------------------------------------------------------------------------------------------------//
// end of path_getFilenameComponent.cc
//------------------------------------------------------------------------------------------------//
