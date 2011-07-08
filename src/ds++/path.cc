//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/path.cc
 * \author Kelly Thompson
 * \date   Thu Jul 08 10:21:36 2011
 * \brief  Encapsulate path information (path separator, etc.)
 * \note   Copyright � 2011 Los Alamos National Security, LLC
 *         All rights reserved.
 * \version $Id$
 */
//---------------------------------------------------------------------------//

#include "path.hh"
#include "Assert.hh"

#include <cerrno>       // errno
#include <cstring>      // strerror
#include <cstdlib>      // realpath
#include <sstream>

#include <sys/param.h>  // MAXPATHLEN
#include <sys/stat.h>   // stat
#include <unistd.h>

namespace rtt_dsxx
{

std::string currentPath(void)
{
    // Identify the current working directory.
    char curr_path[MAXPATHLEN];
    Insist(getcwd(curr_path, MAXPATHLEN) != NULL,
           std::string("getcwd failed: " + std::string(strerror(errno))));
    return std::string( curr_path );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Helper function to extract path information from a filename.
 * \param fqName A fully qualified filename (/path/to/the/unit/test)
 * \return filename only, or path to file only.
 *
 * This function expects a fully qualfied name of a unit test (e.g.:
 * argv[0]).  It strips off the path and returns the name of the unit test.
 *
 * Options:
 *    FC_PATH        Return path portion only for fqName
 *    FC_ABSOLUTE    not implemented.
 *    FC_NAME        Return filename w/o patyh for fqName
 *    FC_EXT         not implemented.
 *    FC_NAME_WE     not implemented.
 *    FC_REALPATH    resolve all symlinks
 *    FC_LASTVALUE
 */
std::string getFilenameComponent( std::string const & fqName,
                                  FilenameComponent fc )
{
    using std::string;
    string retVal;
    string::size_type idx;
    
    switch( fc )
    {
        case FC_PATH :
            idx=fqName.rfind( rtt_dsxx::UnixDirSep );
            if( idx == string::npos ) 
            {
                // Didn't find directory separator, as 2nd chance look for Windows
                // directory separator. 
                idx=fqName.rfind( rtt_dsxx::WinDirSep );
            }
            // If we still cannot find a path separator, return "./"
            if( idx == string::npos )
                retVal = string( string(".") + rtt_dsxx::dirSep );
            else
                retVal = fqName.substr(0,idx+1); 
            break;
            
        case FC_NAME :
            idx=fqName.rfind( UnixDirSep );
            if( idx == string::npos )
            {
                // Didn't find directory separator, as 2nd chance look for Windows
                // directory separator.
                idx=fqName.rfind( WinDirSep );
            }
            // If we still cannot find a path separator, return the whole
            // string as the test name.
            if( idx == string::npos )
                retVal = fqName;
            else
                retVal = fqName.substr(idx+1);
            break;

        case FC_REALPATH :
            struct stat buf;
            if( stat(fqName.c_str(), &buf) )
            {
                retVal = std::string();
            }
            else
            {
                char npath[PATH_MAX]; npath[0] = '\0';
                Insist(realpath(fqName.c_str(), npath) != NULL,
                       string("realpath failed on " + fqName));
                
                retVal = std::string(npath);
            }
            break;

        case FC_ABSOLUTE :
            Insist( false, "case for FC_ABSOLUTE not implemented." );
            break;
        case FC_EXT :
            Insist( false, "case for FC_EXT not implemented." );
            break;
        case FC_NAME_WE :
            Insist( false, "case for FC_NAME_WE not implemented." );
            break;

        default:
            std::ostringstream msg;
            msg << "Unknown mode for rtt_dsxx::setName(). fc = " << fc;
            Insist( false, msg.str() );
    }
    return retVal;
}





} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
// end of path.hh
//---------------------------------------------------------------------------//
