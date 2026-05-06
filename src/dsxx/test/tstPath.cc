//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/test/tstPath.cc
 * \author Kelly Thompson
 * \date   Tue Jul 12 16:00:59 2011
 * \brief  Test functions found in dsxx/path.hh and path.cc
 * \note   Copyright (C) 2024-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/SystemCall.hh"
#include "dsxx/UnitTest.hh"
#include "dsxx/path.hh"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace {

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void test_currentPath(ScalarUnitTest &ut) {
  cout << "\nTesting currentPath() function ... \n\n";

  // currentPath
  string const cp = draco_getcwd();

  // if we got here, currentPath didn't throw.

  // Note, we have no idea where the this test was run from so we can say nothing about what the
  // path string should contain.

  if (std::filesystem::exists(cp))
    PASSMSG(string("Retrieved current path exists. cp = ") + cp);
  else
    FAILMSG("Retrieved current path does not exist. cp = " + cp);

  // Test behavior of fileExist when file does not exist.
  string const fileDoesNotExist("/bounty/hunter/boba_fett");

  if (!std::filesystem::exists(fileDoesNotExist))
    PASSMSG(string("fileExist() correctly returned false. "));
  else
    FAILMSG(string("fileExist() incorrectly returned true. "));

  return;
}

//------------------------------------------------------------------------------------------------//
void test_getFilenameComponent(ScalarUnitTest &ut, string const &fqp) {
  // Convert path to Native format
  std::string const fqpn = getFilenameComponent(fqp, rtt_dsxx::FC_NATIVE);

  cout << "\nTesting getFilenameComponent() function with fqpn = " << fqpn << " ...\n\n";
  ;

  bool usesUnixDirSep = true;

  // test the FC_PATH mode
  // ------------------------------------------------------------

  // 4 possible cases: ./tstPath, ../test/tstPath, tstPath.exe or test/tstPath.exe

  // Does the provided path use unix or windows directory separator?
  string::size_type idx = fqpn.find(rtt_dsxx::UnixDirSep);
  if (idx == string::npos)
    usesUnixDirSep = false;

  // retrieve the path w/o the filename.
  string const mypath = getFilenameComponent(fqpn, rtt_dsxx::FC_PATH);
  if (usesUnixDirSep) {
    // If we are using UnixDirSep, we have 2 cases (./tstPath or ../test/tstPath).  Look for the
    // case with 'test' first:
    idx = mypath.find(string("test") + rtt_dsxx::UnixDirSep);

    // If the return string does not have 'test/' then we also need to check for './' as a pass
    // condition
    if (idx == string::npos)
      idx = mypath.find(rtt_dsxx::UnixDirSep);

    // Report pass/fail
    if (idx != string::npos)
      PASSMSG(string("Found expected partial path. Path = ") + mypath);
    else
      FAILMSG("Did not find expected partial path. Expected path = " + mypath);
  } else {
    // If we are using WinDirSep, we have 2 cases (.\tstPath.exe or ...\test\tstPath.exe).  Look for
    // the case with 'test' first:
    idx = mypath.find(string("test") + rtt_dsxx::WinDirSep);

    // If the return string does not have 'test\' then we also need to check for './' as a pass
    // condition
    if (idx == string::npos)
      idx = mypath.find(rtt_dsxx::WinDirSep);

    // Report pass/fail
    if (idx != string::npos)
      PASSMSG(string("Found expected partial path. Path = ") + mypath);
    else
      FAILMSG("Did not find expected partial path. Expected path = " + mypath);
  }

  // value if not found

  string const mypath2 = getFilenameComponent(string("foobar"), rtt_dsxx::FC_PATH);
  string const expected = string(".") + string(1, rtt_dsxx::dirSep);
  if (mypath2 == expected)
    PASSMSG(string("FC_PATH: name w/o path successfully returned ") + expected);
  else
    FAILMSG("FC_PATH: name w/o path returned incorrect value.");

  // test the FC_NAME mode
  // ------------------------------------------------------------

  string const myname = getFilenameComponent(fqpn, rtt_dsxx::FC_NAME);

  idx = myname.find(string("tstPath"));
  if (idx != string::npos)
    PASSMSG(string("Found expected filename. myname = ") + myname);
  else
    FAILMSG("Did not find expected filename. Expected filename = " + myname);

  if (usesUnixDirSep) {
    if (mypath + myname == fqpn)
      PASSMSG(string("Successfully divided fqp into path+name = ") + mypath + myname);
    else
      FAILMSG(string("mypath+myname != fqpn") + string("\n\tmypath = ") + mypath +
              string("\n\tmyname = ") + myname + string("\n\tfqp    = ") + fqpn);
  }

  // value if not found

  string const myname2 = getFilenameComponent(string("foobar"), rtt_dsxx::FC_NAME);
  if (myname2 == string("foobar"))
    PASSMSG("name w/o path successfully returned");
  else
    FAILMSG("name w/o path returned incorrect value.");

  // test the FC_REALPATH
  // ------------------------------------------------------------
  string const realpath = getFilenameComponent(fqp, rtt_dsxx::FC_REALPATH);

#if defined(_MSC_VER)
  { // The binary should exist.  Windows does not provide an execute bit.

    if (realpath.size() > 0)
      PASSMSG("FC_REALPATH has length > 0.");
    else
      FAILMSG("FC_REALPATH has length <= 0.");

    // draco_getstat rpstatus( exeExists );

    if (std::ifstream(realpath.c_str()))
      PASSMSG("FC_REALPATH points to a valid executable.");
    else
      FAILMSG(string("FC_REALPATH is invalid or not executable.") + string("  realpath = ") +
              realpath);
  }
#else
  { // The binary should exist and marked by the filesystem as executable.

    if (usesUnixDirSep) {
      if (realpath.size() > 0)
        PASSMSG("FC_REALPATH has length > 0.");
      else
        FAILMSG("FC_REALPATH has length <= 0.");

      draco_getstat rpstatus(realpath);

      // string realpath2 = getFilenameComponent( "Makefile", rtt_dsxx::FC_REALPATH );
      // draco_getstat rpstatus2( realpath2 );
      // bool exebit = rpstatus2.has_permission_bit( 0100 );

      if (rpstatus.has_permission_bit(0100))
        PASSMSG(string("FC_REALPATH points to a valid executable. Path = ") + realpath);
      else
        FAILMSG(string("FC_REALPATH is invalid or not executable. Path = ") + realpath);
    } else {
      if (realpath.size() == 0)
        PASSMSG("FC_REALPATH has length <= 0.");
      else
        FAILMSG("FC_REALPATH has length > 0.");
    }
  }
#endif

  string const absolutepath = getFilenameComponent(fqp, rtt_dsxx::FC_ABSOLUTE);

  // These are not implemented yet
  // ------------------------------------------------------------

  bool caught = false;
  try {
    string const extension = getFilenameComponent(fqp, rtt_dsxx::FC_EXT);
  } catch (...) {
    caught = true;
    PASSMSG("FC_EXT throws.");
  }
  if (!caught)
    FAILMSG("FC_EXT failed to throw.");

  caught = false;
  try {
    string const extension = getFilenameComponent(fqp, rtt_dsxx::FC_NAME_WE);
  } catch (...) {
    caught = true;
    PASSMSG("FC_NAME_WE throws.");
  }
  if (!caught)
    FAILMSG("FC_NAME_WE failed to throw.");

  // FC_LASTVALUE should always throw.
  caught = false;
  try {
    string const foo = getFilenameComponent(fqp, rtt_dsxx::FC_LASTVALUE);
  } catch (...) {
    caught = true;
    PASSMSG("FC_LASTVALUE throws.");
  }
  if (!caught)
    FAILMSG("FC_LASTVALUE failed to throw.");

  return;
}

//------------------------------------------------------------------------------------------------//
void test_draco_remove(rtt_dsxx::ScalarUnitTest &ut) {
  std::cout << "\nBegin test tstPath::test_draco_remove.\n\n";
  {
    // create a file.
    std::cout << "Creating a file...\n";
    std::string const dummyFile("dummyFile.txt");
    std::ofstream outfile(dummyFile.c_str());
    outfile.close();

    // Did we create the file?
    if (std::filesystem::exists(dummyFile))
      PASSMSG(string("Successfully created file = ") + dummyFile);
    else {
      FAILMSG(string("Failed to create file = ") + dummyFile);
      // no reason to continue.
      return;
    }

    // Remove the file
    draco_dir_print(dummyFile);
    draco_remove_dir(dummyFile);
    if (std::filesystem::exists(dummyFile))
      FAILMSG(string("Failed to remove file = ") + dummyFile);
    else
      PASSMSG(string("Successfully removed file = ") + dummyFile);
  }

  {
    // Test a more complex sytem of directories and files.

    std::cout << "\nCreating files in a directory structure...\n\n";

    std::string const dummyFile1("dummydir/d1/dummyFile1.txt");
    std::string const dummyFile2("dummydir/d1/dummyFile2.txt");
    std::string const dummyFile3("dummydir/dummyFile3.txt");
    std::string const dummyDir1(getFilenameComponent(dummyFile1, FC_PATH));
    std::string const dummyDir2(getFilenameComponent(dummyDir1, FC_PATH));

    // Create directories.

    draco_mkdir(dummyDir2); // "dummydir"
    if (std::filesystem::is_directory(dummyDir2))
      PASSMSG(std::string("Successfully created directory = ") + dummyDir2);
    else
      FAILMSG(std::string("Failed to create directory = ") + dummyDir2);
    draco_mkdir(dummyDir1); // "dummydir/d1"
    if (std::filesystem::is_directory(dummyDir1))
      PASSMSG(std::string("Successfully created directory = ") + dummyDir1);
    else
      FAILMSG(std::string("Failed to create directory = ") + dummyDir1);

    // Create files

    std::ofstream outfile1(dummyFile1.c_str());
    std::ofstream outfile2(dummyFile2.c_str());
    std::ofstream outfile3(dummyFile3.c_str());
    outfile1.close();
    outfile2.close();
    outfile3.close();

    // Did we create the files?
    if (std::filesystem::exists(dummyFile1) && std::filesystem::exists(dummyFile2) &&
        std::filesystem::exists(dummyFile3)) {
      PASSMSG(string("Successfully created files = ") + dummyFile1 + std::string(", ") +
              dummyFile2 + std::string(" and ") + dummyFile3);
    } else {
      FAILMSG(string("Failed to create files = ") + dummyFile1 + std::string(", ") + dummyFile2 +
              std::string(" and ") + dummyFile3);
      // no reason to continue.
      return;
    }
    // Print the directory tree
    std::cout << "The directory tree contains: \n";
    draco_dir_print(dummyDir2);

    // Recursively remove the file and subdirectory.
    std::cout << "Removing all entries in " << dummyDir2 << '\n';
    draco_remove_dir(dummyDir2);
    if (std::filesystem::exists(dummyDir2))
      FAILMSG(string("Failed to remove directory = ") + dummyDir2);
    else
      PASSMSG(string("Successfully removed directory = ") + dummyDir2);
  }
  return;
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, release);
  try {
    test_currentPath(ut);
    test_getFilenameComponent(ut, string(argv[0]));
    test_draco_remove(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstPath.cc
//------------------------------------------------------------------------------------------------//
