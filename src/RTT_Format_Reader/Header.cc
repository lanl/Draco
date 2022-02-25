//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/Header.cc
 * \author B.T. Adams
 * \date   WED Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/Header class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Header.hh"
#include "ds++/Assert.hh"

namespace rtt_RTT_Format_Reader {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief ParseS the header data block from the mesh file via calls to private member functions.
 * \param meshfile Mesh file name.
 */
void Header::readHeader(ifstream &meshfile) {
  readKeyword(meshfile);
  readData(meshfile);
  readEndKeyword(meshfile);
}
/*!
 * \brief Reads and validates the header block keyword.
 * \param meshfile Mesh file name.
 */
void Header::readKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "header", "Invalid mesh file: Header block missing");
  std::getline(meshfile, dummyString); // read and discard blank line.
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the header block data.
 * \param meshfile Mesh file name.
 */
void Header::readData(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "version", "Invalid mesh file: Header block missing version");
  std::getline(meshfile, dummyString);
  // strip leading blanks and trailing comments and blanks from the version.
  if (dummyString.rfind('!') != string::npos)
    dummyString.erase(dummyString.rfind('!'));
  version = dummyString.substr(dummyString.find_first_not_of(' '),
                               dummyString.find_last_not_of(' ') -
                                   dummyString.find_first_not_of(' ') + 1);
  Insist(version == "v1.0.0", "Invalid mesh file: Wrong version");

  meshfile >> dummyString;
  Insist(dummyString == "title", "Invalid mesh file: Header block missing title");
  std::getline(meshfile, dummyString);
  // strip leading blanks and trailing comments and blanks from the title.
  if (dummyString.rfind('!') != string::npos)
    dummyString.erase(dummyString.rfind('!'));
  title = dummyString.substr(dummyString.find_first_not_of(' '),
                             dummyString.find_last_not_of(' ') -
                                 dummyString.find_first_not_of(' ') + 1);

  meshfile >> dummyString;
  Insist(dummyString == "date", "Invalid mesh file: Header block missing date");
  std::getline(meshfile, dummyString);
  // strip leading blanks and trailing comments and blanks from the date.
  if (dummyString.rfind('!') != string::npos)
    dummyString.erase(dummyString.rfind('!'));
  date = dummyString.substr(dummyString.find_first_not_of(' '),
                            dummyString.find_last_not_of(' ') - dummyString.find_first_not_of(' ') +
                                1);

  meshfile >> dummyString >> cycle;
  Insist(dummyString == "cycle", "Invalid mesh file: Header block missing cycle");
  std::getline(meshfile, dummyString);

  meshfile >> dummyString >> time;
  Insist(dummyString == "time", "Invalid mesh file: Header block missing time");
  std::getline(meshfile, dummyString);

  meshfile >> dummyString >> ncomments;
  Insist(dummyString == "ncomments", "Invalid mesh file: Header block missing ncomments");
  std::getline(meshfile, dummyString);
  comments.resize(ncomments);
  for (int i = 0; i < ncomments; ++i) {
    std::getline(meshfile, comments[i]);
    // strip leading blanks and trailing comments and blanks from the comment lines.
    if (comments[i].rfind('!') != string::npos)
      comments[i].erase(comments[i].rfind('!'));
    comments[i] = comments[i].substr(comments[i].find_first_not_of(' '),
                                     comments[i].find_last_not_of(' ') -
                                         comments[i].find_first_not_of(' ') + 1);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the end_header block keyword.
 * \param meshfile Mesh file name.
 */
void Header::readEndKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "end_header", "Invalid mesh file: Header block missing end");
  std::getline(meshfile, dummyString); // read and discard blank line.
}

} // end namespace rtt_RTT_Format_Reader

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/Header.cc
//------------------------------------------------------------------------------------------------//
