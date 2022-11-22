//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/Header.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/Header class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_Header_hh
#define rtt_RTT_Format_Reader_Header_hh

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace rtt_RTT_Format_Reader {

//================================================================================================//
/*!
 * \brief Controls parsing, storing, and accessing the data contained in the header block of the
 *        mesh file.
 */
//================================================================================================//
class Header {

private:
  void readKeyword(std::ifstream &meshfile);
  void readData(std::ifstream &meshfile);
  void readEndKeyword(std::ifstream &meshfile);

  // typedefs
  using ifstream = std::ifstream;
  using string = std::string;
  using vector_str = std::vector<string>;

  string version;
  string title;
  string date;
  int cycle{0};
  double time{0.0};
  int ncomments{0};
  vector_str comments;

public:
  Header()
      : version(std::string("")), title(std::string("")), date(std::string("")),
        comments(std::vector<string>()) { /*empty*/
  }
  ~Header() = default;
  Header(Header const &rhs) = delete;
  Header(Header &&rhs) noexcept = delete;
  Header &operator=(Header const &rhs) = delete;
  Header &operator=(Header &&rhs) noexcept = delete;

  void readHeader(ifstream &meshfile);

  /*!
   * \brief Returns the mesh file version number.
   * \return Version number.
   */
  string get_version() const { return version; }
  /*!
   * \brief Returns the mesh file title.
   * \return Title.
   */
  string get_title() const { return title; }
  /*!
   * \brief Returns the mesh file date.
   * \return Date the mesh file was generated.
   */
  string get_date() const { return date; }
  /*!
   * \brief Returns the mesh file cycle number.
   * \return Cycle number.
   */
  int get_cycle() const { return cycle; }
  /*!
   * \brief Returns the mesh file problem time.
   * \return Problem time.
   */
  double get_time() const { return time; }
  /*!
   * \brief Returns the number of comment lines in the mesh file.
   * \return The number of comment lines.
   */
  int get_ncomments() const { return ncomments; }
  /*!
   * \brief Returns the specified comment line from the mesh file.
   * \param i Line number of the comment to be returned.
   * \return The comment line.
   */
  string get_comments(size_t const i) const { return comments[i]; }
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_Header_hh

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/Header.hh
//------------------------------------------------------------------------------------------------//
