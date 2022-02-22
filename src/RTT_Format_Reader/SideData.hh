//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/SideData.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/SideData class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_SideData_hh
#define rtt_RTT_Format_Reader_SideData_hh

#include "Sides.hh"

namespace rtt_RTT_Format_Reader {

//================================================================================================//
/*!
 * \brief Controls parsing, storing, and accessing the data specific to the sidedata block of the
 *        mesh file.
 */
//================================================================================================//
class SideData {
  // typedefs
  using ifstream = std::ifstream;
  using string = std::string;
  using vector_dbl = std::vector<double>;
  using vector_vector_dbl = std::vector<std::vector<double>>;

  const Dims &dims;
  vector_vector_dbl data;

public:
  SideData(const Dims &dims_)
      : dims(dims_), data(dims.get_nsides(), vector_dbl(dims.get_nside_data())) { /* empty */
  }
  ~SideData() = default;

  void readSideData(ifstream &meshfile);

private:
  void readKeyword(ifstream &meshfile);
  void readData(ifstream &meshfile);
  void readEndKeyword(ifstream &meshfile);

public:
  /*!
   * \brief Returns all of the data field values for each of the sides.
   * \return The data field values for each of the sides.
   */
  vector_vector_dbl get_data() const { return data; }

  /*!
   * \brief Returns all of the data field values for the specified side.
   * \param side_numb Side number.
   * \return The side data field values.
   */
  vector_dbl get_data(size_t side_numb) const { return data[side_numb]; }

  /*!
   * \brief Returns the specified data field value for the specified side.
   * \param side_numb Side number.
   * \param data_index Data field.
   * \return The side data field value.
   */
  double get_data(size_t side_numb, size_t data_index) const { return data[side_numb][data_index]; }
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_SideData

//------------------------------------------------------------------------------------------------//
// ned of RTT_Format_Reader/SideData.hh
//------------------------------------------------------------------------------------------------//
