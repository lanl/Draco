//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/test/TestRTTFormatReader.hh
 * \author B.T. Adams
 * \date   Tue Mar 14 09:48:00 2000
 * \brief  Header file for the RTT_Format_Reader class unit test.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_test_TestRTT_Format_Reader_hh
#define rtt_test_TestRTT_Format_Reader_hh

#include "RTT_Format_Reader/RTT_Format_Reader.hh"
#include "ds++/UnitTest.hh"
#include <map>

using rtt_dsxx::UnitTest;
using RTT_Format_Reader = rtt_RTT_Format_Reader::RTT_Format_Reader;

enum Meshes { DEFINED, MESHES_LASTENTRY };

// All function tests with the exception of check_header and check_dims require that the Dims data
// has been properly processed, and use the verify_Dims function to query this data to determine if
// this is true.
extern std::map<Meshes, bool> Dims_validated;

bool verify_Dims(const RTT_Format_Reader &mesh, const Meshes &meshtype);

bool check_header(const RTT_Format_Reader &mesh, const Meshes &meshtype, UnitTest &ut);
bool check_dims(const RTT_Format_Reader &mesh, const Meshes &meshtype, UnitTest &ut);
bool check_node_flags(const RTT_Format_Reader &mesh, const Meshes &meshtype, UnitTest &ut);
bool check_side_flags(const RTT_Format_Reader &mesh, const Meshes &meshtype, UnitTest &ut);
bool check_cell_flags(const RTT_Format_Reader &mesh, const Meshes &meshtype, UnitTest &ut);
bool check_node_data_ids(const RTT_Format_Reader &mesh, const Meshes &meshtype, UnitTest &ut);
bool check_side_data_ids(const RTT_Format_Reader &mesh, const Meshes &meshtype, UnitTest &ut);
bool check_cell_data_ids(const RTT_Format_Reader &mesh, const Meshes &meshtype, UnitTest &ut);
bool check_cell_defs(const RTT_Format_Reader &mesh, const Meshes &meshtype, UnitTest &ut);
bool check_nodes(const RTT_Format_Reader &mesh, const Meshes &meshtype, UnitTest &ut);
bool check_sides(const RTT_Format_Reader &mesh, const Meshes &meshtype, UnitTest &ut);
bool check_cells(const RTT_Format_Reader &mesh, const Meshes &meshtype, UnitTest &ut);
bool check_node_data(const RTT_Format_Reader &mesh, const Meshes &meshtype, UnitTest &ut);
bool check_side_data(const RTT_Format_Reader &mesh, const Meshes &meshtype, UnitTest &ut);
bool check_cell_data(const RTT_Format_Reader &mesh, const Meshes &meshtype, UnitTest &ut);

#endif // _test_TestRTT_Format_Reader_hh

//------------------------------------------------------------------------------------------------//
// end of meshReaders/test/TestRTTFormatReader.hh
//------------------------------------------------------------------------------------------------//
