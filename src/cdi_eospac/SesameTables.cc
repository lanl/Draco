//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_eospac/SesameTables.cc
 * \author Kelly Thompson
 * \date   Fri Apr  6 08:57:48 2001
 * \brief  Implementation file for SesameTables (mapping material IDs to Sesame table indexes).
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "SesameTables.hh"
#include "ds++/Assert.hh"
#include "ds++/Packing_Utils.hh"
#include <iostream>
#include <sstream>

namespace rtt_cdi_eospac {

//------------------------------------------------------------------------------------------------//
// Construct from packed data stream
SesameTables::SesameTables(std::vector<char> const &packed)
    : numReturnTypes(EOS_M_DT + 1), //  EOS_M_DT = 305 (see eos_Interface.h)
      matMap(), rtMap(), tableName(initializeTableNames(numReturnTypes)),
      tableDescription(initializeTableDescriptions(numReturnTypes)) {
  // At least 3 integers (even if the data arrays are empty).
  // min size = numReturnTypes(unsigned) + packed_matmap_size(size_t) +
  //            packed_rtmap_size(size_t)
  Require(packed.size() >= sizeof(unsigned) + 2 * sizeof(size_t));

  // make an unpacker
  rtt_dsxx::Unpacker unpacker;
  unpacker.set_buffer(packed.size(), &packed[0]);

  // unpack and check the number of tables
  unsigned nrt(0);
  unpacker >> nrt;
  Check(nrt == numReturnTypes);

  // unpack and check matMap
  size_t packed_matmap_size(0);
  unpacker >> packed_matmap_size;
  Check(packed_matmap_size > 0);

  // provide container for unpacking and then unpack the data
  //! \bug update this for map<usinged, int>
  std::vector<char> packed_matmap(packed_matmap_size);
  for (size_t i = 0; i < packed_matmap_size; ++i)
    unpacker >> packed_matmap[i];
  rtt_dsxx::unpack_data(matMap, packed_matmap);

  // unpack and check rtMap
  size_t packed_rtmap_size(0);
  unpacker >> packed_rtmap_size;
  Check(packed_rtmap_size > 0);

  // provide container for unpacking and then unpack the data
  std::vector<char> packed_rtmap(packed_rtmap_size);
  for (size_t i = 0; i < packed_rtmap_size; ++i)
    unpacker >> packed_rtmap[i];
  rtt_dsxx::unpack_data(rtMap, packed_rtmap);

  return;
}

//------------------------------------------------------------------------------------------------//
// Set functions
//------------------------------------------------------------------------------------------------//

SesameTables &SesameTables::addTable(EOS_INTEGER const tableID, unsigned const matID) {
  // insert a new entry into the matMap.
  matMap[tableID] = matID;
  if (std::find(rtMap[matID].begin(), rtMap[matID].end(), tableID) == rtMap[matID].end())
    rtMap[matID].push_back(tableID);
  return *this;
}

SesameTables &SesameTables::Ue_DT(unsigned matID) { return addTable(EOS_Ue_DT, matID); }
SesameTables &SesameTables::Zfc_DT(unsigned matID) { return addTable(EOS_Zfc_DT, matID); }
SesameTables &SesameTables::Ktc_DT(unsigned matID) { return addTable(EOS_Ktc_DT, matID); }
SesameTables &SesameTables::Uic_DT(unsigned matID) { return addTable(EOS_Uic_DT, matID); }
SesameTables &SesameTables::Pt_DT(unsigned matID) { return addTable(EOS_Pt_DT, matID); }
SesameTables &SesameTables::Dv_T(unsigned matID) { return addTable(EOS_Dv_T, matID); }
SesameTables &SesameTables::Ogb(unsigned matID) { return addTable(EOS_Ogb, matID); }
SesameTables &SesameTables::T_DUe(unsigned matID) { return addTable(EOS_T_DUe, matID); }
SesameTables &SesameTables::T_DUic(unsigned matID) { return addTable(EOS_T_DUic, matID); }

// Move functions out as they are needed by new code and add unit tests for each.
#if 0

SesameTables& SesameTables::T_DUiz( unsigned matID )
{
    return addTable( EOS_T_DUiz, matID );
}
SesameTables& SesameTables::Ut_DT( unsigned matID )
{
    return addTable( EOS_Ut_DT, matID );
}
SesameTables& SesameTables::T_DPt( unsigned matID )
{
    return addTable( EOS_T_DPt, matID );
}
SesameTables& SesameTables::T_DUt( unsigned matID )
{
    return addTable( EOS_T_DUt, matID );
}
SesameTables& SesameTables::Pt_DUt( unsigned matID )
{
    return addTable( EOS_Pt_DUt, matID );
}
SesameTables& SesameTables::Ut_DPt( unsigned matID )
{
    return addTable( EOS_Ut_DPt, matID );
}
SesameTables& SesameTables::Pic_DT( unsigned matID )
{
    return addTable( EOS_Pic_DT, matID );
}
SesameTables& SesameTables::T_DPic( unsigned matID )
{
    return addTable( EOS_T_DPic, matID );
}
SesameTables& SesameTables::Pic_DUic( unsigned matID )
{
    return addTable( EOS_Pic_DUic, matID );
}
SesameTables& SesameTables::Uic_DPic( unsigned matID )
{
    return addTable( EOS_Uic_DPic, matID );
}
SesameTables& SesameTables::Pe_DT( unsigned matID )
{
    return addTable( EOS_Pe_DT, matID );
}
SesameTables& SesameTables::T_DPe( unsigned matID )
{
    return addTable( EOS_T_DPe, matID );
}

SesameTables& SesameTables::Pe_DUe( unsigned matID )
{
    return addTable( EOS_Pe_DUe, matID );
}
SesameTables& SesameTables::Ue_DPe( unsigned matID )
{
    return addTable( EOS_Ue_DPe, matID );
}
SesameTables& SesameTables::Pc_D( unsigned matID )
{
    return addTable( EOS_Pc_D, matID );
}
SesameTables& SesameTables::Uc_D( unsigned matID )
{
    return addTable( EOS_Uc_D, matID );
}
SesameTables& SesameTables::Kr_DT( unsigned matID )
{
    return addTable( EOS_Kr_DT, matID );
}
SesameTables& SesameTables::Keo_DT( unsigned matID )
{
    return addTable( EOS_Keo_DT, matID );
}
SesameTables& SesameTables::Zfo_DT( unsigned matID )
{
    return addTable( EOS_Zfo_DT, matID );
}
SesameTables& SesameTables::Kp_DT(  unsigned matID )
{
    return addTable( EOS_Kp_DT, matID );
}
SesameTables& SesameTables::Kec_DT( unsigned matID )
{
    return addTable( EOS_Kec_DT, matID );
}
SesameTables& SesameTables::B_DT( unsigned matID )
{
    return addTable( EOS_B_DT, matID );
}
SesameTables& SesameTables::Kc_DT( unsigned matID )
{
    return addTable( EOS_Kc_DT, matID );
}
SesameTables& SesameTables::Tm_D(  unsigned matID )
{
    return addTable( EOS_Tm_D, matID );
}
SesameTables& SesameTables::Pm_D(  unsigned matID )
{
    return addTable( EOS_Pm_D, matID );
}
SesameTables& SesameTables::Um_D(  unsigned matID )
{
    return addTable( EOS_Um_D, matID );
}
SesameTables& SesameTables::Tf_D( unsigned matID )
{
    return addTable( EOS_Tf_D, matID );
}
SesameTables& SesameTables::Pf_D( unsigned matID )
{
    return addTable( EOS_Pf_D, matID );
}
SesameTables& SesameTables::Uf_D( unsigned matID )
{
    return addTable( EOS_Uf_D, matID );
}
SesameTables& SesameTables::Gs_D( unsigned matID )
{
    return addTable( EOS_Gs_D, matID );
}
#endif

//------------------------------------------------------------------------------------------------//
// Get Functions

// Return the enumerated data type associated with the provided integer index
std::vector<EOS_INTEGER> SesameTables::returnTypes(unsigned const tableIndex) const {
  size_t const found = rtMap.count(tableIndex);
  std::vector<EOS_INTEGER> result;
  // note: map::operator[] is non-const only.
  if (found > 0)
    result = rtMap.find(tableIndex)->second;
  else {
    std::ostringstream msg;
    msg << "Requested tableIndex = " << tableIndex
        << ", does not not exist in the SesameTables object.  You must"
        << " assign a table index before attempting to access it." << std::endl;
    Insist(found > 0, msg.str());
  }
  return result;
}

unsigned SesameTables::matID(EOS_INTEGER returnType) const {
  Require(returnType >= 0);
  Require(matMap.count(returnType) == 1);
  // note: map::operator[] is non-const only.
  return matMap.find(returnType)->second;
}

//------------------------------------------------------------------------------------------------//
/*! Pack a SesameTables object into a vector<char> stream.
 *
 * Packed data stream:
 *
 * unsigned numReturnTypes
 * size_t   matMap.size()
 * int[]    matMap
 * size_t   rtMap.size()
 * map<unsigned,vector<int>> rtMap
 */
std::vector<char> SesameTables::pack() const {
  using std::string;
  using std::vector;

  // Size of packed SesameTables
  size_t packed_SesameTable_size(0);
  // The size of data member 'numReturnTypes'
  packed_SesameTable_size += sizeof(unsigned);

  // pack up the matMap
  vector<char> packed_matmap;
  rtt_dsxx::pack_data(matMap, packed_matmap);
  // packed data is an integer for the length of matMap plus the size of the actual data.
  packed_SesameTable_size += sizeof(size_t) + packed_matmap.size();

  // pack up the rtMap
  vector<char> packed_rtmap;
  rtt_dsxx::pack_data(rtMap, packed_rtmap);
  packed_SesameTable_size += sizeof(size_t) + packed_rtmap.size();

  // make a container to hold the packed data
  vector<char> packed(packed_SesameTable_size);

  // make a packer and set it
  rtt_dsxx::Packer packer;
  packer.set_buffer(packed_SesameTable_size, &packed[0]);

  // pack the numReturnTypes
  packer << numReturnTypes;

  // pack the matMap (size+data)
  packer << packed_matmap.size();
  for (char &c : packed_matmap)
    packer << c;

  // pack the rtMap(data)
  packer << packed_rtmap.size();
  for (char &c : packed_rtmap)
    packer << c;

  Ensure(packer.get_ptr() == &packed[0] + packed_SesameTable_size);
  return packed;
}

//------------------------------------------------------------------------------------------------//
void SesameTables::printEosTableList() const {
  std::cout << "List of EOS Tables:\n" << std::endl;
  for (size_t i = 0; i < tableName.size(); ++i)
    if (tableName[i].size() > 0)
      std::cout << tableName[i] << "\t- " << tableDescription[i] << "\n";
  return;
}

//------------------------------------------------------------------------------------------------//
std::vector<std::string> SesameTables::initializeTableNames(size_t const datasize) {
  // Create a mapping between the Enum and a string name
  std::vector<std::string> tableName(datasize);

  tableName[EOS_NullTable] = std::string("EOS_NullTable");
  tableName[EOS_Comment] = std::string("EOS_Comment");
  tableName[EOS_Info] = std::string("EOS_Info");
  tableName[EOS_Pt_DT] = std::string("EOS_Pt_DT");
  tableName[EOS_D_PtT] = std::string("EOS_D_PtT");
  tableName[EOS_T_DPt] = std::string("EOS_T_DPt");
  tableName[EOS_Pt_DUt] = std::string("EOS_Pt_DUt");
  tableName[EOS_Pt_DAt] = std::string("EOS_Pt_DAt");
  tableName[EOS_Pt_DSt] = std::string("EOS_Pt_DSt");
  tableName[EOS_Ut_DT] = std::string("EOS_Ut_DT");
  tableName[EOS_T_DUt] = std::string("EOS_T_DUt");
  tableName[EOS_Ut_DPt] = std::string("EOS_Ut_DPt");
  tableName[EOS_Ut_DAt] = std::string("EOS_Ut_DAt");
  tableName[EOS_Ut_DSt] = std::string("EOS_Ut_DSt");
  tableName[EOS_Ut_PtT] = std::string("EOS_Ut_PtT");
  tableName[EOS_At_DT] = std::string("EOS_At_DT");
  tableName[EOS_T_DAt] = std::string("EOS_T_DAt");
  tableName[EOS_At_DPt] = std::string("EOS_At_DPt");
  tableName[EOS_At_DUt] = std::string("EOS_At_DUt");
  tableName[EOS_At_DSt] = std::string("EOS_At_DSt");
  tableName[EOS_St_DT] = std::string("EOS_St_DT");
  tableName[EOS_T_DSt] = std::string("EOS_T_DSt");
  tableName[EOS_St_DPt] = std::string("EOS_St_DPt");
  tableName[EOS_St_DUt] = std::string("EOS_St_DUt");
  tableName[EOS_St_DAt] = std::string("EOS_St_DAt");
  tableName[EOS_Pic_DT] = std::string("EOS_Pic_DT");
  tableName[EOS_T_DPic] = std::string("EOS_T_DPic");
  tableName[EOS_Pic_DUic] = std::string("EOS_Pic_DUic");
  tableName[EOS_Pic_DAic] = std::string("EOS_Pic_DAic");
  tableName[EOS_Pic_DSic] = std::string("EOS_Pic_DSic");
  tableName[EOS_Uic_DT] = std::string("EOS_Uic_DT");
  tableName[EOS_T_DUic] = std::string("EOS_T_DUic");
  tableName[EOS_Uic_DPic] = std::string("EOS_Uic_DPic");
  tableName[EOS_Uic_DAic] = std::string("EOS_Uic_DAic");
  tableName[EOS_Uic_DSic] = std::string("EOS_Uic_DSic");
  tableName[EOS_Aic_DT] = std::string("EOS_Aic_DT");
  tableName[EOS_T_DAic] = std::string("EOS_T_DAic");
  tableName[EOS_Aic_DPic] = std::string("EOS_Aic_DPic");
  tableName[EOS_Aic_DUic] = std::string("EOS_Aic_DUic");
  tableName[EOS_Aic_DSic] = std::string("EOS_Aic_DSic");
  tableName[EOS_Sic_DT] = std::string("EOS_Sic_DT");
  tableName[EOS_T_DSic] = std::string("EOS_T_DSic");
  tableName[EOS_Sic_DPic] = std::string("EOS_Sic_DPic");
  tableName[EOS_Sic_DUic] = std::string("EOS_Sic_DUic");
  tableName[EOS_Sic_DAic] = std::string("EOS_Sic_DAic");
  tableName[EOS_Pe_DT] = std::string("EOS_Pe_DT");
  tableName[EOS_T_DPe] = std::string("EOS_T_DPe");
  tableName[EOS_Pe_DUe] = std::string("EOS_Pe_DUe");
  tableName[EOS_Pe_DAe] = std::string("EOS_Pe_DAe");
  tableName[EOS_Pe_DSe] = std::string("EOS_Pe_DSe");
  tableName[EOS_Ue_DT] = std::string("EOS_Ue_DT");
  tableName[EOS_T_DUe] = std::string("EOS_T_DUe");
  tableName[EOS_Ue_DPe] = std::string("EOS_Ue_DPe");
  tableName[EOS_Ue_DAe] = std::string("EOS_Ue_DAe");
  tableName[EOS_Ue_DSe] = std::string("EOS_Ue_DSe");
  tableName[EOS_Ae_DT] = std::string("EOS_Ae_DT");
  tableName[EOS_T_DAe] = std::string("EOS_T_DAe");
  tableName[EOS_Ae_DPe] = std::string("EOS_Ae_DPe");
  tableName[EOS_Ae_DUe] = std::string("EOS_Ae_DUe");
  tableName[EOS_Ae_DSe] = std::string("EOS_Ae_DSe");
  tableName[EOS_Se_DT] = std::string("EOS_Se_DT");
  tableName[EOS_T_DSe] = std::string("EOS_T_DSe");
  tableName[EOS_Se_DPe] = std::string("EOS_Se_DPe");
  tableName[EOS_Se_DUe] = std::string("EOS_Se_DUe");
  tableName[EOS_Se_DAe] = std::string("EOS_Se_DAe");
  tableName[EOS_Piz_DT] = std::string("EOS_Piz_DT");
  tableName[EOS_T_DPiz] = std::string("EOS_T_DPiz");
  tableName[EOS_Piz_DUiz] = std::string("EOS_Piz_DUiz");
  tableName[EOS_Piz_DAiz] = std::string("EOS_Piz_DAiz");
  tableName[EOS_Piz_DSiz] = std::string("EOS_Piz_DSiz");
  tableName[EOS_Uiz_DT] = std::string("EOS_Uiz_DT");
  tableName[EOS_T_DUiz] = std::string("EOS_T_DUiz");
  tableName[EOS_Uiz_DPiz] = std::string("EOS_Uiz_DPiz");
  tableName[EOS_Uiz_DAiz] = std::string("EOS_Uiz_DAiz");
  tableName[EOS_Uiz_DSiz] = std::string("EOS_Uiz_DSiz");
  tableName[EOS_Aiz_DT] = std::string("EOS_Aiz_DT");
  tableName[EOS_T_DAiz] = std::string("EOS_T_DAiz");
  tableName[EOS_Aiz_DPiz] = std::string("EOS_Aiz_DPiz");
  tableName[EOS_Aiz_DUiz] = std::string("EOS_Aiz_DUiz");
  tableName[EOS_Aiz_DSiz] = std::string("EOS_Aiz_DSiz");
  tableName[EOS_Siz_DT] = std::string("EOS_Siz_DT");
  tableName[EOS_T_DSiz] = std::string("EOS_T_DSiz");
  tableName[EOS_Siz_DPiz] = std::string("EOS_Siz_DPiz");
  tableName[EOS_Siz_DUiz] = std::string("EOS_Siz_DUiz");
  tableName[EOS_Siz_DAiz] = std::string("EOS_Siz_DAiz");
  tableName[EOS_Pc_D] = std::string("EOS_Pc_D");
  tableName[EOS_Uc_D] = std::string("EOS_Uc_D");
  tableName[EOS_Ac_D] = std::string("EOS_Ac_D");
  tableName[EOS_Pv_T] = std::string("EOS_Pv_T");
  tableName[EOS_T_Pv] = std::string("EOS_T_Pv");
  tableName[EOS_Pv_Dv] = std::string("EOS_Pv_Dv");
  tableName[EOS_Pv_Dls] = std::string("EOS_Pv_Dls");
  tableName[EOS_Pv_Uv] = std::string("EOS_Pv_Uv");
  tableName[EOS_Pv_Uls] = std::string("EOS_Pv_Uls");
  tableName[EOS_Pv_Av] = std::string("EOS_Pv_Av");
  tableName[EOS_Pv_Als] = std::string("EOS_Pv_Als");
  tableName[EOS_Dv_T] = std::string("EOS_Dv_T");
  tableName[EOS_T_Dv] = std::string("EOS_T_Dv");
  tableName[EOS_Dv_Pv] = std::string("EOS_Dv_Pv");
  tableName[EOS_Dv_Dls] = std::string("EOS_Dv_Dls");
  tableName[EOS_Dv_Uv] = std::string("EOS_Dv_Uv");
  tableName[EOS_Dv_Uls] = std::string("EOS_Dv_Uls");
  tableName[EOS_Dv_Av] = std::string("EOS_Dv_Av");
  tableName[EOS_Dv_Als] = std::string("EOS_Dv_Als");
  tableName[EOS_Dls_T] = std::string("EOS_Dls_T");
  tableName[EOS_T_Dls] = std::string("EOS_T_Dls");
  tableName[EOS_Dls_Pv] = std::string("EOS_Dls_Pv");
  tableName[EOS_Dls_Dv] = std::string("EOS_Dls_Dv");
  tableName[EOS_Dls_Uv] = std::string("EOS_Dls_Uv");
  tableName[EOS_Dls_Uls] = std::string("EOS_Dls_Uls");
  tableName[EOS_Dls_Av] = std::string("EOS_Dls_Av");
  tableName[EOS_Dls_Als] = std::string("EOS_Dls_Als");
  tableName[EOS_Uv_T] = std::string("EOS_Uv_T");
  tableName[EOS_T_Uv] = std::string("EOS_T_Uv");
  tableName[EOS_Uv_Pv] = std::string("EOS_Uv_Pv");
  tableName[EOS_Uv_Dv] = std::string("EOS_Uv_Dv");
  tableName[EOS_Uv_Dls] = std::string("EOS_Uv_Dls");
  tableName[EOS_Uv_Uls] = std::string("EOS_Uv_Uls");
  tableName[EOS_Uv_Av] = std::string("EOS_Uv_Av");
  tableName[EOS_Uv_Als] = std::string("EOS_Uv_Als");
  tableName[EOS_Uls_T] = std::string("EOS_Uls_T");
  tableName[EOS_T_Uls] = std::string("EOS_T_Uls");
  tableName[EOS_Uls_Pv] = std::string("EOS_Uls_Pv");
  tableName[EOS_Uls_Dv] = std::string("EOS_Uls_Dv");
  tableName[EOS_Uls_Dls] = std::string("EOS_Uls_Dls");
  tableName[EOS_Uls_Uv] = std::string("EOS_Uls_Uv");
  tableName[EOS_Uls_Av] = std::string("EOS_Uls_Av");
  tableName[EOS_Uls_Als] = std::string("EOS_Uls_Als");
  tableName[EOS_Av_T] = std::string("EOS_Av_T");
  tableName[EOS_T_Av] = std::string("EOS_T_Av");
  tableName[EOS_Av_Pv] = std::string("EOS_Av_Pv");
  tableName[EOS_Av_Dv] = std::string("EOS_Av_Dv");
  tableName[EOS_Av_Dls] = std::string("EOS_Av_Dls");
  tableName[EOS_Av_Uv] = std::string("EOS_Av_Uv");
  tableName[EOS_Av_Uls] = std::string("EOS_Av_Uls");
  tableName[EOS_Av_Als] = std::string("EOS_Av_Als");
  tableName[EOS_Als_T] = std::string("EOS_Als_T");
  tableName[EOS_T_Als] = std::string("EOS_T_Als");
  tableName[EOS_Als_Pv] = std::string("EOS_Als_Pv");
  tableName[EOS_Als_Dv] = std::string("EOS_Als_Dv");
  tableName[EOS_Als_Dls] = std::string("EOS_Als_Dls");
  tableName[EOS_Als_Uv] = std::string("EOS_Als_Uv");
  tableName[EOS_Als_Uls] = std::string("EOS_Als_Uls");
  tableName[EOS_Als_Av] = std::string("EOS_Als_Av");
  tableName[EOS_Tm_D] = std::string("EOS_Tm_D");
  tableName[EOS_D_Tm] = std::string("EOS_D_Tm");
  tableName[EOS_Tm_Pm] = std::string("EOS_Tm_Pm");
  tableName[EOS_Tm_Um] = std::string("EOS_Tm_Um");
  tableName[EOS_Tm_Am] = std::string("EOS_Tm_Am");
  tableName[EOS_Pm_D] = std::string("EOS_Pm_D");
  tableName[EOS_D_Pm] = std::string("EOS_D_Pm");
  tableName[EOS_Pm_Tm] = std::string("EOS_Pm_Tm");
  tableName[EOS_Pm_Um] = std::string("EOS_Pm_Um");
  tableName[EOS_Pm_Am] = std::string("EOS_Pm_Am");
  tableName[EOS_Um_D] = std::string("EOS_Um_D");
  tableName[EOS_D_Um] = std::string("EOS_D_Um");
  tableName[EOS_Um_Tm] = std::string("EOS_Um_Tm");
  tableName[EOS_Um_Pm] = std::string("EOS_Um_Pm");
  tableName[EOS_Um_Am] = std::string("EOS_Um_Am");
  tableName[EOS_Am_D] = std::string("EOS_Am_D");
  tableName[EOS_D_Am] = std::string("EOS_D_Am");
  tableName[EOS_Am_Tm] = std::string("EOS_Am_Tm");
  tableName[EOS_Am_Pm] = std::string("EOS_Am_Pm");
  tableName[EOS_Am_Um] = std::string("EOS_Am_Um");
  tableName[EOS_Tf_D] = std::string("EOS_Tf_D");
  tableName[EOS_D_Tf] = std::string("EOS_D_Tf");
  tableName[EOS_Tf_Pf] = std::string("EOS_Tf_Pf");
  tableName[EOS_Tf_Uf] = std::string("EOS_Tf_Uf");
  tableName[EOS_Tf_Af] = std::string("EOS_Tf_Af");
  tableName[EOS_Pf_D] = std::string("EOS_Pf_D");
  tableName[EOS_D_Pf] = std::string("EOS_D_Pf");
  tableName[EOS_Pf_Tf] = std::string("EOS_Pf_Tf");
  tableName[EOS_Pf_Uf] = std::string("EOS_Pf_Uf");
  tableName[EOS_Pf_Af] = std::string("EOS_Pf_Af");
  tableName[EOS_Uf_D] = std::string("EOS_Uf_D");
  tableName[EOS_D_Uf] = std::string("EOS_D_Uf");
  tableName[EOS_Uf_Tf] = std::string("EOS_Uf_Tf");
  tableName[EOS_Uf_Pf] = std::string("EOS_Uf_Pf");
  tableName[EOS_Uf_Af] = std::string("EOS_Uf_Af");
  tableName[EOS_Af_D] = std::string("EOS_Af_D");
  tableName[EOS_D_Af] = std::string("EOS_D_Af");
  tableName[EOS_Af_Tf] = std::string("EOS_Af_Tf");
  tableName[EOS_Af_Pf] = std::string("EOS_Af_Pf");
  tableName[EOS_Af_Uf] = std::string("EOS_Af_Uf");
  tableName[EOS_Gs_D] = std::string("EOS_Gs_D");
  tableName[EOS_D_Gs] = std::string("EOS_D_Gs");
  tableName[EOS_Ogb] = std::string("EOS_Ogb");
  tableName[EOS_Kr_DT] = std::string("EOS_Kr_DT");
  tableName[EOS_Keo_DT] = std::string("EOS_Keo_DT");
  tableName[EOS_Zfo_DT] = std::string("EOS_Zfo_DT");
  tableName[EOS_Kp_DT] = std::string("EOS_Kp_DT");
  tableName[EOS_Zfc_DT] = std::string("EOS_Zfc_DT");
  tableName[EOS_Kec_DT] = std::string("EOS_Kec_DT");
  tableName[EOS_Ktc_DT] = std::string("EOS_Ktc_DT");
  tableName[EOS_B_DT] = std::string("EOS_B_DT");
  tableName[EOS_Kc_DT] = std::string("EOS_Kc_DT");
  tableName[EOS_V_PtT] = std::string("EOS_V_PtT");
  tableName[EOS_M_DT] = std::string("EOS_M_DT");

  return tableName;
}

//------------------------------------------------------------------------------------------------//
// Initialize the tableDescriptions database
//------------------------------------------------------------------------------------------------//
std::vector<std::string> SesameTables::initializeTableDescriptions(size_t const datasize) {
  std::vector<std::string> tableDescription(datasize);

  tableDescription[EOS_NullTable] = std::string("null table");
  tableDescription[EOS_Comment] = std::string("Descriptive Comments");
  tableDescription[EOS_Info] = std::string("Atomic Number, Atomic Mass, Normal Density, Solid Bulk "
                                           "Modulus, Exchange Coefficient");
  tableDescription[EOS_Pt_DT] = std::string("Total Pressure (Density- and Temperature-dependent)");
  tableDescription[EOS_D_PtT] = std::string("Density (Total Pressure- and Temperature-dependent)");
  tableDescription[EOS_T_DPt] = std::string("Temperature (Density- and Total Pressure-dependent)");
  tableDescription[EOS_Pt_DUt] =
      std::string("Total Pressure (Density- and Total Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Pt_DAt] =
      std::string("Total Pressure (Density- and Total Specific-Free-Energy-dependent)");
  tableDescription[EOS_Pt_DSt] =
      std::string("Total Pressure (Density- and Total Specific-Entropy-dependent)");
  tableDescription[EOS_Ut_DT] =
      std::string("Total Specific-Internal-Energy (Density- and Temperature-dependent)");
  tableDescription[EOS_T_DUt] =
      std::string("Temperature (Density- and Total Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Ut_DPt] =
      std::string("Total Specific-Internal-Energy (Density- and Total Pressure-dependent)");
  tableDescription[EOS_Ut_DAt] = std::string("Total Specific-Internal-Energy (Density- and Total "
                                             "Specific-Free-Energy-dependent)");
  tableDescription[EOS_Ut_DSt] = std::string("Total Specific-Internal-Energy (Density- and Total "
                                             "Specific-Entropy-dependent)");
  tableDescription[EOS_Ut_PtT] = std::string("Total Specific-Internal-Energy (Total Pressure- and "
                                             "Temperature-dependent)");
  tableDescription[EOS_At_DT] =
      std::string("Total Specific-Free-Energy (Density- and Temperature-dependent)");
  tableDescription[EOS_T_DAt] =
      std::string("Temperature (Density- and Total Specific-Free-Energy-dependent)");
  tableDescription[EOS_At_DPt] =
      std::string("Total Specific-Free-Energy (Density- and Total Pressure-dependent)");
  tableDescription[EOS_At_DUt] = std::string("Total Specific-Free-Energy (Density- and Total "
                                             "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_At_DSt] = std::string("Total Specific-Free-Energy (Density- and Total "
                                             "Specific-Entropy-dependent)");
  tableDescription[EOS_St_DT] =
      std::string("Total Specific-Entropy (Density- and Temperature-dependent)");
  tableDescription[EOS_T_DSt] =
      std::string("Temperature (Density- and Total Specific-Entropy-dependent)");
  tableDescription[EOS_St_DPt] =
      std::string("Total Specific-Entropy (Density- and Total Pressure-dependent)");
  tableDescription[EOS_St_DUt] = std::string("Total Specific-Entropy (Density- and Total "
                                             "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_St_DAt] = std::string("Total Specific-Entropy (Density- and Total "
                                             "Specific-Free-Energy-dependent)");
  tableDescription[EOS_Pic_DT] = std::string("Ion Pressure plus Cold Curve Pressure (Density- and "
                                             "Temperature-dependent)");
  tableDescription[EOS_T_DPic] =
      std::string("Temperature (Density- and Ion Pressure plus Cold Curve "
                  "Pressure-dependent)");
  tableDescription[EOS_Pic_DUic] =
      std::string("Ion Pressure plus Cold Curve Pressure (Density- and Ion "
                  "Specific-Internal-Energy plus Cold Curve "
                  "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Pic_DAic] =
      std::string("Ion Pressure plus Cold Curve Pressure (Density- and Ion "
                  "Specific-Free-Energy plus Cold Curve Specific-Free-Energy-dependent)");
  tableDescription[EOS_Pic_DSic] =
      std::string("Ion Pressure plus Cold Curve Pressure (Density- and Ion "
                  "Pressure plus Cold Curve Specific-Entropy-dependent)");
  tableDescription[EOS_Uic_DT] =
      std::string("Ion Specific-Internal-Energy plus Cold Curve Specific-Internal-Energy "
                  "(Density- and Temperature-dependent)");
  tableDescription[EOS_T_DUic] =
      std::string("Temperature (Density- and Ion Specific-Internal-Energy plus "
                  "Cold Curve Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Uic_DPic] =
      std::string("Ion Specific-Internal-Energy plus Cold Curve Specific-Internal-Energy "
                  "(Density- and Ion Pressure plus Cold Curve Pressure-dependent)");
  tableDescription[EOS_Uic_DAic] =
      std::string("Ion Specific-Internal-Energy plus Cold Curve Specific-Internal-Energy "
                  "(Density- and Ion Specific-Free-Energy plus Cold Curve "
                  "Specific-Free-Energy-dependent)");
  tableDescription[EOS_Uic_DSic] =
      std::string("Ion Specific-Internal-Energy plus Cold Curve Specific-Internal-Energy "
                  "(Density- and Ion Pressure plus Cold Curve Specific-Entropy-dependent)");
  tableDescription[EOS_Aic_DT] =
      std::string("Ion Specific-Free-Energy plus Cold Curve "
                  "Specific-Free-Energy (Density- and Temperature-dependent)");
  tableDescription[EOS_T_DAic] =
      std::string("Temperature (Density- and Ion Specific-Free-Energy plus "
                  "Cold Curve Specific-Free-Energy-dependent)");
  tableDescription[EOS_Aic_DPic] =
      std::string("Ion Specific-Free-Energy plus Cold Curve Specific-Free-Energy (Density- "
                  "and Ion Pressure plus Cold Curve Pressure-dependent)");
  tableDescription[EOS_Aic_DUic] =
      std::string("Ion Specific-Free-Energy plus Cold Curve Specific-Free-Energy (Density- "
                  "and Ion Specific-Internal-Energy plus Cold Curve "
                  "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Aic_DSic] =
      std::string("Ion Specific-Free-Energy plus Cold Curve Specific-Free-Energy (Density- "
                  "and Ion Pressure plus Cold Curve Specific-Entropy-dependent)");
  tableDescription[EOS_Sic_DT] =
      std::string("Ion Pressure plus Cold Curve Specific-Entropy (Density- and "
                  "Temperature-dependent)");
  tableDescription[EOS_T_DSic] =
      std::string("Temperature (Density- and Ion Pressure plus Cold Curve "
                  "Specific-Entropy-dependent)");
  tableDescription[EOS_Sic_DPic] =
      std::string("Ion Pressure plus Cold Curve Specific-Entropy (Density- and "
                  "Ion Pressure plus Cold Curve Pressure-dependent)");
  tableDescription[EOS_Sic_DUic] =
      std::string("Ion Pressure plus Cold Curve Specific-Entropy (Density- and "
                  "Ion Specific-Internal-Energy plus Cold Curve "
                  "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Sic_DAic] =
      std::string("Ion Pressure plus Cold Curve Specific-Entropy (Density- and Ion "
                  "Specific-Free-Energy plus Cold Curve Specific-Free-Energy-dependent)");
  tableDescription[EOS_Pe_DT] =
      std::string("Electron Pressure (Density- and Temperature-dependent)");
  tableDescription[EOS_T_DPe] =
      std::string("Temperature (Density- and Electron Pressure-dependent)");
  tableDescription[EOS_Pe_DUe] = std::string("Electron Pressure (Density- and Electron "
                                             "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Pe_DAe] = std::string("Electron Pressure (Density- and Electron "
                                             "Specific-Free-Energy-dependent)");
  tableDescription[EOS_Pe_DSe] =
      std::string("Electron Pressure (Density- and Electron Specific-Entropy-dependent)");
  tableDescription[EOS_Ue_DT] =
      std::string("Electron Specific-Internal-Energy (Density- and Temperature-dependent)");
  tableDescription[EOS_T_DUe] =
      std::string("Temperature (Density- and Electron Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Ue_DPe] =
      std::string("Electron Specific-Internal-Energy (Density- and Electron "
                  "Pressure-dependent)");
  tableDescription[EOS_Ue_DAe] =
      std::string("Electron Specific-Internal-Energy (Density- and Electron "
                  "Specific-Free-Energy-dependent)");
  tableDescription[EOS_Ue_DSe] =
      std::string("Electron Specific-Internal-Energy (Density- and Electron "
                  "Specific-Entropy-dependent)");
  tableDescription[EOS_Ae_DT] =
      std::string("Electron Specific-Free-Energy (Density- and Temperature-dependent)");
  tableDescription[EOS_T_DAe] =
      std::string("Temperature (Density- and Electron Specific-Free-Energy-dependent)");
  tableDescription[EOS_Ae_DPe] = std::string("Electron Specific-Free-Energy (Density- and Electron "
                                             "Pressure-dependent)");
  tableDescription[EOS_Ae_DUe] = std::string("Electron Specific-Free-Energy (Density- and Electron "
                                             "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Ae_DSe] = std::string("Electron Specific-Free-Energy (Density- and Electron "
                                             "Specific-Entropy-dependent)");
  tableDescription[EOS_Se_DT] =
      std::string("Electron Specific-Entropy (Density- and Temperature-dependent)");
  tableDescription[EOS_T_DSe] =
      std::string("Temperature (Density- and Electron Specific-Entropy-dependent)");
  tableDescription[EOS_Se_DPe] =
      std::string("Electron Specific-Entropy (Density- and Electron Pressure-dependent)");
  tableDescription[EOS_Se_DUe] = std::string("Electron Specific-Entropy (Density- and Electron "
                                             "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Se_DAe] = std::string("Electron Specific-Entropy (Density- and Electron "
                                             "Specific-Free-Energy-dependent)");
  tableDescription[EOS_Piz_DT] =
      std::string("Ion Pressure Including Zero Point (Density- and Temperature-dependent)");
  tableDescription[EOS_T_DPiz] =
      std::string("Temperature (Density- and Ion Pressure Including Zero Point-dependent)");
  tableDescription[EOS_Piz_DUiz] =
      std::string("Ion Pressure Including Zero Point (Density- and Ion "
                  "Specific-Internal-Energy Including Zero Point-dependent)");
  tableDescription[EOS_Piz_DAiz] =
      std::string("Ion Pressure Including Zero Point (Density- and Ion "
                  "Specific-Free-Energy Including Zero Point-dependent)");
  tableDescription[EOS_Piz_DSiz] =
      std::string("Ion Pressure Including Zero Point (Density- and Ion "
                  "Pressure Including Zero Specific-Entropy-dependent)");
  tableDescription[EOS_Uiz_DT] =
      std::string("Ion Specific-Internal-Energy Including Zero Point (Density- "
                  "and Temperature-dependent)");
  tableDescription[EOS_T_DUiz] =
      std::string("Temperature (Density- and Ion Specific-Internal-Energy "
                  "Including Zero Point-dependent)");
  tableDescription[EOS_Uiz_DPiz] =
      std::string("Ion Specific-Internal-Energy Including Zero Point (Density- "
                  "and Ion Pressure Including Zero Point-dependent)");
  tableDescription[EOS_Uiz_DAiz] =
      std::string("Ion Specific-Internal-Energy Including Zero Point (Density- and Ion "
                  "Specific-Free-Energy Including Zero Point-dependent)");
  tableDescription[EOS_Uiz_DSiz] =
      std::string("Ion Specific-Internal-Energy Including Zero Point (Density- and Ion "
                  "Pressure Including Zero Specific-Entropy-dependent)");
  tableDescription[EOS_Aiz_DT] =
      std::string("Ion Specific-Free-Energy Including Zero Point (Density- and "
                  "Temperature-dependent)");
  tableDescription[EOS_T_DAiz] = std::string("Temperature (Density- and Ion Specific-Free-Energy "
                                             "Including Zero Point-dependent)");
  tableDescription[EOS_Aiz_DPiz] =
      std::string("Ion Specific-Free-Energy Including Zero Point (Density- and "
                  "Ion Pressure Including Zero Point-dependent)");
  tableDescription[EOS_Aiz_DUiz] =
      std::string("Ion Specific-Free-Energy Including Zero Point (Density- and Ion "
                  "Specific-Internal-Energy Including Zero Point-dependent)");
  tableDescription[EOS_Aiz_DSiz] =
      std::string("Ion Specific-Free-Energy Including Zero Point (Density- and "
                  "Ion Pressure Including Zero Specific-Entropy-dependent)");
  tableDescription[EOS_Siz_DT] =
      std::string("Ion Pressure Including Zero Specific-Entropy (Density- and "
                  "Temperature-dependent)");
  tableDescription[EOS_T_DSiz] =
      std::string("Temperature (Density- and Ion Pressure Including Zero "
                  "Specific-Entropy-dependent)");
  tableDescription[EOS_Siz_DPiz] =
      std::string("Ion Pressure Including Zero Specific-Entropy (Density- and "
                  "Ion Pressure Including Zero Point-dependent)");
  tableDescription[EOS_Siz_DUiz] =
      std::string("Ion Pressure Including Zero Specific-Entropy (Density- and Ion "
                  "Specific-Internal-Energy Including Zero Point-dependent)");
  tableDescription[EOS_Siz_DAiz] =
      std::string("Ion Pressure Including Zero Specific-Entropy (Density- and "
                  "Ion Specific-Free-Energy Including Zero Point-dependent)");
  tableDescription[EOS_Pc_D] = std::string("Pressure Cold Curve (Density-dependent)");
  tableDescription[EOS_Uc_D] =
      std::string("Specific-Internal-Energy Cold Curve (Density-dependent)");
  tableDescription[EOS_Ac_D] = std::string("Specific-Free-Energy Cold Curve (Density-dependent)");
  tableDescription[EOS_Pv_T] = std::string("Vapor Pressure (Temperature-dependent)");
  tableDescription[EOS_T_Pv] = std::string("Temperature (Vapor Pressure-dependent)");
  tableDescription[EOS_Pv_Dv] =
      std::string("Vapor Pressure (Vapor Density on coexistence line-dependent)");
  tableDescription[EOS_Pv_Dls] =
      std::string("Vapor Pressure (Liquid or Solid Density on coexistence line-dependent)");
  tableDescription[EOS_Pv_Uv] =
      std::string("Vapor Pressure (Vapor Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Pv_Uls] =
      std::string("Vapor Pressure (Liquid or Solid Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Pv_Av] =
      std::string("Vapor Pressure (Vapor Specific-Free-Energy-dependent)");
  tableDescription[EOS_Pv_Als] =
      std::string("Vapor Pressure (Liquid or Solid Specific-Free-Energy-dependent)");
  tableDescription[EOS_Dv_T] =
      std::string("Vapor Density on coexistence line (Temperature-dependent)");
  tableDescription[EOS_T_Dv] =
      std::string("Temperature (Vapor Density on coexistence line-dependent)");
  tableDescription[EOS_Dv_Pv] =
      std::string("Vapor Density on coexistence line (Vapor Pressure-dependent)");
  tableDescription[EOS_Dv_Dls] =
      std::string("Vapor Density on coexistence line (Liquid or Solid Density "
                  "on coexistence line-dependent)");
  tableDescription[EOS_Dv_Uv] = std::string("Vapor Density on coexistence line (Vapor "
                                            "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Dv_Uls] = std::string("Vapor Density on coexistence line (Liquid or Solid "
                                             "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Dv_Av] = std::string("Vapor Density on coexistence line (Vapor "
                                            "Specific-Free-Energy-dependent)");
  tableDescription[EOS_Dv_Als] = std::string("Vapor Density on coexistence line (Liquid or Solid "
                                             "Specific-Free-Energy-dependent)");
  tableDescription[EOS_Dls_T] =
      std::string("Liquid or Solid Density on coexistence line (Temperature-dependent)");
  tableDescription[EOS_T_Dls] =
      std::string("Temperature (Liquid or Solid Density on coexistence line-dependent)");
  tableDescription[EOS_Dls_Pv] =
      std::string("Liquid or Solid Density on coexistence line (Vapor Pressure-dependent)");
  tableDescription[EOS_Dls_Dv] =
      std::string("Liquid or Solid Density on coexistence line (Vapor Density "
                  "on coexistence line-dependent)");
  tableDescription[EOS_Dls_Uv] = std::string("Liquid or Solid Density on coexistence line (Vapor "
                                             "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Dls_Uls] =
      std::string("Liquid or Solid Density on coexistence line (Liquid or "
                  "Solid Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Dls_Av] = std::string("Liquid or Solid Density on coexistence line (Vapor "
                                             "Specific-Free-Energy-dependent)");
  tableDescription[EOS_Dls_Als] =
      std::string("Liquid or Solid Density on coexistence line (Liquid or "
                  "Solid Specific-Free-Energy-dependent)");
  tableDescription[EOS_Uv_T] =
      std::string("Vapor Specific-Internal-Energy (Temperature-dependent)");
  tableDescription[EOS_T_Uv] =
      std::string("Temperature (Vapor Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Uv_Pv] =
      std::string("Vapor Specific-Internal-Energy (Vapor Pressure-dependent)");
  tableDescription[EOS_Uv_Dv] = std::string("Vapor Specific-Internal-Energy (Vapor Density on "
                                            "coexistence line-dependent)");
  tableDescription[EOS_Uv_Dls] =
      std::string("Vapor Specific-Internal-Energy (Liquid or Solid Density on "
                  "coexistence line-dependent)");
  tableDescription[EOS_Uv_Uls] = std::string("Vapor Specific-Internal-Energy (Liquid or Solid "
                                             "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Uv_Av] =
      std::string("Vapor Specific-Internal-Energy (Vapor Specific-Free-Energy-dependent)");
  tableDescription[EOS_Uv_Als] = std::string("Vapor Specific-Internal-Energy (Liquid or Solid "
                                             "Specific-Free-Energy-dependent)");
  tableDescription[EOS_Uls_T] =
      std::string("Liquid or Solid Specific-Internal-Energy (Temperature-dependent)");
  tableDescription[EOS_T_Uls] =
      std::string("Temperature (Liquid or Solid Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Uls_Pv] =
      std::string("Liquid or Solid Specific-Internal-Energy (Vapor Pressure-dependent)");
  tableDescription[EOS_Uls_Dv] =
      std::string("Liquid or Solid Specific-Internal-Energy (Vapor Density on "
                  "coexistence line-dependent)");
  tableDescription[EOS_Uls_Dls] =
      std::string("Liquid or Solid Specific-Internal-Energy (Liquid or Solid "
                  "Density on coexistence line-dependent)");
  tableDescription[EOS_Uls_Uv] = std::string("Liquid or Solid Specific-Internal-Energy (Vapor "
                                             "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Uls_Av] = std::string("Liquid or Solid Specific-Internal-Energy (Vapor "
                                             "Specific-Free-Energy-dependent)");
  tableDescription[EOS_Uls_Als] =
      std::string("Liquid or Solid Specific-Internal-Energy (Liquid or Solid "
                  "Specific-Free-Energy-dependent)");
  tableDescription[EOS_Av_T] = std::string("Vapor Specific-Free-Energy (Temperature-dependent)");
  tableDescription[EOS_T_Av] = std::string("Temperature (Vapor Specific-Free-Energy-dependent)");
  tableDescription[EOS_Av_Pv] =
      std::string("Vapor Specific-Free-Energy (Vapor Pressure-dependent)");
  tableDescription[EOS_Av_Dv] =
      std::string("Vapor Specific-Free-Energy (Vapor Density on coexistence "
                  "line-dependent)");
  tableDescription[EOS_Av_Dls] =
      std::string("Vapor Specific-Free-Energy (Liquid or Solid Density on "
                  "coexistence line-dependent)");
  tableDescription[EOS_Av_Uv] =
      std::string("Vapor Specific-Free-Energy (Vapor Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Av_Uls] = std::string("Vapor Specific-Free-Energy (Liquid or Solid "
                                             "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Av_Als] = std::string("Vapor Specific-Free-Energy (Liquid or Solid "
                                             "Specific-Free-Energy-dependent)");
  tableDescription[EOS_Als_T] =
      std::string("Liquid or Solid Specific-Free-Energy (Temperature-dependent)");
  tableDescription[EOS_T_Als] =
      std::string("Temperature (Liquid or Solid Specific-Free-Energy-dependent)");
  tableDescription[EOS_Als_Pv] =
      std::string("Liquid or Solid Specific-Free-Energy (Vapor Pressure-dependent)");
  tableDescription[EOS_Als_Dv] =
      std::string("Liquid or Solid Specific-Free-Energy (Vapor Density on "
                  "coexistence line-dependent)");
  tableDescription[EOS_Als_Dls] =
      std::string("Liquid or Solid Specific-Free-Energy (Liquid or Solid "
                  "Density on coexistence line-dependent)");
  tableDescription[EOS_Als_Uv] = std::string("Liquid or Solid Specific-Free-Energy (Vapor "
                                             "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Als_Uls] =
      std::string("Liquid or Solid Specific-Free-Energy (Liquid or Solid "
                  "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Als_Av] = std::string("Liquid or Solid Specific-Free-Energy (Vapor "
                                             "Specific-Free-Energy-dependent)");
  tableDescription[EOS_Tm_D] = std::string("Melt Temperature (Density-dependent)");
  tableDescription[EOS_D_Tm] = std::string("Density (Melt Temperature-dependent)");
  tableDescription[EOS_Tm_Pm] = std::string("Melt Temperature (Melt Pressure-dependent)");
  tableDescription[EOS_Tm_Um] =
      std::string("Melt Temperature (Melt Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Tm_Am] =
      std::string("Melt Temperature (Melt Specific-Free-Energy-dependent)");
  tableDescription[EOS_Pm_D] = std::string("Melt Pressure (Density-dependent)");
  tableDescription[EOS_D_Pm] = std::string("Density (Melt Pressure-dependent)");
  tableDescription[EOS_Pm_Tm] = std::string("Melt Pressure (Melt Temperature-dependent)");
  tableDescription[EOS_Pm_Um] =
      std::string("Melt Pressure (Melt Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Pm_Am] = std::string("Melt Pressure (Melt Specific-Free-Energy-dependent)");
  tableDescription[EOS_Um_D] = std::string("Melt Specific-Internal-Energy (Density-dependent)");
  tableDescription[EOS_D_Um] = std::string("Density (Melt Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Um_Tm] =
      std::string("Melt Specific-Internal-Energy (Melt Temperature-dependent)");
  tableDescription[EOS_Um_Pm] =
      std::string("Melt Specific-Internal-Energy (Melt Pressure-dependent)");
  tableDescription[EOS_Um_Am] =
      std::string("Melt Specific-Internal-Energy (Melt Specific-Free-Energy-dependent)");
  tableDescription[EOS_Am_D] = std::string("Melt Specific-Free-Energy (Density-dependent)");
  tableDescription[EOS_D_Am] = std::string("Density (Melt Specific-Free-Energy-dependent)");
  tableDescription[EOS_Am_Tm] =
      std::string("Melt Specific-Free-Energy (Melt Temperature-dependent)");
  tableDescription[EOS_Am_Pm] = std::string("Melt Specific-Free-Energy (Melt Pressure-dependent)");
  tableDescription[EOS_Am_Um] =
      std::string("Melt Specific-Free-Energy (Melt Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Tf_D] = std::string("Freeze Temperature (Density-dependent)");
  tableDescription[EOS_D_Tf] = std::string("Density (Freeze Temperature-dependent)");
  tableDescription[EOS_Tf_Pf] = std::string("Freeze Temperature (Freeze Pressure-dependent)");
  tableDescription[EOS_Tf_Uf] =
      std::string("Freeze Temperature (Freeze Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Tf_Af] =
      std::string("Freeze Temperature (Freeze Specific-Free-Energy-dependent)");
  tableDescription[EOS_Pf_D] = std::string("Freeze Pressure (Density-dependent)");
  tableDescription[EOS_D_Pf] = std::string("Density (Freeze Pressure-dependent)");
  tableDescription[EOS_Pf_Tf] = std::string("Freeze Pressure (Freeze Temperature-dependent)");
  tableDescription[EOS_Pf_Uf] =
      std::string("Freeze Pressure (Freeze Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Pf_Af] =
      std::string("Freeze Pressure (Freeze Specific-Free-Energy-dependent)");
  tableDescription[EOS_Uf_D] = std::string("Freeze Specific-Internal-Energy (Density-dependent)");
  tableDescription[EOS_D_Uf] = std::string("Density (Freeze Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Uf_Tf] =
      std::string("Freeze Specific-Internal-Energy (Freeze Temperature-dependent)");
  tableDescription[EOS_Uf_Pf] =
      std::string("Freeze Specific-Internal-Energy (Freeze Pressure-dependent)");
  tableDescription[EOS_Uf_Af] = std::string("Freeze Specific-Internal-Energy (Freeze "
                                            "Specific-Free-Energy-dependent)");
  tableDescription[EOS_Af_D] = std::string("Freeze Specific-Free-Energy (Density-dependent)");
  tableDescription[EOS_D_Af] = std::string("Density (Freeze Specific-Free-Energy-dependent)");
  tableDescription[EOS_Af_Tf] =
      std::string("Freeze Specific-Free-Energy (Freeze Temperature-dependent)");
  tableDescription[EOS_Af_Pf] =
      std::string("Freeze Specific-Free-Energy (Freeze Pressure-dependent)");
  tableDescription[EOS_Af_Uf] = std::string("Freeze Specific-Free-Energy (Freeze "
                                            "Specific-Internal-Energy-dependent)");
  tableDescription[EOS_Gs_D] = std::string("Shear Modulus (Density-dependent)");
  tableDescription[EOS_D_Gs] = std::string("Density (Shear Modulus-dependent)");
  tableDescription[EOS_Ogb] = std::string(" Calculated versus Interpolated Opacity Grid Boundary");
  tableDescription[EOS_Kr_DT] =
      std::string("Rosseland Mean Opacity (Density- and Temperature-dependent)");
  tableDescription[EOS_Keo_DT] =
      std::string("Electron Conductive Opacity (Opacity Model) (Density- and "
                  "Temperature-dependent)");
  tableDescription[EOS_Zfo_DT] =
      std::string("Mean Ion Charge (Opacity Model) (Density- and Temperature-dependent)");
  tableDescription[EOS_Kp_DT] =
      std::string("Planck Mean Opacity (Density- and Temperature-dependent)");
  tableDescription[EOS_Zfc_DT] = std::string("Mean Ion Charge (Conductivity Model) (Density- and "
                                             "Temperature-dependent)");
  tableDescription[EOS_Kec_DT] =
      std::string("Electrical Conductivity (Density- and Temperature-dependent)");
  tableDescription[EOS_Ktc_DT] =
      std::string("Thermal Conductivity (Density- and Temperature-dependent)");
  tableDescription[EOS_B_DT] =
      std::string("Thermoelectric Coefficient (Density- and Temperature-dependent)");
  tableDescription[EOS_Kc_DT] =
      std::string("Electron Conductive Opacity (Conductivity Model) (Density- "
                  "and Temperature-dependent)");
  tableDescription[EOS_V_PtT] =
      std::string("Specific-Volume (Pressure- and Temperature-dependent)");
  tableDescription[EOS_M_DT] = std::string("Mass Fraction (Density- and Temperature-dependent)");

  return tableDescription;
}

} // end namespace rtt_cdi_eospac

//------------------------------------------------------------------------------------------------//
// end of SesameTables.cc
//------------------------------------------------------------------------------------------------//
