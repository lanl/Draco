//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_eospac/SesameTables.hh
 * \author Kelly Thompson
 * \date   Fri Apr  6 08:57:48 2001
 * \brief  Header file for SesameTables (mapping material IDs to Sesame table indexes).
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_eospac_SesameTables_hh
#define rtt_cdi_eospac_SesameTables_hh

#include "eos_Interface.h"
#include <map>
#include <string>
#include <vector>

namespace rtt_cdi_eospac {

//================================================================================================//
/*!
 * \class SesameTables
 *
 * \brief This is a helper class for Eospac.  It tells Eospac what Sesame data is being requested
 *        and what lookup tables to use.
 *
 * \sa The web page for <a href="http://xweb.lanl.gov/PROJECTS/DATA/eos/">EOSPAC</a>.
 *
 * Each sesame material definition has 16 data tables (actually material identifiers) that define
 * its state.  At least one table must be defined for this to be a valid object.  This list of
 * tables is used by the Eospac constructor to determine what Sesame table data to cache.  There are
 * 37 return types defined by EOSPAC.  Some data tables provide information for more than one return
 * type.
 *
 * \sa cdi_eospac/test/tEospac.cc
 */
//================================================================================================//

class SesameTables {
  // DATA

  //! There are 305 return types defined by EOSPAC (see eos_Interface.h)
  unsigned const numReturnTypes; // should be 305

  /*!
   * \brief Map from EOSPAC data type to material identifier.
   *
   * Each of the enumerated EOSPAC data types can have a different SesameTable material identifier.
   * This vector contains a list of these material IDs.
   */
  std::map<EOS_INTEGER, unsigned> matMap;

  //! Toggle list to identify which data types have been requested (materialID, data type)
  std::map<unsigned, std::vector<EOS_INTEGER>> rtMap;

public:
  // DATA

  //! Return the table name based on the enum.
  std::vector<std::string> const tableName;

  //! Return a description of the table based on the enum
  std::vector<std::string> const tableDescription;

  // CREATORS

  //! Default constructor
  SesameTables()
      : numReturnTypes(EOS_M_DT + 1), //  EOS_M_DT = 305 (see eos_Interface.h)
        matMap(), rtMap(), tableName(initializeTableNames(numReturnTypes)),
        tableDescription(initializeTableDescriptions(numReturnTypes)) { /* empty */
  }

  //! Construct by unpacking a vector<char> stream.
  explicit SesameTables(std::vector<char> const &packed);

  // ACCESSORS

  // set functions: Uncomment when these are needed and unit tests are added.

  SesameTables &addTable(EOS_INTEGER const tableID, unsigned const matID);

  //! Thermal Conductivity (1/cm/s)
  SesameTables &Ktc_DT(unsigned matID);
  //! Electron Specific-Internal-Energy (MJ/kg)
  SesameTables &Ue_DT(unsigned matID);
  //! Ion Specific-Internal-Energy plus Cold Curve Specific-Internal-Energy (MJ/kg)
  SesameTables &Uic_DT(unsigned matID);
  //! Mean Ion Charge (Conductivity Model) (free electrons per atom)
  SesameTables &Zfc_DT(unsigned matID);
  //! Total Pressure (GPa)
  SesameTables &Pt_DT(unsigned matID);
  //! Vapor Density on coexistence line (Mg/m^3) from table 401
  SesameTables &Dv_T(unsigned matID);
  //! Calculated versus Interpolated Opacity Grid Boundary (from table 501)
  SesameTables &Ogb(unsigned matID);
  //! Temperature (K)
  SesameTables &T_DUe(unsigned matID);
  //! Temperature (K)
  SesameTables &T_DUic(unsigned matID);

#if 0
  //! Temperature (K)
  SesameTables& T_DUiz( unsigned matID );
  //! Thermoelectric Coefficient (1/cm^2/s)
  SesameTables& B_DT( unsigned matID );
  //! Shear Modulus (Gpa)
  SesameTables& Gs_D( unsigned matID );
  //! Electron Conductive Opacity (cm^2/g)
  SesameTables& Kc_DT( unsigned matID );
  //! Electrical Conductivity (1/s)
  SesameTables& Kec_DT( unsigned matID );
  //! Electron Conductive Opacity (cm^2/g)
  SesameTables& Keo_DT( unsigned matID );
  //! Planck Mean Opacity (cm^2/g)
  SesameTables& Kp_DT(  unsigned matID );
  //! Rosseland Mean Opacity (cm^2/g)
  SesameTables& Kr_DT( unsigned matID );
  //! Pressure Cold Curve (GPa)
  SesameTables& Pc_D( unsigned matID );
  //! Electron Pressure (GPa)
  SesameTables& Pe_DT( unsigned matID );
  //! Electron Pressure (GPa)
  SesameTables& Pe_DUe( unsigned matID );
  //! Freeze Pressure (GPa)
  SesameTables& Pf_D( unsigned matID );
  //! Ion Pressure plus Cold Curve Pressure (GPa)
  SesameTables& Pic_DT( unsigned matID );
  //! Ion Pressure plus Cold Curve Pressure (GPa)
  SesameTables& Pic_DUic( unsigned matID );
  //! Melt Pressure (GPa)
  SesameTables& Pm_D(  unsigned matID );
  //! Total Pressure (Gpa)
  SesameTables& Pt_DUt( unsigned matID );
  //! Temperature (K)
  SesameTables& T_DPe( unsigned matID );
  //! Temperature (K)
  SesameTables& T_DPic( unsigned matID );
  //! Temperature (K)
  SesameTables& T_DPt( unsigned matID );
  //! Temperature (K)
  SesameTables& T_DUt( unsigned matID );
  //! Freeze Temperature (eV)
  SesameTables& Tf_D( unsigned matID );
  //! Melt Temperature (K)
  SesameTables& Tm_D(  unsigned matID );
  //! Specific-Internal-Energy Cold Curve (MJ/kg)
  SesameTables& Uc_D( unsigned matID );
  //! Electron Specific-Internal-Energy (MJ/kg)
  SesameTables& Ue_DPe( unsigned matID );
  //! Freeze Specific-Internal-Energy (MJ/kg)
  SesameTables& Uf_D( unsigned matID );
  /*! \brief Ion Specific-Internal-Energy plus Cold Curve
   * Specific-Internal-Energy (MJ/kg) */
  SesameTables& Uic_DPic( unsigned matID );
  //! Melt Specific-Internal-Energy (MJ/kg)
  SesameTables& Um_D(  unsigned matID );
  //! Melt Specific-Internal-Energy (MJ/kg)
  SesameTables& Ut_DPt( unsigned matID );
  //! Melt Specific-Internal-Energy (MJ/kg)
  SesameTables& Ut_DT( unsigned matID );
  //! Mean Ion Charge (Opacity Model) (free electrons per atom)
  SesameTables& Zfo_DT( unsigned matID );
#endif

  // Get functions

  //! Return the material identifier associated with a Sesame return type.
  unsigned matID(EOS_INTEGER returnType) const;

  //! Return the enumerated data type associated with the integer index.
  std::vector<EOS_INTEGER> returnTypes(unsigned const index) const;

  //! Return a list of table/material identifiers associated with this SesameTables object.
  std::vector<unsigned> matList() const {
    std::vector<unsigned> mlist;
    for (auto const &it : rtMap)
      mlist.push_back(it.first);
    return mlist;
  }

  //! Pack a SesameTables object into a vector<char> stream.
  std::vector<char> pack() const;

  //! Print a list of EOS tables that EOSPAC knows about
  void printEosTableList() const;

  // implementation
  static std::vector<std::string> initializeTableNames(size_t const datasize);
  static std::vector<std::string> initializeTableDescriptions(size_t const datasize);

}; // end class SesameTables

} // end namespace rtt_cdi_eospac

#endif // rtt_cdi_eospac_SesameTables_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_eospac/SesameTables.hh
//------------------------------------------------------------------------------------------------//
