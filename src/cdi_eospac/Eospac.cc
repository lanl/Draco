//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_eospac/Eospac.cc
 * \author Kelly Thompson
 * \date   Mon Apr  2 14:14:29 2001
 * \brief  Implementation for Eospac class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Eospac.hh"
#include "EospacException.hh"
#include "ds++/Assert.hh"
#include "ds++/Packing_Utils.hh"
#include <array>
#include <iostream>
#include <sstream>
#include <string>

namespace rtt_cdi_eospac {

// ------------ //
// Constructors //
// ------------ //

//------------------------------------------------------------------------------------------------//
/*!
 * \brief The constructor for Eospac.
 *
 * \sa The definition of rtt_cdi_eospac::SesameTables.
 *
 * \param[in] in_SesTabs A rtt_cdi_eospac::SesameTables object that defines what data tables will be
 *        available for queries from the Eospac object.
 */
Eospac::Eospac(SesameTables in_SesTabs)
    : SesTabs(std::move(in_SesTabs)), matIDs(), returnTypes(), tableHandles(),
      infoItems(initializeInfoItems()), infoItemDescriptions(initializeInfoItemDescriptions()) {
  // Eospac can only be instantiated if SesameTables is provided.  If SesameTables is invalid this
  // will be caught in expandEosTable();

  // PreCache the default data type
  expandEosTable();

} // end Eospac::Eospac()

//------------------------------------------------------------------------------------------------//
//! Construct an Eospac by unpacking a vector<char> stream
Eospac::Eospac(std::vector<char> const &packed)
    : SesTabs(SesameTables(packed)), matIDs(), returnTypes(), tableHandles(),
      infoItems(initializeInfoItems()), infoItemDescriptions(initializeInfoItemDescriptions()) {
  // Use the initializer list to build the needed SesTabs.  Now initialize libeospac using the data
  // found in SesTabs just like in the default ctor.
  expandEosTable();
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Default Eospac() destructor.
 *
 * This is required to correctly release memeroy when an Eospac object is destroyed.  We define the
 * destructor in the implementation file to avoid including the unnecessary header files.
 */
Eospac::~Eospac() {
  // Destroy all data objects:
  EOS_INTEGER errorCode(EOS_OK);
  eos_DestroyAll(&errorCode);
  if (errorCode != EOS_OK) {
    std::ostringstream outputString;
    for (size_t i = 0; i < returnTypes.size(); ++i) {
      std::array<EOS_CHAR, EOS_MaxErrMsgLen> errorMessage{};
      EOS_INTEGER tableHandleErrorCode = EOS_OK;
      eos_GetErrorCode(&tableHandles[i], &tableHandleErrorCode);
      eos_GetErrorMessage(&tableHandleErrorCode, errorMessage.data());

      outputString << "\n\tAn unsuccessful request was made to destroy the "
                   << "EOSPAC table area by ~Eospac().\n"
                   << "\tThe error code returned by eos_DestroyAll(...) was \""
                   << tableHandleErrorCode << "\".\n"
                   << "\tThe associated error message is:\n\t\"" << std::string(errorMessage.data())
                   << "\"\n";
    }
    // Never throw an exception from the destructor.  This can cause confusion during stack
    // unwiding.
    std::cerr << outputString.str() << std::endl;
  }
}

// --------- //
// Accessors //
// --------- //

void Eospac::printTableInformation(EOS_INTEGER const tableType, std::ostream &out) const {
  // Obtain the table handle for this type
  EOS_INTEGER tableHandle(tableHandles[tableIndex(tableType)]);

  EOS_INTEGER one(1);
  EOS_INTEGER errorCode(EOS_OK);
  EOS_INTEGER ok(EOS_OK);
  EOS_INTEGER invalid_info_flag(EOS_INVALID_INFO_FLAG);
  EOS_REAL infoVal;
  EOS_BOOLEAN match1;
  EOS_BOOLEAN match2;

  out << "\nEOSPAC information for Table " << SesTabs.tableName[tableType] << " ("
      << SesTabs.tableDescription[tableType] << ")\n"
      << "------------------------------------------------------------------------------------\n";

  // There are 11 descriptions available for all tables.
  size_t numItems(infoItems.size());
  // There are extra descriptions available for non inverted tables (infoItems.size()).
  if (tableType == EOS_T_DUe || tableType == EOS_T_DUic)
    numItems = 11;

  for (size_t i = 0; i < numItems; ++i) {
    eos_GetTableInfo(&tableHandle, &one, &infoItems[i], &infoVal, &errorCode);

    eos_ErrorCodesEqual(&errorCode, &ok, &match1);
    eos_ErrorCodesEqual(&errorCode, &invalid_info_flag, &match2);
    if (match1 == EOS_TRUE) {
      out << std::setiosflags(std::ios::fixed) << std::setw(70) << std::left
          << infoItemDescriptions[i] << ": " << std::setprecision(6)
          << std::setiosflags(std::ios::fixed) << std::setw(13) << std::right << infoVal
          << std::endl;
    } else if (match2 == EOS_FALSE) {
      std::ostringstream outputString;
      std::array<EOS_CHAR, EOS_MaxErrMsgLen> errorMessage{};
      // Ignore EOS_INVALID_INFO_FLAG since not all infoItems are currently applicable to a specific
      // tableHandle.
      eos_GetErrorMessage(&errorCode, errorMessage.data());
      outputString << "\n\tAn unsuccessful request for EOSPAC table information "
                   << "was made by eos_GetTableInfo().\n"
                   << "\tThe requested infoType was \"" << infoItems[i]
                   << "\" (see eos_Interface.h for type)\n"
                   << "\tThe error code returned was \"" << errorCode << "\".\n"
                   << "\tThe associated error message is:\n\t\"" << std::string(errorMessage.data())
                   << "\"\n";
      throw EospacException(outputString.str());
    }
  }

  return;
}

//------------------------------------------------------------------------------------------------//
double Eospac::getSpecificElectronInternalEnergy(double temperature, double density) const {
  EOS_INTEGER const returnType = EOS_Ue_DT; // ES4enelc;
  // Convert temperatures from keV to degrees Kelvin.
  double vtempsKelvin = keV2K(temperature);
  return getF(dbl_v1(density), dbl_v1(vtempsKelvin), returnType, ETDD_VALUE)[0];
}

//------------------------------------------------------------------------------------------------//
std::vector<double>
Eospac::getSpecificElectronInternalEnergy(std::vector<double> const &vtemperature,
                                          std::vector<double> const &vdensity) const {
  EOS_INTEGER const returnType = EOS_Ue_DT; // ES4enelc;
  // Convert temperatures from keV to degrees Kelvin.
  std::vector<double> vtempsKelvin = vtemperature;
  std::transform(vtemperature.begin(), vtemperature.end(), vtempsKelvin.begin(), keV2K);
  return getF(vdensity, vtempsKelvin, returnType, ETDD_VALUE);
}

//------------------------------------------------------------------------------------------------//
double Eospac::getElectronHeatCapacity(double temperature, double density) const {
  // specific Heat capacity is dE/dT at constant pressure.  To obtain the specific electron heat
  // capacity we load the specific electron internal energy (E) and it's first derivative w.r.t
  // temperature.

  // Convert temperatures from keV to degrees Kelvin.
  double vtempsKelvin = keV2K(temperature);
  EOS_INTEGER const returnType = EOS_Ue_DT; // ES4enelc
  std::vector<double> Cve = getF(dbl_v1(density), dbl_v1(vtempsKelvin), returnType, ETDD_DFDY);
  // Convert back to Temperature units in keV
  std::transform(Cve.begin(), Cve.end(), Cve.begin(), keV2K);
  return Cve[0];
}

//------------------------------------------------------------------------------------------------//
std::vector<double> Eospac::getElectronHeatCapacity(std::vector<double> const &vtemperature,
                                                    std::vector<double> const &vdensity) const {
  // specific Heat capacity is dE/dT at constant pressure.  To obtain the specific electron heat
  // capacity we load the specific electron internal energy (E) and it's first derivative w.r.t
  // temperature.

  // Convert temperatures from keV to degrees Kelvin.
  std::vector<double> vtempsKelvin = vtemperature;
  std::transform(vtemperature.begin(), vtemperature.end(), vtempsKelvin.begin(), keV2K);
  EOS_INTEGER const returnType = EOS_Ue_DT; // ES4enelc;
  std::vector<double> Cve = getF(vdensity, vtempsKelvin, returnType, ETDD_DFDY);
  // Convert back to Temperature units in keV
  std::transform(Cve.begin(), Cve.end(), Cve.begin(), keV2K);
  return Cve;
}

//------------------------------------------------------------------------------------------------//
double Eospac::getSpecificIonInternalEnergy(double temperature, double density) const {
  EOS_INTEGER const returnType = EOS_Uic_DT; // ES4enion;

  // Convert temperatures from keV to degrees Kelvin.
  double vtempsKelvin = keV2K(temperature);
  return getF(dbl_v1(density), dbl_v1(vtempsKelvin), returnType, ETDD_VALUE)[0];
}

//------------------------------------------------------------------------------------------------//
std::vector<double>
Eospac::getSpecificIonInternalEnergy(std::vector<double> const &vtemperature,
                                     std::vector<double> const &vdensity) const {
  // Convert temperatures from keV to degrees Kelvin.
  std::vector<double> vtempsKelvin = vtemperature;
  std::transform(vtemperature.begin(), vtemperature.end(), vtempsKelvin.begin(), keV2K);
  EOS_INTEGER const returnType = EOS_Uic_DT; //ES4enion;
  return getF(vdensity, vtemperature, returnType, ETDD_VALUE);
}

//------------------------------------------------------------------------------------------------//
double Eospac::getIonHeatCapacity(double temperature, double density) const {
  // specific Heat capacity is dE/dT at constant pressure.  To obtain the specific electron heat
  // capacity we load the specific electron internal energy (E) and it's first derivative w.r.t
  // temperature.

  // Convert temperatures from keV to degrees Kelvin.
  double vtempsKelvin = keV2K(temperature);
  EOS_INTEGER const returnType = EOS_Uic_DT; //ES4enion;
  std::vector<double> Cvi = getF(dbl_v1(density), dbl_v1(vtempsKelvin), returnType, ETDD_DFDY);
  // Convert back to Temperature units in keV
  std::transform(Cvi.begin(), Cvi.end(), Cvi.begin(), keV2K);
  return Cvi[0];
}

//------------------------------------------------------------------------------------------------//
std::vector<double> Eospac::getIonHeatCapacity(std::vector<double> const &vtemperature,
                                               std::vector<double> const &vdensity) const {
  // specific Heat capacity is dE/dT at constant pressure.  To obtain the specific electron heat
  // capacity we load the specific electron internal energy (E) and it's first derivative w.r.t
  // temperature.

  // Convert temperatures from keV to degrees Kelvin.
  std::vector<double> vtempsKelvin = vtemperature;
  std::transform(vtemperature.begin(), vtemperature.end(), vtempsKelvin.begin(), keV2K);
  EOS_INTEGER const returnType = EOS_Uic_DT; //ES4enion;
  std::vector<double> Cvi = getF(vdensity, vtempsKelvin, returnType, ETDD_DFDY);
  // Convert back to Temperature units in keV
  std::transform(Cvi.begin(), Cvi.end(), Cvi.begin(), keV2K);
  return Cvi;
}

//------------------------------------------------------------------------------------------------//
double Eospac::getNumFreeElectronsPerIon(double temperature, double density) const {
  // Convert temperatures from keV to degrees Kelvin.
  double vtempsKelvin = keV2K(temperature);
  EOS_INTEGER const returnType = EOS_Zfc_DT; // ES4zfree3; // (zfree3)
  return getF(dbl_v1(density), dbl_v1(vtempsKelvin), returnType, ETDD_VALUE)[0];
}

//------------------------------------------------------------------------------------------------//
std::vector<double> Eospac::getNumFreeElectronsPerIon(std::vector<double> const &vtemperature,
                                                      std::vector<double> const &vdensity) const {
  // Convert temperatures from keV to degrees Kelvin.
  std::vector<double> vtempsKelvin = vtemperature;
  std::transform(vtemperature.begin(), vtemperature.end(), vtempsKelvin.begin(), keV2K);
  EOS_INTEGER const returnType = EOS_Zfc_DT; //ES4zfree3; // (zfree3)
  return getF(vdensity, vtempsKelvin, returnType, ETDD_VALUE);
}

//------------------------------------------------------------------------------------------------//
double Eospac::getElectronThermalConductivity(double temperature, double density) const {
  // Convert temperatures from keV to degrees Kelvin.
  double vtempsKelvin = keV2K(temperature);
  EOS_INTEGER const returnType = EOS_Ktc_DT; //ES4tconde; // (tconde)
  return getF(dbl_v1(density), dbl_v1(vtempsKelvin), returnType, ETDD_VALUE)[0];
}

//------------------------------------------------------------------------------------------------//
std::vector<double>
Eospac::getElectronThermalConductivity(std::vector<double> const &vtemperature,
                                       std::vector<double> const &vdensity) const {
  // Convert temperatures from keV to degrees Kelvin.
  std::vector<double> vtempsKelvin = vtemperature;
  std::transform(vtemperature.begin(), vtemperature.end(), vtempsKelvin.begin(), keV2K);
  EOS_INTEGER const returnType = EOS_Ktc_DT; //ES4tconde; // (tconde)
  return getF(vdensity, vtempsKelvin, returnType, ETDD_VALUE);
}

//------------------------------------------------------------------------------------------------//
double Eospac::getElectronTemperature(     // keV
    double density,                        // g/cm^3
    double SpecificElectronInternalEnergy, // kJ/g
    double /*Tguess*/) const               // keV
{
  EOS_INTEGER const returnType = EOS_T_DUe;
  double Te_K =
      getF(dbl_v1(density), dbl_v1(SpecificElectronInternalEnergy), returnType, ETDD_VALUE)[0];
  return Te_K / keV2K(1.0); // Convert from K back to keV.
}

//------------------------------------------------------------------------------------------------//
double Eospac::getIonTemperature(     // keV
    double density,                   // g/cm^3
    double SpecificIonInternalEnergy, // kJ/g
    double /*Tguess*/) const          // keV
{
  EOS_INTEGER const returnType = EOS_T_DUic;
  // EOS_INTEGER const returnType = EOS_T_DUiz; - I think I need the DUic version!
  double Te_K = getF(dbl_v1(density), dbl_v1(SpecificIonInternalEnergy), returnType, ETDD_VALUE)[0];
  return Te_K / keV2K(1.0); // Convert from K back to keV.
}

// ------- //
// Packing //
// ------- //

//------------------------------------------------------------------------------------------------//
/*!
 * Pack the Eospac state into a char string represented by a vector<char>. This can be used for
 * persistence, communication, etc. by accessing the char * under the vector (required by
 * implication by the standard) with the syntax &char_string[0]. Note, it is unsafe to use iterators
 * because they are \b not required to be char *.
 *
 * Eospac has no state of its own to pack.  However, this function does call the SesameTable pack()
 * because this data is required to rebuild Eospac.
 */
std::vector<char> Eospac::pack() const { return SesTabs.pack(); }

// -------------- //
// Implementation //
// -------------- //

//------------------------------------------------------------------------------------------------//
/*! \brief Retrieves the EoS data associated with the returnType specified and the given (density,
 *         temperature) tuples.
 */
std::vector<double> Eospac::getF(std::vector<double> const &vdensity,
                                 std::vector<double> const &vtemperature,
                                 EOS_INTEGER const returnType,
                                 EosTableDataDerivative const etdd) const {
  // The density and vector parameters must be a tuple.
  Require(vtemperature.size() == vdensity.size());

  unsigned returnTypeTableIndex(tableIndex(returnType));

  // There is one piece of returned information for each (density, temperature) tuple.
  Check(vtemperature.size() < INT32_MAX);
  auto returnSize = static_cast<EOS_INTEGER>(vtemperature.size());

  std::vector<double> returnVals(returnSize);
  std::vector<double> dFx(returnSize);
  std::vector<double> dFy(returnSize);
  EOS_INTEGER errorCode = 0;
  std::vector<double> nc_vx(vdensity);
  std::vector<double> nc_vy(vtemperature);

  eos_Interpolate(&tableHandles[returnTypeTableIndex], &returnSize, &nc_vx[0], &nc_vy[0],
                  &returnVals[0], &dFx[0], &dFy[0], &errorCode);

  if (errorCode != 0) {
    std::ostringstream outputString;
    std::array<EOS_CHAR, EOS_MaxErrMsgLen> errorMessage{};
    eos_GetErrorMessage(&errorCode, errorMessage.data());

    outputString << "\n\tAn unsuccessful request for EOSPAC data was made by eos_Interpolate() "
                 << "from within getF().\n\tThe requested returnType was \"" << returnType
                 << "\" (see eos_Interface.h for type)\n\tThe error code returned was \""
                 << errorCode << "\".\n\tThe associated error message is:\n\t\""
                 << std::string(errorMessage.data()) << "\"\n";

    if (errorCode == EOS_INTERP_EXTRAPOLATED) {
      // If the EOS_INTERP_EXTRAPOLATED error code is returned by either eos_Interpolate or eos_Mix,
      // then the eos_CheckExtrap routine allows the user to determine which (x,y) pairs caused
      // extrapolation and in which direction (high or low), it occurred. The units of the xVals,
      // and yVals arguments listed below are determined by the units listed for each tableType in
      // APPENDIX B and APPENDIX C.

      std::vector<EOS_INTEGER> xyBounds(returnSize);

      eos_CheckExtrap(&tableHandles[returnTypeTableIndex], &returnSize, &nc_vx[0], &nc_vy[0],
                      &xyBounds[0], &errorCode);

      for (int i = 0; i < returnSize; ++i) {
        if (xyBounds[i] == EOS_OK)
          continue;
        outputString << "\tThe specific extrapolation error for entry "
                     << "i = " << i << " is: ";
        if (xyBounds[i] == EOS_xHi_yHi)
          outputString << "\"Both the x and y arguments were high.\"";
        if (xyBounds[i] == EOS_xHi_yOk)
          outputString << "\"The x argument was high.\"";
        if (xyBounds[i] == EOS_xHi_yLo)
          outputString << "\"The x argument was high, the y argument was low.\"";
        if (xyBounds[i] == EOS_xOk_yLo)
          outputString << "\"The y argument was low.\"";
        if (xyBounds[i] == EOS_xLo_yLo)
          outputString << "\"The x argument was low, the y argument was low.\"";
        if (xyBounds[i] == EOS_xLo_yOk)
          outputString << "\"The x argument was low.\"";
        if (xyBounds[i] == EOS_xLo_yHi)
          outputString << "\"The x argument was low, the y argument was high\"";
        if (xyBounds[i] == EOS_xOk_yHi)
          outputString << "\"The y argument was high\"";
        outputString << "\t(x,y) = ( " << nc_vx[0] << ", " << nc_vy[0] << " )\n";
      }
    }

    // This is a fatal exception right now.  It might be useful to throw a specific exception that
    // is derived from EospacException.  The host code could theoretically catch such an exception,
    // fix the problem and then continue.
    throw EospacException(outputString.str());
  }

  switch (etdd) {
  case ETDD_VALUE: {
    break;
  } // return returnVals
  case ETDD_DFDX: {
    returnVals = dFx;
    break;
  }
  case ETDD_DFDY: {
    returnVals = dFy;
    break;
  }
  default: {
    Insist(etdd == ETDD_VALUE || etdd == ETDD_DFDX || etdd == ETDD_DFDY,
           "Bad value for EosTableDataDerivative.");
  }
  }
  return returnVals;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief This member function examines the contents of the data member "SesTabs" and then calls the
 *        EOSPAC routine to load the required EoS Tables.
 */
void Eospac::expandEosTable() const {
  // loop over all possible EOSPAC data types.  If a matid has been assigned to a table then add
  // this information to the vectors returnTypes[] and matIDs[] which are used by EOSPAC.

  std::vector<unsigned> materialTableList(SesTabs.matList());
  for (unsigned mattabid : materialTableList) {
    std::vector<EOS_INTEGER> tableTypes(SesTabs.returnTypes(mattabid));
    for (EOS_INTEGER tableType : tableTypes) {
      matIDs.push_back(mattabid);
      returnTypes.push_back(tableType);
    }
  }

  // Allocate eosTable.  The length and location of eosTable will be modified by es1tabs() as
  // needed.
  for (size_t i = 0; i < returnTypes.size(); ++i)
    tableHandles.push_back(EOS_NullTable);

  // Initialize eosTable and find it's required length

  EOS_INTEGER errorCode(0);
  Check(returnTypes.size() < INT32_MAX);
  auto nTables(static_cast<EOS_INTEGER>(returnTypes.size()));
  eos_CreateTables(&nTables, &returnTypes[0], &matIDs[0], &tableHandles[0], &errorCode);

  // Check for errors
  if (errorCode != EOS_OK) {
    std::ostringstream outputString;
    std::array<EOS_CHAR, EOS_MaxErrMsgLen> errorMessage{};
    eos_GetErrorMessage(&errorCode, errorMessage.data());
    outputString << "\n   An unsuccessful request was made to initialize the EOSPAC table area by "
                 << "expandEosTable().\n  The error code returned by eos_CreateTables(...) was "
                 << "\"" << errorCode << "\".\n  The associated error message is:\n\t\""
                 << std::string(errorMessage.data()) << ".\"\n";
    for (size_t i = 0; i < returnTypes.size(); ++i) {
      EOS_INTEGER tableHandleErrorCode(EOS_OK);
      eos_GetErrorCode(&tableHandles[i], &tableHandleErrorCode);
      eos_GetErrorMessage(&tableHandleErrorCode, errorMessage.data());

      outputString << "\n   The error code associated with tableHandle = " << tableHandles[i]
                   << " was \"" << tableHandleErrorCode << "\".\n   The associated error message "
                   << "is:\n\t\"" << std::string(errorMessage.data()) << "\"\n";
    }

    // Clean up temporaries before we throw the exception.

    // This is a fatal exception right now.  It might be useful to throw a specific exception that
    // is derived from EospacException.  The host code could theoretically catch such an exception,
    // fix the problem and then continue.
    throw EospacException(outputString.str());
  }

  // Load data into table data objects

  eos_LoadTables(&nTables, &tableHandles[0], &errorCode);

  if (errorCode != EOS_OK) {
    std::ostringstream outputString;
    std::array<EOS_CHAR, EOS_MaxErrMsgLen> errorMessage{};
    eos_GetErrorMessage(&errorCode, errorMessage.data());
    outputString << "\n   An unsuccessful request was made to load the EOSPAC table area by "
                 << "expandEosTable().\n  The error code returned by eos_LoadTables(...) was \""
                 << errorCode << "\".\n  The associated error message is:\n\t\""
                 << std::string(errorMessage.data()) << ".\"\n";
    for (size_t i = 0; i < returnTypes.size(); ++i) {
      EOS_INTEGER tableHandleErrorCode(EOS_OK);
      eos_GetErrorCode(&tableHandles[i], &tableHandleErrorCode);
      eos_GetErrorMessage(&tableHandleErrorCode, errorMessage.data());

      outputString << "\n   The error code associated with tableHandle = " << tableHandles[i]
                   << " was \"" << tableHandleErrorCode << "\".\n\tThe associated error message "
                   << "is:\n\t\"" << std::string(errorMessage.data()) << "\"\n";
    }
    throw EospacException(outputString.str());
  }

  // We don't delete eosTable until ~Eospac() is called.
}

//------------------------------------------------------------------------------------------------//
//! Returns true if the EoS data associated with "returnType" has been loaded.
bool Eospac::typeFound(EOS_INTEGER returnType) const {
  return std::find(returnTypes.begin(), returnTypes.end(), returnType) != returnTypes.end();
}

//------------------------------------------------------------------------------------------------//
unsigned Eospac::tableIndex(EOS_INTEGER returnType) const {
  // Loop over all available types.  If the requested type id matches on in the list then return
  // true.  If we reach the end of the list without a match return false.

  // Throw an exception if the required return type has not been loaded by Eospac.
  if (!typeFound(returnType)) {
    std::ostringstream outputString;
    outputString << "\n\tA request was made for data by getF() for which EOSPAC does not have an\n"
                 << "\tassociated material identifier.\n\tRequested returnType = \""
                 << SesTabs.tableName[returnType] << " (" << SesTabs.tableDescription[returnType]
                 << ")\"\n";
    throw EospacUnknownDataType(outputString.str());
  }

  for (unsigned i = 0; i < returnTypes.size(); ++i)
    if (returnType == returnTypes[i])
      return i;

  return 0;
}

//------------------------------------------------------------------------------------------------//
std::vector<EOS_INTEGER> Eospac::initializeInfoItems() {
  std::vector<EOS_INTEGER> ii;
  ii.push_back(EOS_Cmnt_Len);
  ii.push_back(EOS_Exchange_Coeff);
  ii.push_back(EOS_F_Convert_Factor);
  ii.push_back(EOS_Log_Val);
  ii.push_back(EOS_Material_ID);
  ii.push_back(EOS_Mean_Atomic_Mass);
  ii.push_back(EOS_Mean_Atomic_Num);
  ii.push_back(EOS_Modulus);
  ii.push_back(EOS_Normal_Density);
  ii.push_back(EOS_Table_Type);
  ii.push_back(EOS_X_Convert_Factor);
  ii.push_back(EOS_Y_Convert_Factor);

  // These only work for non-inverted tables.
  ii.push_back(EOS_NR);
  ii.push_back(EOS_NT);
  ii.push_back(EOS_Rmin);
  ii.push_back(EOS_Rmax);
  ii.push_back(EOS_Tmin);
  ii.push_back(EOS_Tmax);
  ii.push_back(EOS_Fmin);
  ii.push_back(EOS_Fmax);
  ii.push_back(EOS_NUM_PHASES);
  return ii;
}

//------------------------------------------------------------------------------------------------//
std::vector<std::string> Eospac::initializeInfoItemDescriptions() {
  // These are taken from Appendix E of the EOSPAC user manual.
  using std::string;

  std::vector<std::string> iid;
  iid.emplace_back("The length in chars of the comments for the specified data table");
  iid.emplace_back("The exchange coefficient");
  iid.emplace_back("The conversion factor corresponding to the dependent variable, F(x,y)");
  iid.emplace_back("Non-zero if the data table is in a log10 format");
  iid.emplace_back("The SESAME material identification number");
  iid.emplace_back("The mean atomic mass");
  iid.emplace_back("The mean atomic number");
  iid.emplace_back("The solid bulk modulus");
  iid.emplace_back("The normal density");
  iid.emplace_back("The type of data table. See APPENDIX B and APPENDIX C");
  iid.emplace_back("The conv. factor corresponding to the primary indep. variable, x");
  iid.emplace_back("The conv. factor corresponding to the secondary indep. variable, y");

  // These only work for non-inverted tables.
  iid.emplace_back("The number of densities points");
  iid.emplace_back("The number of temperature points");
  iid.emplace_back("The minimum density value (g/cc)");
  iid.emplace_back("The maximum density value (g/cc)");
  iid.emplace_back("The minimum Temperature value (K)");
  iid.emplace_back("The maximum Temperature value (K)");
  iid.emplace_back("The minimum F value");
  iid.emplace_back("The maximum F value");
  iid.emplace_back("The number of material phases");
  return iid;
}

} // end namespace rtt_cdi_eospac

//------------------------------------------------------------------------------------------------//
// end of Eospac.cc
//------------------------------------------------------------------------------------------------//
