//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_eospac/Eospac.hh
 * \author Kelly Thompson
 * \date   Mon Apr  2 14:14:29 2001
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_eospac_Eospac_hh
#define rtt_cdi_eospac_Eospac_hh

#include "SesameTables.hh"
#include "cdi/EoS.hh"
#include <iostream>

namespace rtt_cdi_eospac {
//================================================================================================//
/*!
 * \class Eospac
 *
 * \brief Provides access to Equation of State data located in Sesame data files.
 *
 * \sa The web page for <a href="http://xweb.lanl.gov/PROJECTS/DATA/">EOSPAC</a>.
 *
 * Eospac allows the client code to retrieve equation of state (EoS) data for a specified material.
 * The material is specified by the SesameTables object which links a lookup table to each type of
 * data requested.
 *
 * This is a concrete class derived from cdi/EoS.  This class allows the client to access
 * (interpolate) on the EoS tables.
 *
 * This class is designed to be used in conjunction with the CDI package. The client code will need
 * to create a SesameTable object that is used in the construction of Eospac.  The Eospac object is
 * then used in the instantiation of a CDI object.  The CDI object might contain other material data
 * (e.g. Opacity data). A single CDI object should only contain information for a single material
 * (the same is true for SesameTable and Eospac objects).
 *
 * <b>User's environment</b>
 *
 * The equation of state data files live in specific locations on the X-Div LAN and ACL.  If you are
 * not working on one of these LANs you must set the following two system environment variables
 * (SESPATHU and SESPATHC) so that the EOSPAC libraries can find the data tables.  On the CCS Linux
 * LAN you can use the following values:
 *
 * export SESPATHU=/ccs/codes/radtran/physical_data/eos
 * export SESPATHC=/ccs/codes/radtran/physical_data/eos
 *
 * Because of the way this object hooks into EOSPAC, we have chosen to implement it as a Meyers
 * singleton.  This ensures that the loaded EOS data remains available until program termination.
 *
 * \example cdi_eospac/test/tEospac.cc
 *
 * This unit test demonstrates the creation of a SesameTable object for aluminum.  Once the Al
 * SesameTable is created the Eospac object for Al is then created using the SesameTable object in
 * the constructor.  The Al Eospac object is then queried for EoS data such as heat capacity, free
 * electrons per ion and a few other things.
 */
//================================================================================================//

class Eospac : public rtt_cdi::EoS {

  // NESTED CLASSES AND TYPEDEFS

  enum EosTableDataDerivative {
    ETDD_VALUE, //!< Return the table value.
    ETDD_DFDX,  //!< Return the first derivative wrt temperature.
    ETDD_DFDY,  //!< Return the first derivative wrt density.
    ETDD_LAST   //!< Last value (invalid)
  };

  // DATA

  // ----------------------- //
  // Specify unique material //
  // ----------------------- //

  /*!
   * \brief The SesameTables object uniquely defines a material.
   *
   * The SesameTables object uniquely defines a material by linking specific lookup tables (sesame,
   * sesou) to material identifiers.
   *
   * \sa rtt_cdi_eospac::SesameTables class definition.
   *
   * \sa Web page for <a href="http://xweb.lanl.gov/projects/data">EOSPAC Data Types</a>
   */
  SesameTables const SesTabs;

  // -------------------- //
  // Available data types //
  // -------------------- //

  // These next four data members are mutable because they specify what data is cached by the Eospac
  // object.  The cached data set may be changed when the user calls a get... function.

  /*!
   * \brief List of material IDs that are specified by SesTabs.
   *
   * \sa returnTypes data member.
   *
   * These are the materials that are available for querying.  There is a one-to-one correspondence
   * between matIDs and returnTypes.  The returnTypes correspond to data that you can request from
   * the sesame tables (e.g.  electron based internal energy has returnType 12) and the
   * corresponding matID value is the material identifier extracted from the associated SesameTables
   * object.
   */
  mutable std::vector<EOS_INTEGER> matIDs;

  /*!
   * \brief List of available EoS data tables that can be queried.
   *
   * \sa matIDs data member.
   *
   * List of numeric identifiers that specify what EoS data tables are available from this
   * object. (e.g. P(T,rho), internal energy, etc.).  There is a one-to-one correspondence between
   * matIDs and returnTypes.  The returnTypes correspond to data that you can request from the
   * sesame tables (e.g.  electron based internal energy has returnType 12) and the corresponding
   * matID value is the material identifier extracted from the associated SesameTables object.
   */
  mutable std::vector<EOS_INTEGER> returnTypes;

  /*! \brief handles to individual portions of the EOS table.
   *
   * The EOS tables are allocated and controlled by EOSPAC.  These handles act as pointers into the
   * table.  Each handle is associated with a tuple {material identifier, data type}.
   */
  mutable std::vector<EOS_INTEGER> tableHandles;

  /*!
   * \brief A list of information enumerations that can be used to query information about EOS
   *        tables.
   */
  mutable std::vector<EOS_INTEGER> infoItems;
  mutable std::vector<std::string> infoItemDescriptions;

public:
  // ------------ //
  // Constructors //
  // ------------ //

  //! The constructor for Eospac.
  explicit Eospac(SesameTables in_SesTabs);

  //! Create an Eospac by unpacking a vector<char> stream.
  explicit Eospac(std::vector<char> const &packed);

  /*!
   * \brief Default Eospac() destructor.
   *
   * This is required to correctly release memory when an Eospac object is destroyed.  We define the
   * destructor in the implementation file to avoid including the unnecessary header files.
   */
  ~Eospac() override;

  //! Disable copy construction
  Eospac(Eospac const &rhs) = delete;

  //! Disable move construction
  Eospac(Eospac &&rhs) noexcept = delete;

  //! Disable assignment
  Eospac &operator=(Eospac const &rhs) = delete;

  //! Disable move-assignment
  Eospac &operator=(Eospac &&rhs) noexcept = delete;

  // --------- //
  // Accessors //
  // --------- //

  //! Print a summary of EOS table information.
  void printTableInformation(EOS_INTEGER const tableType, std::ostream &out = std::cout) const;

  /*!
   * \brief Retrieve the specific electron internal energy given a temperature and a density for
   *        this material.
   *
   * \param density Density of the material in g/cm^3
   * \param temperature Temperature of the material in keV.
   * \return The specific electron internal energy in kJ/g.
   */
  double getSpecificElectronInternalEnergy(double temperature, double density) const override;

  /*!
   * \brief Retrieve a set of specific electron internal energies that correspond to a tuple list of
   *        temperatures and densities for this material.
   *
   * \param[in] vdensity Density of the material in g/cm^3
   * \param[in] vtemperature Temperature of the material in keV.
   * \return The specific electron internal energy in kJ/g.
   */
  std::vector<double>
  getSpecificElectronInternalEnergy(std::vector<double> const &vtemperature,
                                    std::vector<double> const &vdensity) const override;

  /*!
   * \brief Retrieve the electron based heat capacity for this material at the provided density and
   *        temperature.
   *
   * \param density Density of the material in g/cm^3
   * \param temperature Temperature of the material in keV.
   * \return The electron based heat capacity in kJ/g/keV.
   */
  double getElectronHeatCapacity(double temperature, double density) const override;

  /*!
   * \brief Retrieve a set of electron based heat capacities for this material that correspond to
   *        the tuple list of provided densities and temperatures.
   *
   * \param[in] vdensity Density of the material in g/cm^3
   * \param[in] vtemperature Temperature of the material in keV.
   * \return The electron based heat capacity in kJ/g/keV.
   */
  std::vector<double> getElectronHeatCapacity(std::vector<double> const &vtemperature,
                                              std::vector<double> const &vdensity) const override;

  /*!
   * \brief Retrieve the specific ion internal energy for this material at the provided density and
   *        temperature.
   *
   * \param density Density of the material in g/cm^3
   * \param temperature Temperature of the material in keV.
   * \return The specific ion internal energy in kJ/g.
   */
  double getSpecificIonInternalEnergy(double temperature, double density) const override;

  /*!
   * \brief Retrieve a set of specific ion internal energies for this material that correspond to
   *        the tuple list of provided densities and temperatures.
   *
   * \param vdensity Density of the material in g/cm^3
   * \param vtemperature Temperature of the material in keV.
   * \return A vector of specific ion internal energies in kJ/g.
   */
  std::vector<double>
  getSpecificIonInternalEnergy(std::vector<double> const &vtemperature,
                               std::vector<double> const &vdensity) const override;

  /*!
   * \brief Retrieve the ion based heat capacity for this material at the provided density and
   *        temperature.
   *
   * \param density Density of the material in g/cm^3
   * \param temperature Temperature of the material in keV.
   * \return The ion based heat capacity in kJ/g/keV.
   */
  double getIonHeatCapacity(double temperature, double density) const override;

  /*!
   * \brief Retrieve a set of ion based heat capacities for this material that correspond to the
   *        tuple list of provided densities and temperatures.
   *
   * \param vdensity Density of the material in g/cm^3
   * \param vtemperature Temperature of the material in keV.
   * \return A vector of ion based heat capacities in kJ/g/keV.
   */
  std::vector<double> getIonHeatCapacity(std::vector<double> const &vtemperature,
                                         std::vector<double> const &vdensity) const override;

  /*!
   * \brief Retrieve the number of free electrons per ion for this material at the provided density
   *        and temperature.
   *
   * \param density Density of the material in g/cm^3
   * \param temperature Temperature of the material in keV.
   * \return The number of free electrons per ion.
   */
  double getNumFreeElectronsPerIon(double temperature, double density) const override;

  /*!
   * \brief Retrieve a set of free electrons per ion averages for this material that correspond to
   *        the tuple list of provided densities and temperatures.
   *
   * \param vdensity Density of the material in g/cm^3
   * \param vtemperature Temperature of the material in keV.
   * \return A vector of the number of free electrons per ion.
   */
  std::vector<double> getNumFreeElectronsPerIon(std::vector<double> const &vtemperature,
                                                std::vector<double> const &vdensity) const override;

  /*!
   * \brief Retrieve the electron based thermal conductivity for this material at the provided
   *        density and temperature.
   *
   * \param density Density of the material in g/cm^3
   * \param temperature Temperature of the material in keV.
   * \return The electron based thermal conductivity in 1/s/cm.
   */
  double getElectronThermalConductivity(double temperature, double density) const override;

  /*!
   * \brief Retrieve a set of electron based thermal conductivities for this material that
   *        correspond to the tuple list of provided densities and temperatures.
   *
   * \param vdensity Density of the material in g/cm^3
   * \param vtemperature Temperature of the material in keV.
   * \return A vector of electron based thermal conductivities in 1/s/cm.
   */
  std::vector<double>
  getElectronThermalConductivity(std::vector<double> const &vtemperature,
                                 std::vector<double> const &vdensity) const override;

  /*!
   * \brief Retrieve an electron temperature based on the specific electron internal energy.
   *
   * \param density Density of the material in g/cm^3
   * \param SpecificElectronInternalEnergy in kJ/g.
   * \param Tguess Guess of the result to aid the root finder, K.  This is required by the signature
   *           in cdi/EoS.hh but is not used for cdi_eospac.
   * \return temperature Temperature of the material in K.
   */
  double getElectronTemperature(double density, double SpecificElectronInternalEnergy,
                                double Tguess = 1.0) const override;

  /*!
   * \brief Retrieve an ion temperature based on the specific ion internal energy.
   *
   * \param density Density of the material in g/cm^3
   * \param SpecificIonInternalEnergy in kJ/g.
   * \param Tguess Guess of the result to aid the root finder, K.  This is required by the signature
   *           in cdi/EoS.hh but is not used for cdi_eospac.
   * \return temperature Temperature of the material in K.
   */
  double getIonTemperature(double density, double SpecificIonInternalEnergy,
                           double Tguess = 1.0) const override;

  /*!
   * \brief Interface for packing a derived EoS object.
   *
   * Note, the user hands the return value from this function to a derived EoS constructor.  Thus,
   * even though one can pack a EoS through a base class pointer, the client must know the derived
   * type when unpacking.
   */
  std::vector<char> pack() const override;

private:
  // -------------- //
  // Implementation //
  // -------------- //

  /*!
   * \brief Retrieves the EoS data associated with the returnType specified and the given (density,
   *        temperature) tuples.
   *
   * Each of the public access functions calls either getF() or getdFdT() after assigning the
   * correct value to "returnType".
   *
   * \param vdensity A vector of independent values (e.g. temperature or density).
   * \param vtemperature A vector of independent values (e.g. temperature or density).
   * \param returnType The integer index that corresponds to the type of data being retrieved from
   *           the EoS tables.
   * \param etdd Eos Table Derivative
   */
  std::vector<double> getF(std::vector<double> const &vdensity,
                           std::vector<double> const &vtemperature, EOS_INTEGER const returnType,
                           EosTableDataDerivative const etdd) const;

  /*!
   * \brief This member function examines the contents of the data member "SesTabs" and then calls
   *        the EOSPAC routine to load the required EoS Tables.
   */
  void expandEosTable() const;

  //! Returns true if the EoS data associated with "returnType" has been loaded.
  bool typeFound(EOS_INTEGER returnType) const;
  unsigned tableIndex(EOS_INTEGER returnType) const;

  //--------------------//
  // Static Members     //
  //--------------------//

  //! Initialize list of available table info items
  static std::vector<EOS_INTEGER> initializeInfoItems();

  //! Initialize descriptions of available table info items
  static std::vector<std::string> initializeInfoItemDescriptions();

  //! Converts a double to a length one vector.
  static inline std::vector<double> dbl_v1(double const dbl) { return {dbl}; }

  /*!
   * \brief keV2K converts keV temperatures into degrees Kelvin.  libeospac.a requires input
   *        temperatures to use degrees Kelvin.
   *
   * Boltzmann constant k = R/N_A = 8.6174118e-5 eV/K
   *
   * keV2K = 1.1604412e+7 Kelvin/keV
   *
   * This is only used in getF() and getdFdT().
   */
  static inline double keV2K(double tempKeV) {
    const double c = 1.1604412E+7; // Kelvin per keV
    return c * tempKeV;
  }
};

} // end namespace rtt_cdi_eospac

#endif // rtt_cdi_eospac_Eospac_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_eospac/Eospac.hh
//------------------------------------------------------------------------------------------------//
