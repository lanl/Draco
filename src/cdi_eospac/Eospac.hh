//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/Eospac.hh
 * \author Kelly Thompson
 * \date   Mon Apr  2 14:14:29 2001
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_eospac_Eospac_hh__
#define __cdi_eospac_Eospac_hh__

// cdi_eospac dependencies
#include "SesameTables.hh"
#include <cdi_eospac/config.h>  // defines V_FLOAT - also see configure.in

// Draco dependencies
#include "cdi/EoS.hh"

// C++ standard library dependencies
#include <vector>

namespace rtt_cdi_eospac
{
    //===========================================================================//
    /*!
     * \class Eospac
     *
     * \brief Provides access to Equation of State data located in
     *        Sesame data files.
     *
     * \sa The web page for <a 
     *     href="http://laurel.lanl.gov/XCI/PROJECTS/DATA/eos/eos.html">EOSPAC</a>.
     *
     * \sa The web page for <a 
     *     href="http://int.lanl.gov/projects/sdm/win/materials/">Eos
     *     Material Identifiers</a>.  This web site also does dynamic
     *     plotting of EoS values.
     *
     * \sa The web page for <a href="http://laurel.lanl.gov/PROJECTS/DATA/eos/sesameLocations.html">data table locations</a>. 
     *
     * Eospac allows the client code to retrive equation of state (EoS)
     * data for a specified material.  The material is specified by
     * the SesameTables object which links a lookup table to each type 
     * of data requested.
     *
     * This is a concrete class derived from cdi/EoS.  This class
     * allows the client to access (interpolate) on the EoS tables
     * provided by X-5 (sesame, sescu1, sesou, sescu and sescu9).
     *
     * This class is designed to be used in conjuction with the CDI
     * package. The client code will need to create a SesameTable
     * object that is used in the construction of Eospac.  The Eospac
     * object is then used in the instantiation of a CDI object.  The
     * CDI object might contain other material data (e.g. Opacity
     * data). A single CDI object should only contain information for
     * a single material (the same is true for SesameTable and Eospac
     * objects). 
     *
     * When compiling DRACO with support for this package you must add 
     * the following option on the configure line:<br><br>
     * <tt>   --with-eospac-lib=${VENDORS}/eospac/IRIX64/lib64</tt><br><br>
     * - or -
     * <tt>   --with-eospac-lib=${VENDORS}/eospac/Linux/lib</tt><br><br>
     */

    /*!
     * \example cdi_eospac/test/tEospac.cc
     * 
     * This unit test demonstrates the creation of a SesameTable
     * object for aluminum.  Once the Al SesameTable is created the
     * Eospac object for Al is then created using the SesameTable
     * object in the constructor.  The Al Eospac object is then
     * queried for EoS data such as heat capacity, free electrons per
     * ion and a few other things.
     */

    // Todo:
    // --------------------
    // 1. Add STL like accessors.

    // revision history:
    // -----------------
    // 0) original

    //===========================================================================//
    
    class Eospac : public rtt_cdi::EoS
    {
	
	// NESTED CLASSES AND TYPEDEFS
	
	// DATA

	// ----------------------- //
	// Specify unique material //
	// ----------------------- //
	
	/*!
	 * \brief The SesameTables object uniquely defines a material.
	 *
	 * The SesameTables object uniquely defines a material by
	 * linking specific lookup tables (sesame, sescu1, sesou,
	 * sescu and sescu9) to material identifiers. 
	 *
	 * \sa rtt_cdi_eospac::SesameTables class definition.
	 *
	 * \sa Web page for <a
	 * href="http://laurel.lanl.gov/XCI/PROJECTS/DATA/eos/UsersDocument/HTML/EOSPAC.html#5.4">EOSPAC
	 * Data Types</a> and the web page for <a 
	 * href="http://int.lanl.gov/projects/sdm/win/materials/">EoS
	 * material and table identifiers</a>.
	 */
	const SesameTables SesTabs;

	/*!
	 * \brief The number of material regions.  
	 *
	 * For a single Eospac object there is always exactly one
	 * region.       
	 */
	const int numRegions;

	/*!
	 * \brief regionIndex is uniquely "1" because there will only
	 *        be one material region per eospac object.
	 */
	const int regionIndex;

	// -------------------- //
	// Available data types //
	// -------------------- //

	// These next four data members are mutalbe because they
	// specify what data is cached by the Eospac object.  The
	// cached data set may be changed when the user calls a
	// get... function.
	
	/*!
	 * \brief List of materierial IDs that are specified by
	 *        SesTabs.
	 *
	 * \sa returnTypes data member.
	 *
	 * These are the materials that are available for querying.
	 * There is a one-to-one correspondence between matIDs and
	 * returnTypes.  The returnTypes correspond to data that you
	 * can request from the sesame tables (e.g. electron based
	 * interal energy has returnType 12) and the corresponding
	 * matID value is the material identifier extracted from the
	 * associated SesameTables object.
	 */
	mutable std::vector< int > matIDs;

	/*!
	 * \brief List of available EoS data tables that can be
	 *        queried.
	 *
	 * \sa matIDs data member.
	 *
	 * List of numeric identifiers that specify what EoS data
	 * tables are available from this object. (e.g. P(T,rho),
	 * internal energy, etc.).  There is a one-to-one
	 * correspondence between matIDs and returnTypes.  The
	 * returnTypes correspond to data that you can request from
	 * the sesame tables (e.g. electron based interal energy has
	 * returnType 12) and the corresponding matID value is the
	 * material identifier extracted from the associated
	 * SesameTables object.
	 */
// 	mutable std::vector< int > returnTypes;
	mutable std::vector< ES4DataType > returnTypes;

	/*!
	 * \brief The eos tables are cached with this pointer. 
	 *
	 * The length and contents are set in the contructor by the
	 * EOSPAC routine es1tabs_(). 
	 *
	 * EOSPAC uses different data types on different
	 * architectures.  Most of these differences are dealt with
	 * directly by the EospacWrapper class.  Unfortunatley, the
	 * actual eosTable must be controlled by the host (Eospac) and 
	 * must match the data type expected by libeospac.a.
	 */
	mutable V_FLOAT *eosTable;

	/*!
	 * \brief The length (in words) of the eosTable.
	 *
	 * \sa eosTable data member.
	 */
	mutable int eosTableLength;
	
      public:
	
	// ------------ //
	// Constructors //
	// ------------ //

	/*!
	 * \brief The constructor for Eospac.
	 *
	 * \sa The definition of rtt_cdi_eospac::SesameTables.
	 *
	 * \param SesTabs A rtt_cdi_eospac::SesameTables object that
	 * defines what data tables will be available for queries from
	 * the Eospac object. 
	 */
	Eospac( const SesameTables& SesTabs );

	// (defaulted) Eospac(const Eospac &rhs);

	/*!
	 * \brief Default Eospac() destructor.
	 *
	 * This is required to correctly release memeroyt when an
	 * Eospac object is destroyed.  We define the destructor in
	 * the implementation file to avoid including the unnecessary
	 * header files.
	 */
	~Eospac();
	
	// MANIPULATORS
	
	// (defaulted ) Eospac& operator=(const Eospac &rhs);
	
	// --------- //
	// Accessors //
	// --------- //

	/*!
	 * \brief Retrieve the specific electron internal energy given 
	 *        a temperature and a density for this material.
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in keV.
	 * \return The specific electron internal energy in kJ/g.
	 */
	double getSpecificElectronInternalEnergy(
	    double temperature, double density ) const;

	/*!
	 * \brief Retrieve a set of specific electron internal
	 *        energies that correspond to a tuple list of
	 *        temperatures and densities for this material.
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in keV.
	 * \return The specific electron internal energy in kJ/g.
	 */
	std::vector< double > getSpecificElectronInternalEnergy(
	    const std::vector< double >& vtemperature, 
	    const std::vector< double >& vdensity ) const;
	    
	/*!
	 * \brief Retrieve the electron based heat capacity for this
	 *        material at the provided density and temperature.
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in keV.
	 * \return The electron based heat capacity in kJ/g/keV.
	 */
	double getElectronHeatCapacity(
	    double temperature, double density ) const;

	/*!
	 * \brief Retrieve a set of electron based heat capacities for
	 *        this material that correspond to the tuple list of
	 *        provided densities and temperatures. 
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in keV.
	 * \return The electron based heat capacity in kJ/g/keV.
	 */
	std::vector< double > getElectronHeatCapacity(
	    const std::vector< double >& vtemperature, 
	    const std::vector< double >& vdensity ) const;

	/*!
	 * \brief Retrieve the specific ion internal energy for this
	 *        material at the provided density and temperature.
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in keV.
	 * \return The specific ion internal energy in kJ/g.
	 */
	double getSpecificIonInternalEnergy(
	    double temperature, double density ) const;

	/*!
	 * \brief Retrieve a set of specific ion internal energies for
	 *        this material that correspond to the tuple list of
	 *        provided densities and temperatures. 
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in keV.
	 * \return A vector of specific ion internal energies in kJ/g.
	 */
	std::vector< double > getSpecificIonInternalEnergy(
	    const std::vector< double >& vtemperature, 
	    const std::vector< double >& vdensity ) const;

	/*!
	 * \brief Retrieve the ion based heat capacity for this
	 *        material at the provided density and temperature.
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in keV.
	 * \return The ion based heat capacity in kJ/g/keV.
	 */
	double getIonHeatCapacity(
	    double temperature, double density ) const;

	/*!
	 * \brief Retrieve a set of ion based heat capacities for
	 *        this material that correspond to the tuple list of
	 *        provided densities and temperatures. 
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in keV.
	 * \return A vector of ion based heat capacities in kJ/g/keV.
	 */
	std::vector< double > getIonHeatCapacity(
	    const std::vector< double >& vtemperature,
	    const std::vector< double >& vdensity ) const;

	/*!
	 * \brief Retrieve the number of free electrons per ion for this
	 *        material at the provided density and temperature.
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in keV.
	 * \return The number of free electrons per ion.
	 */
	double getNumFreeElectronsPerIon(
	    double temperature, double density ) const;

	/*!
	 * \brief Retrieve a set of free electrons per ion averages for
	 *        this material that correspond to the tuple list of
	 *        provided densities and temperatures. 
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in keV.
	 * \return A vector of the number of free electrons per ion.
	 */
	std::vector< double > getNumFreeElectronsPerIon(
	    const std::vector< double >& vtemperature,
	    const std::vector< double >& vdensity ) const;

	/*!
	 * \brief Retrieve the electron based thermal conductivity for this
	 *        material at the provided density and temperature.
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in keV.
	 * \return The electron based thermal conductivity in 1/s/cm.
	 */
	double getElectronThermalConductivity(
	    double temperature, double density ) const;

	/*!
	 * \brief Retrieve a set of electron based thermal conductivities for
	 *        this material that correspond to the tuple list of
	 *        provided densities and temperatures. 
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in keV.
	 * \return A vector of electron based thermal conductivities
	 * in 1/s/cm.
	 */
	std::vector< double > getElectronThermalConductivity(
	    const std::vector< double >& vtemperature,
	    const std::vector< double >& vdensity ) const;

	/*!
	 * \brief Interface for packing a derived EoS object.
	 *
	 * Note, the user hands the return value from this function to a
	 * derived EoS constructor.  Thus, even though one can pack a EoS
	 * through a base class pointer, the client must know the derived
	 * type when unpacking.
	 */
	std::vector<char> pack() const;

      private:
	
	// -------------- //
	// Implementation //
	// -------------- //

	/*!
	 * \brief Retrieves the EoS data associated with the returnType 
	 *        specified and the given (density, temperature) tuples.
	 *
	 * Each of the public access functions calls either getF() or
	 * getdFdT() after assigning the correct value to
	 * "returnType".
	 *
	 * \param vdensity A vector of density values (g/cm^3).
	 * \param vtemperature A vector of temperature values (K).
	 * \param returnType The integer index that corresponds to the 
	 *        type of data being retrieved from the EoS tables.
	 */
	std::vector< double > getF( 
	    const std::vector< double >& vdensity, 
	    const std::vector< double >& vtemperature, 
	    ES4DataType returnType ) const;
	    
	/*!
	 * \brief Retrieves the EoS data associated with the returnType 
	 *        specified and the given (density, temperature) tuples.
	 *
	 * Each of the public access functions calls either getF() or
	 * getdFdT() after assigning the correct value to
	 * "returnType".
	 *
	 * \param vdensity A vector of density values (g/cm^3).
	 * \param vtemperature A vector of temperature values (K).
	 * \param returnType The integer index that corresponds to the 
	 *        type of data being retrieved from the EoS tables.
	 */
	std::vector< double > getdFdT( 
	    const std::vector< double >& vdensity, 
	    const std::vector< double >& vtemperature, 
	    ES4DataType returnType ) const;

	/*!
	 * \brief This member function examines the contents of
	 *        the data member "SesTabs" and then calls the EOSPAC
	 *        routine to load the required EoS Tables.
	 */
 	void expandEosTable () const;

	/*!
	 * \brief Returns true if the EoS data associated with
	 *        "returnType" has been loaded.
	 */
	bool typeFound( ES4DataType returnType ) const;
	
	/*!
	 * \brief Converts a double to a length one vector.
	 */
	std::vector< double > dbl_v1( const double dbl ) const;

	/*!
	 * \brief keV2K converts keV temperatures into degrees
	 *        Kelvin.  libeospac.a requires input temperatures to
	 *        use degrees Kelvin.
	 *
	 * Boltzmann constant k = R/N_A = 8.6174118e-5 eV/K
	 *
	 * keV2K = 1.1604412e+7 Kelvin/keV
	 *
	 * This is only used in getF() and getdFdT().
	 */
	static inline double keV2K( double tempKeV )
	{
	    const double c = 1.1604412E+7; // Kelven per keV
	    return c*tempKeV;
	}
    };
    
} // end namespace rtt_cdi_eospac

#endif // __cdi_eospac_Eospac_hh__

//---------------------------------------------------------------------------//
// end of cdi_eospac/Eospac.hh
//---------------------------------------------------------------------------//
