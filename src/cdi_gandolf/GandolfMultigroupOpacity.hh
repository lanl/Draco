//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfMultigroupOpacity.hh
 * \author Kelly Thompson
 * \date   Mon Jan 22 13:56:01 2001
 * \brief  GandolfMultigroupOpacity class header file (derived from
 *         cdi/MultigroupOpacity) 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_GandolfMultigroupOpacity_hh__
#define __cdi_gandolf_GandolfMultigroupOpacity_hh__

#include <vector>
#include <string>

#include "ds++/SP.hh"
#include "cdi/MultigroupOpacity.hh"

namespace rtt_cdi_gandolf
{
    // -------------------- //
    // Forward declarations //
    // -------------------- //

    class GandolfFile;
    class GandolfDataTable;

//===========================================================================//
/*!
 * \class GandolfMultigroupOpacity
 *
 * \brief GandolfMultigroupOpacity allows the client code to retrieve
 *        opacity data for a particular material.  Each GandolfOpacity
 *        object represents a specific type of data defined by five
 *        attributes: an IPCRESS File (via a GandolfFile object), a
 *        material identifier, an energy policy (already selecte since
 *        this is a Multigroup Opacity class), a physics model and a
 *        reaction type.
 *
 * \sa  This is a concrete class derived from cdi/MultigroupOpacity.
 *      This class allows to client to access the data in IPCRESS
 *      files via the Gandolf libraries.
 * <p>
 *      This class is designed to be used in conjuction with the CDI.
 *      The client code will create a GandolfMultigroupOpacity object
 *      and use this object as an argument during the CDI
 *      instantiation.  The purpose of this class is to provide a
 *      mechanism for accessing data in IPCRESS files and works by
 *      calling the Gandolf library provided by X-5.  The
 *      GandolfMultigroupOpacity constructor expects four arguments: a
 *      hook to IPCRESS data file (spGandolfFile), a material
 *      identifier, an opacity model (Rosseland or Plank) and an
 *      opacity reaction specifier (total, scattering or absorption).
 *      Once constructed, this object allows the client to access any
 *      data found in the IPCRESS file for that one material.  The
 *      client code will need to create a separate
 *      GandolfMultigroupOpacity object for each material that it
 *      needs information about. Multiple opacity objects can exist
 *      per IPCRESS file.
 * <p> 
 *      This class only provides access to multigroup opacity data.
 *      If the user needs gray opacity IPCRESS data he/she should use
 *      the cdi_gandolf/GandolfGrayOpacity class.
 * <p>
 *      When instantiated, the GandolfMultigroupOpacity object creates
 *      a GandolfDataTable object.  The IPCRESS data is cached in this
 *      table object.  When the client requests an opacity value at a
 *      specified temperature and density the GandolfMultigroupOpcity
 *      object calls the appropriate GANDOLF library routine, which in
 *      turn, interpolates on the data cached in the GandolfDataTable
 *      object.
 */

/*!
 * \example cdi_gandolf/test/tGandolfOpacity.cc
 *
 * \sa  Example of GandolfMultigroupOpacity usage independent of CDI.  In
 *      this example we construct a GandolfMultigroupOpacity object for
 *      the material Aluminum (matID=10001 in our example IPCRESS
 *      file).  We then use the GandolfMultigroupOpacity object to compute
 *      Rosseland Multigroup opacity values for a specified material,
 *      temperature and density.  Other forms of the getOpacity()
 *      accessor are tested along with accessors that return
 *      information about the data set and the cached data table.
 */
//===========================================================================//

class GandolfMultigroupOpacity : public rtt_cdi::MultigroupOpacity
{

    // DATA

    // ----------------------- //
    // Specify unique material //
    // ----------------------- //

    /*!
     * \brief DS++ Smart Pointer to a GandolfFile object.
     *     spGandolfFile acts as a hook to link this object to an
     *     IPCRESS file.
     */
    const rtt_dsxx::SP< GandolfFile > spGandolfFile;

    /*!
     * \brief Identification number for one of the materials found in
     *     the IPCRESS file pointed to by spGandolfFile.
     */
    const int materialID;

    // -------------------- //
    // Available data types //
    // -------------------- //
    
    // The IPCRESS file only holds specific data for each of its materials.

    /*!
     * \brief Number of types of data found in the IPCRESS file.
     */
    int numKeys;

    /*!
     * \brief A list of keys known by the IPCRESS file.
     */
    std::vector< std::string > vKnownKeys;

    // --------------- //
    // Data specifiers //
    // --------------- //

    /*!
     * \brief The physics model that the current data set is based on.
     *        { Rosseland, Plank }.  This enumeration is defined
     *        in cdi/OpacityCommon.hh.
     */
    const rtt_cdi::Model opacityModel;

    /*!
     * \brief The type of reaction rates that the current data set
     *        represents { Total, Scattering, Absorption }. This
     *        enumeration is defined in cdi/OpacityCommon.hh.
     */
    const rtt_cdi::Reaction opacityReaction;

    /*!
     * \brief A string that identifies the energy policy for this
     *     class.
     */
    const std::string energyPolicyDescriptor;

    // -------------------- //
    // Opacity lookup table //
    // -------------------- //

    /*!
     * \brief spGandolfDataTable contains a cached copy of the
     *        requested IPCRESS opacity lookup table.
     *
     * There is a one-to-one relationship between GandolfOpacity and
     * GandolfDataTable. 
     */
    rtt_dsxx::SP< GandolfDataTable > spGandolfDataTable;

  public:

    // ------------ //
    // Constructors //
    // ------------ //

    /*!
     * \brief This is the default GandolfMultigroupOpacity
     *     constructor.  It requires four arguments plus the energy
     *     policy (this class) to be instantiated.
     * 
     * \sa The combiniation of a data file and a material ID uniquely 
     *     specifies a material.  If we add the Model, Reaction and
     *     EnergyPolicy the opacity table is uniquely defined.
     *
     * \parameter _spGandolfFile This smart pointer links an IPCRESS
     *     file (via the GandolfFile object) to a GandolfOpacity
     *     object. There may be many GandolfOpacity objects per
     *     GandolfFile object but only one GandolfFile object for each 
     *     GandolfOpacity object.
     *
     * \parameter _materialID An identifier that links the
     *     GandolfOpacity object to a single material found in the
     *     specified IPCRESS file.
     *
     * \parameter _opacityModel The physics model that the current
     *     data set is based on.
     *
     * \parameter _opacityReaction The type of reaction rate that the
     *     current data set represents. 
     */
    GandolfMultigroupOpacity( const rtt_dsxx::SP< GandolfFile > _spGandolfFile,
			      const int _materialID, 
			      const rtt_cdi::Model _opacityModel,
			      const rtt_cdi::Reaction _opacityReaction );

    /*!
     * \brief Default GandolfOpacity() destructor.
     *
     * \sa This is required to correctly release memory when a 
     *     GandolfMultigroupOpacity is destroyed.  This constructor's
     *     * definition must be declared in the implementation file so
     *     that * we can avoid including too many header files
     */
    ~GandolfMultigroupOpacity();

    // --------- //
    // Accessors //
    // --------- //

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of (temperature,density) tuples.
     *     An opacity value will be returned for each tuple.  The
     *     temperatureIterator and density iterators are required to
     *     be the same length.  The opacity iterator should also have
     *     this same length for gray data or this length times the
     *     number of energy groups for multigroup data.
     *
     * The InputIterator and OutputIterator classes must be
     *     instantiated for each STL container used.  This has already
     *     been done for a few STL containers in GandolfOpacity_pt.cc.
     * 
     * \parameter temperatureIterator The beginning position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter temperatureIteratorEnd The end position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter densityIterator The beginning position of a STL
     *     container that holds a list of densities.
     *
     * \parameter densityIteratorEnd The end position of a STL
     *     container that holds a list of temperatures.
     * 
     * \parameter opacityIterator The beginning position of a STL
     *     container into which opacity values corresponding to the
     *     given (temperature,density) tuple will be stored.
     *
     * \return A list (of type OutputIterator) of opacities are
     *     returned.  These opacities correspond to the temperature
     *     and density values provied in the two InputIterators.
     */
//     template < class InputIterator, class OutputIterator >
//     OutputIterator getOpacity( InputIterator temperatureIterator, 
// 			       InputIterator temperatureIteratorEnd,
// 			       InputIterator densityIterator, 
// 			       InputIterator densityIteratorEnd,
// 			       OutputIterator opacityIterator ) const;

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of temperatures in an STL container.
     *     An opacity value will be returned for each temperature
     *     provided.  The opacity iterator should be the same length
     *     as the temperatureIterator for gray data or the length of
     *     the temperatureIterator times the number of energy groups
     *     for multigroup data.
     *
     * The InputIterator and OutputIterator classes must be
     *     instantiated for each STL container used.  This has already
     *     been done for a few STL containers in GandolfOpacity_pt.cc.
     *
     * \parameter temperatureIterator The beginning position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter temperatureIteratorEnd The end position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter targetDensity The single density value used when
     *     computing opacities for each given temperature.
     * 
     * \parameter opacityIterator The beginning position of a STL
     *     container into which opacity values (corresponding to the
     *     provided temperature and density values) will be stored.
     *
     * \return A list (of type OutputIterator) of opacities are
     *     returned.  These opacities correspond to the temperature
     *     and density values provied.
     */
//     template < class InputIterator, class OutputIterator >
//     OutputIterator getOpacity( InputIterator temperatureIterator,
// 			       InputIterator temperatureIteratorEnd,
// 			       const double targetDensity,
// 			       OutputIterator opacityIterator ) const;

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of densities in an STL container.
     *     An opacity value will be returned for each density
     *     provided.  The opacity iterator should be the same length
     *     as the densityIterator for gray data or the length of the
     *     densityIterator times the number of energy groups for
     *     multigroup data.
     *
     * The InputIterator and OutputIterator classes must be
     *     instantiated for each STL container used.  This has already
     *     been done for a few STL containers in GandolfOpacity_pt.cc.
     *
     * \parameter targetTemperature The single temperature value used when
     *     computing opacities for each given density.
     * 
     * \parameter densityIterator The beginning position of a STL
     *     container that holds a list of densities.
     *
     * \parameter densityIteratorEnd The end position of a STL
     *     container that holds a list of densities.
     *
     * \parameter opacityIterator The beginning position of a STL
     *     container into which opacity values (corresponding to the
     *     provided temperature and density values) will be stored.
     *
     * \return A list (of type OutputIterator) of opacities are
     *     returned.  These opacities correspond to the temperature
     *     and density values provied.
     */
//     template < class InputIterator, class OutputIterator >
//     OutputIterator getOpacity( const double targetTemperature,
// 			       InputIterator densityIterator, 
// 			       InputIterator densityIteratorEnd,
// 			       OutputIterator opacityIterator ) const;

    /*!
     * \brief Opacity accessor that returns a single opacity (or a
     *     vector of opacities for the multigroup EnergyPolicy) that 
     *     corresponds to the provided temperature and density.
     *
     * \parameter targetTemperature The temperature value for which an
     *     opacity value is being requested.
     *
     * \parameter targetDensity The density value for which an opacity 
     *     value is being requested.
     *
     * \return A single opacity (or a vector of opacities for the
     *     multigroup EnergyPolicy).
     */
    std::vector< double > getOpacity( const double targetTemperature,
				      const double targetDensity ) const; 
    
    /*!
     * \brief Opacity accessor that returns a vector of opacities (or a
     *     vector of vectors of opacities for the multigroup
     *     EnergyPolicy) that correspond to the provided vector of
     *     temperatures and a single density value.
     *
     * \parameter targetTemperature A vector of temperature values for
     *     which opacity values are being requested.
     *
     * \parameter targetDensity The density value for which an opacity 
     *     value is being requested.
     *
     * \return A vector of opacities (or a vector of vectors of
     *     opacities for the multigroup EnergyPolicy).
     */
    std::vector< std::vector< double > > getOpacity( 
	const std::vector<double>& targetTemperature,
	const double targetDensity ) const; 

    /*!
     * \brief Opacity accessor that returns a vector of opacities (or a
     *     vector of vectors of opacities for the multigroup
     *     EnergyPolicy) that correspond to the provided vector of
     *     densities and a single temperature value.
     *
     * \parameter targetTemperature The temperature value for which an 
     *     opacity value is being requested.
     *
     * \parameter targetDensity A vector of density values for which
     *     opacity values are being requested.
     *
     * \return A vector of opacities (or a vector of vectors of
     *     opacities for the multigroup EnergyPolicy).
     */
    std::vector< std::vector< double > > getOpacity( 
	const double targetTemperature,
	const std::vector<double>& targetDensity ) const; 

    // It is not clear how to assume order of opacity(temp,dens) when
    // accessed in this manner --> for now use the STL-style accessor
    // or a loop over one of the other vector-accessors.

//     std::vector< std::vector< double > > getOpacity( 
// 	const std::vector<double>& targetTemperature,
// 	const std::vector<double>& targetDensity ) const;

    /*!
     * \brief Returns a string that describes the templated
     *     EnergyPolicy.  Currently this will return either "mg" or
     *     "gray."
     */ 
    const std::string& getEnergyPolicyDescriptor() const {
	return energyPolicyDescriptor; };

    /*!
     * \brief Returns a "plain English" description of the opacity
     *     data that this class references. (e.g. "Multigroup Rosseland
     *     Scattering".) 
     *
     * The definition of this function is not included here to prevent 
     *     the inclusion of the GandolfFile.hh definitions within this 
     *     header file.
     */
    const std::string& getDataDescriptor() const;

    /*!
     * \brief Returns the name of the associated IPCRESS file.
     *
     * The definition of this function is not included here to prevent 
     *     the inclusion of the GandolfFile.hh definitions within this 
     *     header file.
     */
    const std::string& getDataFilename() const;

    /*!
     * \brief Returns a vector of temperatures that define the cached
     *     opacity data table.
     * 
     * We do not return a const reference because this function
     * must construct this information from more fundamental tables.
     */
    const std::vector< double >& getTemperatureGrid() const;

    /*!
     * \brief Returns a vector of densities that define the cached
     *     opacity data table.
     * 
     * We do not return a const reference because this function
     * must construct this information from more fundamental tables.
     */
    const std::vector< double >& getDensityGrid() const;

    /*!
     * \brief Returns a vector of energy values (keV) that define the
     *     energy boundaries of the cached multigroup opacity data
     *     table.  (This accessor is only valid for the Multigroup
     *     EnergyPolicy version of GandolfOpacity.)
     */
    const std::vector< double >& getGroupBoundaries() const;
    
    /*!
     * \brief Returns the size of the temperature grid.
     */
    int getNumTemperatures() const;

    /*! 
     * \brief Returns the size of the density grid.
     */
    int getNumDensities() const;

    /*!
     * \brief Returns the number of group boundaries found in the
     *     current multigroup data set.
     */
    int getNumGroupBoundaries() const;

    /*!
     * \brief Returns the number of gruops found in the current
     *     multigroup data set.
     */
    int getNumGroups() const {
	return getNumGroupBoundaries() - 1; };

}; // end of class GandolfMultigroupOpacity

} // end namespace rtt_cdi_gandolf

#endif // __cdi_gandolf_GandolfMultigroupOpacity_hh__

//---------------------------------------------------------------------------//
//                end of cdi_gandolf/GandolfMultigroupOpacity.hh
//---------------------------------------------------------------------------//
