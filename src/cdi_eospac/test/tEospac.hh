//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/test/tEospac.hh
 * \author Kelly Thompson
 * \date   Mon Apr 2 14:15:57 2001
 * \brief  Header file for tEospac.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_eospac_tEospac_hh__
#define __cdi_eospac_tEospac_hh__

#include "UnitTestFrame/TestAppNoC4.hh"

namespace rtt_cdi_eospac_test
{

//===========================================================================//
/*!
 * \class tEospac
 *
 * \brief 
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class tEospac : public rtt_UnitTestFrame::TestApp
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    tEospac( int argc, char *argv[], std::ostream &os_in );
    //Defaulted: tCDI(const tCDI &rhs);
    //Defaulted: ~tCDI();

    // MANIPULATORS
    
    //Defaulted: tCDI& operator=(const tCDI &rhs);

    // ACCESSORS
    std::string name() const { return "tGandolfOpacity"; }
    std::string version() const;

  protected:
    
    std::string runTest();

  private:
    
    // IMPLEMENTATION
    
    /*!
     * \brief Returns true if the two values are identical to 10
     *        decimal places. 
     */
//     bool match( const double computedValue,
// 		const double referenceValue ) const;

    /*!
     * \brief Returns true if the elements of the two vectors are
     *        identical to 10 decimal places.
     */
//     bool match( const std::vector< double >& computedValue, 
// 		const std::vector< double >& refereneceValue ) const;

    /*!
     * \brief Returns true if the two values are identical to 10
     *        decimal places. 
     */
//     bool match( const std::vector< std::vector< double > >& computedValue,
//  		const std::vector< std::vector< double > >& referenceValue ) const;

    
};

} // end namespace rtt_cdi_eospac_test

#endif // __cdi_eospac_tEospac_hh__

//---------------------------------------------------------------------------//
// end of cdi_eospac/test/tEospac.hh
//---------------------------------------------------------------------------//

