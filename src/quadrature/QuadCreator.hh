//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/QuadCreator.hh
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:46:17 2000
 * \brief  Quadrature Creator class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __quadrature_QuadCreator_hh__
#define __quadrature_QuadCreator_hh__

// Quadrature.hh defines the "Quadrature" class.
#include "Quadrature.hh"

namespace rtt_quadrature
{

static double PI = 3.14159265358979323846;
 
//===========================================================================//
/*!
 * \class QuadCreator
 *
 * \brief A class to instantiate Quadrature objects.
 *
 * The generation of Quadrature objects can be implemented as a
 * "Parameterized Factory Method" (Booch).  This class acts as the virtual
 * creator for all quadrature objects.
 *
 * \example quadrature/test/it_quad.cc
 * 
 * Example of Quadrature usage.  In this example the test code requests input 
 * from the user via stdin, creates a quadrature set, displays some
 * information about the set and then destroys itself.
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class QuadCreator 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

  public:

    // DATA

    // Identifiers for each type of quadrature that may be generated by this
    // creator. This list will expand as more types of quadrature are added
    // to this package.
    enum Qid { GaussLeg, LevelSym };

    // CREATORS

    /*!
     * \brief QuadCreate constructs a Quadrature object.
     *
     * I'm not sure if this needs to be virtual or not.
     * I'm not sure if we need the assignment, copy and destructor member
     * functions. 
     */
    virtual Quadrature* QuadCreate( Qid quad_type, int sn_order = 4, 
				    double norm = 4*PI );
    //    QuadCreator(const QuadCreator &rhs);
    //   ~QuadCreator();

    // MANIPULATORS
    
    //    QuadCreator& operator=(const QuadCreator &rhs);

    // ACCESSORS

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_quadrature

#endif                          // __quadrature_QuadCreator_hh__

//---------------------------------------------------------------------------//
//                              end of quadrature/QuadCreator.hh
//---------------------------------------------------------------------------//
