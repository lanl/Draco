//----------------------------------*-C++-*----------------------------------//
// XYZCoord_sys.hh
// Thomas M. Evans
// Fri Jan 30 16:45:36 1998
/*! 
 * \file   amr_mesh/XYZCoord_sys.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 30 16:52:13 1998
 * \brief  Header file for XYZCoord_sys class library.
 */
//---------------------------------------------------------------------------//
// @> XYZCoord_sys derived class header file
//---------------------------------------------------------------------------//

#ifndef __amr_XYZCoord_sys_hh__
#define __amr_XYZCoord_sys_hh__

//===========================================================================//
// class XYZCoord_sys - 
//
// Purpose : XYZ geometry coordinate system functions, derived
//           class of Coord_sys
//
// revision history:
// -----------------
//  0) original
//  1)  2-18-98 : added getCoord() function for debugging
//  2)  3-12-98 : moved Calc and Set_omega functions to Coord_sys because
//                they are the same for XY and XYZ transport
//  3)  3-16-98 : reserve calc_normal function for later if need be
//  4)   5-5-98 : added sample_pos virtual function
//  5)  6-10-98 : added sample_pos_on_face virtual function
//  6)  6-12-98 : changed interface to sample_pos()
//  7)  4-13-99 : moved to mc package
//  9)  5-21-99 : Modified for AMR mesh topology. Moved to amr_mesh package.
// 
//===========================================================================//

#include "mc/Coord_sys.hh"
#include "rng/Sprng.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <string>
#include <cmath>

namespace rtt_amr 
{
    
using std::vector;
using std::string;
using std::sqrt;

using rtt_rng::Sprng;

/*!
 * \brief  XYZCoord_sys is a base class that is used to define a 
 *         three-dimensional Cartesian coordinate system and provide some 
 *         basic cell sampling member functions needed for implicit monte 
 *         carlo (imc). The derived mc/Coord_sys class inherits functionality 
 *         from the base XYZCoord_sys class.
 *
 *\sa The XYZCoord_sys class is used by both the OS_Mesh and CAR_CU_Mesh 
 *    classes. An \ref amr_overview is provided to describe the basic 
 *    functionality of that particular mesh class (which is also very similar 
 *    to the functionality of the OS_Mesh class from which it was derived).
 */     
class XYZCoord_sys : public Coord_sys
{
  // Begin_Doc xyzcoord_sys-int.tex
  // Begin_Verbatim 

public:
  // default constructor for 3D meshes
/*!
 * \brief Constructs an XYZCoord_sys class object and sets the number of 
 *        spatial dimensions for the derived Coord_sys class object.
 */
    XYZCoord_sys() : Coord_sys(3) {}

  // virtual functions
/*!
 * \brief Returns the coordinate system (i.e, xyz).
 */
    virtual string get_Coord() const { string c = "xyz"; return c; }

/*!
 * \brief Randomly selects a spatial position within a region of space.
 * \param vmin Minimum coordinate values.
 * \param vmax Maximum coordinate values.
 * \param random Random number.
 * \return Spatial position coordinate values. 
 */
    inline virtual vector<double> 
    sample_pos(vector<double> & vmin, vector<double> & vmax, 
	       Sprng & random) const;
    
/*!
 * \brief Randomly selects a spatial position within a region of space with 
 *        a given linear function.
 * \param vmin Minimum coordinate values.
 * \param vmax Maximum coordinate values.
 * \param random Random number.
 * \param slope Linear function gradient.
 * \param center_pt Linear function "intercept" at the cell center-point.
 * \return Spatial position coordinate values. 
 */
    inline virtual vector<double> 
    sample_pos(vector<double> & vmin, vector<double> & vmax, Sprng & random, 
	       vector<double> & slope, double center_pt) const;
    
/*!
 * \brief Randomly selects a spatial position within a planar region of space.
 * \param vmin Minimum coordinate values.
 * \param vmax Maximum coordinate values.
 * \param vmax Cell face (plane) number.
 * \param random Random number.
 * \return Spatial position coordinate values. 
 */
    inline virtual 
    vector<double> sample_pos_on_face(vector<double> & vmin, 
				      vector<double> & vmax, int face, 
				      Sprng & random) const; 
    
  // End_Verbatim 
  // End_Doc 
};

//---------------------------------------------------------------------------//
// INLINE Functions
//---------------------------------------------------------------------------//
// sample the position in an XYZ cell

inline vector<double> 
XYZCoord_sys::sample_pos(vector<double> & min, vector<double> & max,
			 Sprng & random) const
{
  // make return vector
    vector<double> r(3);

  // some assertions
    Check (min.size() == 3);
    Check (max.size() == 3);

    for (int d = 0; d < 3; d++)
    {
      // do uniform sampling
	r[d] = (max[d] - min[d]) * random.ran() + min[d];
    }

  // return assigned array
    return r;
}

//---------------------------------------------------------------------------//
// sample the position in a cell from a linear function

inline vector<double> 
XYZCoord_sys::sample_pos(vector<double> & min, vector<double> & max,
			 Sprng & random, vector<double> & slope, 
			 double center_pt) const
{
  // make return vector
    vector<double> r(3);

  // some assertions
    Check (min.size() == 3);
    Check (max.size() == 3);

    for (int d = 0; d < 3; d++)
    {
      // sample the linear function using linear-linear decomposition of
      // y = mx + b (b is intercept on low side of cell)
	double b = center_pt - slope[d] * (max[d] - min[d]) * 0.5;

      // prob is the fractional area of the negative slope line
	double prob = 0.5 * b / center_pt;

      // sample the dimension
	if (random.ran() <= prob)
	    r[d] = max[d] - (max[d] - min[d]) * sqrt(random.ran());
	else
	    r[d] = min[d] + (max[d] - min[d]) * sqrt(random.ran());
    }

  // return assigned array
    return r;
}

//---------------------------------------------------------------------------//
// sample the position on an XYZ face

inline vector<double> 
XYZCoord_sys::sample_pos_on_face(vector<double> & min, vector<double> & max, 
				 int face, Sprng & random) const
{
  // make return vector
    vector<double> r(3);

  // some assertions
    Check (min.size() == 3);
    Check (max.size() == 3);
    Check (face >= 1 && face <= 6);

    // distribute uniformly over face
    // -z face
    if (face == 1)
    {
	r[0] = (max[0] - min[0]) * random.ran() + min[0];
	r[1] = (max[1] - min[1]) * random.ran() + min[1];
	r[2] = min[2];
    }
    // -y face
    else if (face == 2)
    {
	r[0] = (max[0] - min[0]) * random.ran() + min[0];
	r[1] = min[1];
	r[2] = (max[2] - min[2]) * random.ran() + min[2];
    }
    // -x face
    else if (face == 3)
    {
	r[0] = min[0];
	r[1] = (max[1] - min[1]) * random.ran() + min[1];
	r[2] = (max[2] - min[2]) * random.ran() + min[2];
    }
    // +x face
    else if (face == 4)
    {
	r[0] = max[0];
	r[1] = (max[1] - min[1]) * random.ran() + min[1];
	r[2] = (max[2] - min[2]) * random.ran() + min[2];
    }
    // +y face
    else if (face == 5)
    {
	r[0] = (max[0] - min[0]) * random.ran() + min[0];
	r[1] = max[1];
	r[2] = (max[2] - min[2]) * random.ran() + min[2];
    }
    // +z face
    else if (face == 6)
    {
	r[0] = (max[0] - min[0]) * random.ran() + min[0];
	r[1] = (max[1] - min[1]) * random.ran() + min[1];
	r[2] = max[2];
    }

  // return assigned array
    return r;
}
} // end namespace rtt_amr

#endif                          // __amr_XYZCoord_sys_hh__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/XYZCoord_sys.hh
//---------------------------------------------------------------------------//
