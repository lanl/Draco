//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Sphere.hh
 * \author Mike Buksas
 * \date   Mon Jun 16 16:14:46 2003
 * \brief  Implements a spherical surface for surface tallies
 * \note   Copyright � 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_mc_Sphere_hh
#define rtt_mc_Sphere_hh

#include "Surface.hh"
#include <vector>

namespace rtt_mc
{

//===========================================================================//
/*!
 * \class Sphere
 *
 * \brief Implement a sphere for surface tracking vis-a-vis the abstract
 * Surface interface
 *
 * Class Sphere implements the abstract interface of class Surface, making it
 * appropiate for surface tracking
 *
 * The spheres are assumed to lie along the z axis. Hence they are described
 * by the z value for the center and their radius.
 *
 * \sa Sphere.cc for detailed descriptions.
 *
 */
/*! 
 * \example mc/test/tstSphere.cc
 * 
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Sphere : public Surface
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS
    
    //! constructor
    Sphere(double center, double radius);

    //! copy constructor
    Sphere(const Sphere &rhs);

    //! destructor
    ~Sphere() { /* ... */ }

    // MANIPULATORS
    
    //! Assignment operator for Sphere
    Sphere& operator=(const Sphere &rhs);

    // ACCESSORS

    double distance_to(std::vector<double> position,
		       const std::vector<double>& direction) const;

    double distance_to(std::vector<double> position,
		       const std::vector<double>& direction,
		       bool is_inside) const;

    bool is_inside(std::vector<double> position) const;
    bool is_inside(std::vector<double> position,
		   const std::vector<double> direction) const;

    double surface_area() const; //!< compute and return surface area of the sphere
    double volume() const;       //!< compute and return volume of the sphere

  private:

    // IMPLEMENTATION

    // DATA

    double center;
    double radius, radius_2;


};

} // end namespace rtt_mc

#endif // rtt_mc_Sphere_hh

//---------------------------------------------------------------------------//
//              end of mc/Sphere.hh
//---------------------------------------------------------------------------//
