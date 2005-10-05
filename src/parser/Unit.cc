//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file Unit.cc
 * \author Kent G. Budge
 * \date Wed Jan 22 15:18:23 MST 2003
 * \brief Definitions of Unit methods.
 * \note   Copyright � 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include "Unit.hh"

namespace
{
using namespace std;

//---------------------------------------------------------------------------//
/*! 
 * \brief Write text for a component of a unit.
 * 
 * This is a helper function for operator<<(ostream, const Unit &).  It has
 * the logic for ensuring that dashes are inserted where necessary, e.g., the
 * dash in kg-m^2/s^2.
 *
 * \param str Stream to which to write the unit component.
 * \param dash A dash is required after the last component written to the
 * stream.  May be modified on return.
 * \param value Unit dimension value
 * \param unit Unit name
 */

void dash_insert(ostream &str,
		 bool &dash,
		 double const value,
		 char const *const name)
{
    if (value!=0.0)
    {
	if (dash)
	{
	    str << '-';
	}
	if (value!=1.0)
	{
	    str << name << '^' << value;
	}
	else
	{
	    str << name;
	}
	dash = true;
    }
}

}

namespace rtt_parser 
{
using namespace std;

//---------------------------------------------------------------------------//
/*! 
 * \param str Stream to which to write the text description.
 * \param u Unit to write the text description for.
 * \return A reference to s.
 */

std::ostream &operator<<(std::ostream &str, const Unit &u)
{
    str << u.conv << ' ';
    bool dash = false;

    dash_insert(str, dash, u.m, "m");
    dash_insert(str, dash, u.kg, "kg");
    dash_insert(str, dash, u.s, "s");
    dash_insert(str, dash, u.A, "A");
    dash_insert(str, dash, u.K, "K");
    dash_insert(str, dash, u.mol, "mol");
    dash_insert(str, dash, u.cd, "cd");
    dash_insert(str, dash, u.rad, "rad");
    dash_insert(str, dash, u.sr, "sr");
    return str;
}

} // namespace rtt_parser

//---------------------------------------------------------------------------//
//                      end of Unit.cc
//---------------------------------------------------------------------------//
