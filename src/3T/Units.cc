//----------------------------------*-C++-*----------------------------------//
// Units.cc
// Randy M. Roberts
// Tue Mar 17 14:52:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/Units.hh"
#include <limits>

double Units::minConversion()
{
    static double minConv = std::numeric_limits<double>::min();
    return minConv;
}

//---------------------------------------------------------------------------//
//                              end of Units.cc
//---------------------------------------------------------------------------//
