//----------------------------------*-C++-*----------------------------------//
// SMrng.hh
// Thomas M. Evans
// Thu Feb  5 13:34:55 1998
//---------------------------------------------------------------------------//
// @> SMrng class header
//---------------------------------------------------------------------------//

#ifndef __rng_SMrng_hh__
#define __rng_SMrng_hh__

//===========================================================================//
// class SMrng - 
//
// Purpose : Random number class for use with IMC and MC applications;
//           origin, RAN3 function, pg.283 Numerical Recipes in C;
//           subtractive method by Knuth
//
// revision history:
// -----------------
//  0) original
//  1)   4-6-98 : removed explicit from constructor, but added assertion to 
//                assure that the conversion argument to the random number 
//                object is a long int < 0
// 
//===========================================================================//

#include "rng/Names.hh"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cassert>

RNGSPACE

using std::vector;
using std::fill;
using std::ostream;

class SMrng
{
private:
  // generator const variables
    const long mbig;
    const long mseed;
    const int mz;
    const double fac;
  // variables used by Ran()
    int inext, inextp;
    vector<long> ma;
    int iff;
  // seed value, negative integer initializes the sequence
    long idum;
  // counter
    long count;
  // original seed value
    const long seed;

public:
  // must give seed to Random number object
    inline SMrng(long);
 
  // overloaded assignment operators
    inline const SMrng& operator=(const SMrng &rhs);

  // get Random number function
    double ran();

  // determine the number of random numbers used in this object
    long get_count() { return count; }

  // random number diagnostics
    long get_seed() { return seed; }
    void print_values(); 
    double test_avg(int num);
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

// overloaded stream operator
inline ostream& operator<<(ostream &output, SMrng &object)
{
    output << object.ran();
    return output;
}

//---------------------------------------------------------------------------//

inline const SMrng& SMrng::operator=(const SMrng &rhs)
{
  // overloaded assignment operator
    inext  = rhs.inext;
    inextp = rhs.inextp;
    ma     = rhs.ma;
    iff    = rhs.iff;
    idum   = rhs.idum;
    count  = rhs.count;

    return *this;
}

//---------------------------------------------------------------------------//
// inline functions for SMrng
//---------------------------------------------------------------------------//

// constructor
inline SMrng::SMrng(long idum_)
    : mbig(1000000000), mseed(161803398), mz(0), fac(1.0/mbig),
      inext(0), inextp(0), ma(56), iff(0), idum(idum_), count(0),
      seed(idum_)
{     
    assert (idum_ < 0);
    fill(ma.begin(), ma.end(), 0); 
}

CSPACE

#endif                          // __rng_SMrng_hh__

//---------------------------------------------------------------------------//
//                              end of rng/SMrng.hh
//---------------------------------------------------------------------------//
