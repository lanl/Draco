//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   roots/test/tstlnsrch.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:39:20 2004
 * \brief  
 * \note   Copyright 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>

#include "ds++/ScalarUnitTest.hh"
#include "c4/global.hh"
#include "../Release.hh"
#include "linear/fnorm.hh"
#include "../lnsrch.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;
using namespace rtt_roots;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void func(const vector<double> &x, vector<double> &fvec)
{
    fvec.resize(2);
    fvec[0] = 7.2*x[0] + 3.5*x[1] + 2.3;
    fvec[1] = -2.2*x[0] + 2.7*x[1] + 5.4;
}

//---------------------------------------------------------------------------//
void tstlnsrch(UnitTest &ut)
{
    vector<double> x(2);
    x[0] = 0;
    x[1] = 0;

    vector<double> xold = x;
    vector<double> fvec;
    double fold = fnorm(x, fvec, &func);

    vector<double> g(2), p(2);
    g[0] = 7.2*fvec[0] - 2.2*fvec[1];
    g[1] = 3.5*fvec[0] + 2.7*fvec[1];

    p[0] = -0.4675753;
    p[1] =  1.6190125;

    double f;
    bool check;
    lnsrch(xold,
	   fold,
	   g, 
	   p,
	   x,
	   f,
	   check, 
	   fvec,
	   &func,
	   0.01,
           0.0);

    if (!check)
    {
	ut.failure("lnsrch did NOT bomb gracefully");
    }
    else
    {
	ut.passes("lnsrch bombed gracefully");
    }

    for (unsigned i=0; i<p.size(); ++i) p[i] = -p[i];
    lnsrch(xold,
	   fold,
	   g, 
	   p,
	   x,
	   f,
	   check, 
	   fvec,
	   &func,
	   0.01,
           0.0);

    if (check || f>1.0e-10)
    {
	ut.failure("lnsrch did NOT succeed");
    }
    else
    {
	ut.passes("lnsrch successful");
    }
   
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        UnitTest ut( argc,argv,release);
	tstlnsrch(ut);
    }
    catch (exception &err)
    {
	cout << "ERROR: While testing tstlnsrch, " << err.what() << endl;
	return 1;
    }
    catch( ... )
    {
	cout << "ERROR: While testing tstlnsrch, " 
             << "An unknown exception was thrown." << endl;
	return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstlnsrch.cc
//---------------------------------------------------------------------------//
