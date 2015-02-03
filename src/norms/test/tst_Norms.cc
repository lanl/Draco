//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   norms/test/tst_Norms.cc
 * \author Rob Lowrie
 * \date   Fri Jan 14 09:12:16 2005
 * \brief  Tests Norms.
 * \note   Copyright (C) 2004-2015 Los Alamos National Security, LLC. 
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../Norms.hh"
#include "../Norms_Labeled.hh"
#include "../Norms_Proc.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/Release.hh"
#include <sstream>

#define UNIT_TEST(A) if(!(A)) ITFAILS

using namespace std;

using rtt_norms::Norms;
using rtt_norms::Norms_Labeled;
using rtt_norms::Norms_Proc;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

// Test Norms
void test_Norms(rtt_dsxx::UnitTest & ut)
{
    Norms norms;

    // All norms should be the same if we add just a single quantity, in this
    // case, the processor number + 1.

    size_t x   = rtt_c4::node();
    double xp1 = x + 1.0;

    norms.add(xp1, x);

    UNIT_TEST(soft_equiv(norms.L1(), xp1));
    UNIT_TEST(soft_equiv(norms.L2(), xp1));
    UNIT_TEST(soft_equiv(norms.Linf(), xp1));
    UNIT_TEST(norms.index_Linf() == x);

    // Check accumulating the results back to the host proc.

    norms.comm(0);

    if ( x == 0 )
    {
        size_t num_nodes = rtt_c4::nodes();

        UNIT_TEST(norms.index_Linf() == num_nodes-1);
        UNIT_TEST(soft_equiv(norms.Linf(), double(num_nodes)));

        double tL1 = 0.0;
        double tL2 = 0.0;

        for ( size_t i = 0; i < num_nodes; ++i )
        {
            double xt = i + 1.0;
            tL1 += xt;
            tL2 += xt * xt;
        }

        tL1 /= num_nodes;
        tL2 = std::sqrt(tL2 / num_nodes);

        UNIT_TEST(soft_equiv(norms.L1(), tL1));
        UNIT_TEST(soft_equiv(norms.L2(), tL2));
    }

    { // test assignment

        Norms nc;
        nc = norms;
        UNIT_TEST(nc == norms);
    }

    { // test copy ctor.

        Norms nc(norms);
        UNIT_TEST(nc == norms);
    }

    // Done testing

    if ( ut.numFails==0 )
        PASSMSG("test_Norms() ok.");    
    return;
}

//---------------------------------------------------------------------------//
// Test Norms_Labeled
void test_Norms_Labeled(rtt_dsxx::UnitTest & ut)
{
    Norms_Labeled norms;

    // All norms should be the same if we add just a single quantity, in this
    // case, the processor number + 1.

    size_t x   = rtt_c4::node();
    double xp1 = x + 1.0;

    Norms_Labeled::Index indx(x);
    ostringstream o;
    o << "proc " << x;
    indx.label = o.str();

    norms.add(xp1, indx);

    UNIT_TEST(soft_equiv(norms.L1(), xp1));
    UNIT_TEST(soft_equiv(norms.L2(), xp1));
    UNIT_TEST(soft_equiv(norms.Linf(), xp1));
    UNIT_TEST(norms.index_Linf() == indx);

    // Check accumulating the results back to the host proc.

    norms.comm(0);

    if ( x == 0 )
    {
        size_t num_nodes = rtt_c4::nodes();

        //UNIT_TEST(norms.index_Linf() == num_nodes-1);
        UNIT_TEST(soft_equiv(norms.Linf(), double(num_nodes)));

        double tL1 = 0.0;
        double tL2 = 0.0;

        for ( size_t i = 0; i < num_nodes; ++i )
        {
            double xt = i + 1.0;
            tL1 += xt;
            tL2 += xt * xt;
        }

        tL1 /= num_nodes;
        tL2 = std::sqrt(tL2 / num_nodes);

        UNIT_TEST(soft_equiv(norms.L1(), tL1));
        UNIT_TEST(soft_equiv(norms.L2(), tL2));
    }

    { // test assignment

        Norms_Labeled nc;
        nc = norms;
        UNIT_TEST(nc == norms);
    }

    { // test copy ctor.

        Norms_Labeled nc(norms);
        UNIT_TEST(nc == norms);
    }

    // Done testing

    if ( ut.numFails == 0 )
        PASSMSG("test_Norms_Labeled() ok.");
    return;
}

//---------------------------------------------------------------------------//
// Test Norms_Proc
void test_Norms_Proc(rtt_dsxx::UnitTest & ut)
{
    Norms_Proc norms;

    // All norms should be the same if we add just a single quantity, in this
    // case, the processor number + 1.

    size_t x   = rtt_c4::node();
    double xp1 = x + 1.0;

    Norms_Proc::Index indx(x);

    norms.add(xp1, indx);

    UNIT_TEST(soft_equiv(norms.L1(), xp1));
    UNIT_TEST(soft_equiv(norms.L2(), xp1));
    UNIT_TEST(soft_equiv(norms.Linf(), xp1));
    UNIT_TEST(norms.index_Linf() == indx);

    // Check accumulating the results back to the host proc.

    norms.comm(0);

    if ( x == 0 )
    {
        size_t num_nodes = rtt_c4::nodes();

        //UNIT_TEST(norms.index_Linf() == num_nodes-1);
        UNIT_TEST(soft_equiv(norms.Linf(), double(num_nodes)));

        double tL1 = 0.0;
        double tL2 = 0.0;

        for ( size_t i = 0; i < num_nodes; ++i )
        {
            double xt = i + 1.0;
            tL1 += xt;
            tL2 += xt * xt;
        }

        tL1 /= num_nodes;
        tL2 = std::sqrt(tL2 / num_nodes);

        UNIT_TEST(soft_equiv(norms.L1(), tL1));
        UNIT_TEST(soft_equiv(norms.L2(), tL2));
    }

    { // test assignment

        Norms_Proc nc;
        nc = norms;
        UNIT_TEST(nc == norms);
    }

    { // test copy ctor.

        Norms_Proc nc(norms);
        UNIT_TEST(nc == norms);
    }

    // Done testing

    if ( ut.numFails==0 )
        PASSMSG("test_Norms_Proc() ok.");
    return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
    try
    {
        test_Norms(ut);
        test_Norms_Labeled(ut);
        test_Norms_Proc(ut);
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tst_Norms.cc
//---------------------------------------------------------------------------//
