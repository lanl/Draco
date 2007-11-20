//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstInvert_Comm_Map.cc
 * \author Mike Buksas
 * \date   Mon Nov 19 16:33:08 2007
 * \brief  
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "../Release.hh"
#include "../ParallelUnitTest.hh"

#include "../Invert_Comm_Map.hh"

using namespace std;
using namespace rtt_c4;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void test2(rtt_c4::ParallelUnitTest& ut)
{

    const int node  = rtt_c4::node();
    const int nodes = rtt_c4::nodes();

    std::vector<int> to_nodes;

    if (node == 0)
    {
        to_nodes.resize(1);
        to_nodes[0] = 1;
    }

    if (node == 1)
    {
        // No communication from node 1.
        to_nodes.resize(0);
    }


    std::vector<int> from_nodes(0);

    invert_comm_map(to_nodes, from_nodes);

    if (node == 0)
    {
        if (from_nodes.size() != 0) ut.failure("Incorrect map size on node 0.");
    }

    if (node == 1)
    {
        if (from_nodes.size() != 1) ut.failure("Incorrect size of map on node 1.");
        if (from_nodes[0] != 0)     ut.failure("Incorrect map contents on node 1.");
    }

}


void test4(rtt_c4::ParallelUnitTest& ut)
{

    const int node = rtt_c4::node();
    const int nodes = rtt_c4::nodes();

    std::list<int> to_nodes;

    if (node == 0)
    {
        to_nodes.push_back(1);
        to_nodes.push_back(2);
        to_nodes.push_back(3);
    }
    if (node == 1)
    {
        to_nodes.push_back(0);
    }
    if (node == 2)
    {
        to_nodes.push_back(0);
    }
    if (node == 3)
    {
        to_nodes.push_back(0);
    }

    std::list<int> from_nodes(0);

    invert_comm_map(to_nodes, from_nodes);

    if (node == 0)
    {
        if (from_nodes.size() != 3) ut.failure("Incorrect map size on node 0");
        for (int i = 1; i<=3; ++i)
        {
            if (from_nodes.front() != i)
                ut.failure("Incorrent map contents on node 0");
            from_nodes.pop_front();
        }
    }
    else
    {
        if (from_nodes.size() != 1)  ut.failure("Incorrect map size.");
        if (from_nodes.front() != 0) ut.failure("Incorrect map contents.");
    }
            

}


void test_n_to_n(rtt_c4::ParallelUnitTest& ut)
{

    const int node  = rtt_c4::node();
    const int nodes = rtt_c4::nodes();

    std::list<int> to_nodes;
    for (int i=0; i < nodes; ++i) to_nodes.push_back(i);

    std::list<int> from_nodes;
    invert_comm_map(to_nodes, from_nodes);

    for (int i=0; i < nodes; ++i)
    {
        if (to_nodes.front() != from_nodes.front())
            ut.failure("Incorrect data in map.");

        to_nodes.pop_front();
        from_nodes.pop_front();
    }

}

void test_cyclic(rtt_c4::ParallelUnitTest& ut)
{

    const int node = rtt_c4::node();
    const int nodes = rtt_c4::nodes();

    std::vector<int> to_nodes(1);
    to_nodes[0] = (node + 1) % nodes;

    std::vector<int> from_nodes;
    invert_comm_map(to_nodes, from_nodes);

    if (from_nodes.size() != 1) ut.failure("Incorrect map size.");
    if (from_nodes[0] != (node+nodes-1) % nodes)
        ut.failure("Incorrect map contents in cyclc test.");

}


void test_empty(rtt_c4::ParallelUnitTest& ut)
{

    std::vector<int> to_nodes;
    std::vector<int> from_nodes;

    invert_comm_map(to_nodes, from_nodes);

    if (from_nodes.size() != 0) ut.failure("Incorrect map size in empty test.");

}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut(argc, argv, release);
    try
    {
        if (nodes() == 2)
            test2(ut);

        if (nodes() == 4)
            test4(ut);

        test_n_to_n(ut);

        test_cyclic(ut);

        test_empty(ut);

        ut.passes("Just Because.");
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstInvert_Comm_Map, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstInvert_Comm_Map, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstInvert_Comm_Map.cc
//---------------------------------------------------------------------------//
