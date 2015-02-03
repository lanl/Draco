//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/scatterv.t.hh
 * \author Kent G. Budge
 * \date   Thu Mar 21 16:56:17 2002
 * \brief  C4 MPI template implementation.
 * \note   Copyright (C) 2002-2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef c4_scatterv_t_hh
#define c4_scatterv_t_hh

#include "c4/config.h"
#include "scatterv.hh"
#include "C4_Functions.hh"
#include "ds++/Assert.hh"
#include <algorithm>

namespace rtt_c4
{
using std::vector;
using std::copy;

//---------------------------------------------------------------------------//
// SCATTER
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
template<class T>
DLL_PUBLIC void indeterminate_scatterv(vector<vector<T> >  &outgoing_data,
                                       vector<T>           &incoming_data)
{
#ifdef C4_MPI
    { // This block is a no-op for with-c4=scalar 

        unsigned const N = rtt_c4::nodes();

        if (rtt_c4::node()==0)
        {
            vector<int> counts(N), displs(N);
            unsigned total_count = 0;
            for (unsigned p=0; p<N; ++p)
            {
                unsigned const n = outgoing_data[p].size();
                counts[p] = n;
                displs[p] = total_count;
                total_count += n;
            }
            int count;
            Remember(int check =)
            scatter(&counts[0], &count, 1);
            Check(check==MPI_SUCCESS);
            incoming_data.resize(count);
            
            vector<T> sendbuf(total_count);
            for (unsigned p=0; p<N; ++p)
            {
                Check(outgoing_data[p].size()+displs[p] <= sendbuf.size());
                copy(outgoing_data[p].begin(),
                     outgoing_data[p].end(),
                     sendbuf.begin()+displs[p]);
            }
            Remember(check =)
                rtt_c4::scatterv( (sendbuf.size()>0? &sendbuf[0] : NULL),
                                  (counts.size()>0? &counts[0] : NULL),
                                  (displs.size()>0? &displs[0] : NULL),
                                  (incoming_data.size()>0 ? &incoming_data[0] : NULL),
                                  count);
            Check(check == MPI_SUCCESS);
        }
        else
        {
            int count;
            Remember(int check =)
            scatter(static_cast<int*>(NULL), &count, 1);
            Check(check==MPI_SUCCESS);
            incoming_data.resize(count);
            Remember(check =)
                rtt_c4::scatterv(static_cast<T*>(NULL),
                                 static_cast<int*>(NULL),
                                 static_cast<int*>(NULL),
                                 (incoming_data.size()>0? &incoming_data[0] : NULL),
                                 count);
            Check(check == MPI_SUCCESS);
        }
    }
#else
    {
        // Only need to copy outgoing to incoming
        incoming_data = outgoing_data[0];
    }
#endif // C4_MPI

    return;
}

//---------------------------------------------------------------------------//
template<class T>
DLL_PUBLIC void determinate_scatterv(vector<vector<T> >  &outgoing_data,
                                     vector<T>           &incoming_data)
{
    Require(static_cast<int>(outgoing_data.size())==rtt_c4::nodes());
    
#ifdef C4_MPI
    { // This block is a no-op for with-c4=scalar 

        unsigned const N = rtt_c4::nodes();

        if (rtt_c4::node()==0)
        {
            vector<int> counts(N), displs(N);
            unsigned total_count = 0;
            for (unsigned p=0; p<N; ++p)
            {
                unsigned const n = outgoing_data[p].size();
                counts[p] = n;
                displs[p] = total_count;
                total_count += n;
            }
            int count = counts[0];
            Check(static_cast<int>(incoming_data.size())==count);
            
            vector<T> sendbuf(total_count);
            for (unsigned p=0; p<N; ++p)
            {
                copy(outgoing_data[p].begin(),
                     outgoing_data[p].end(),
                     sendbuf.begin()+displs[p]);
            }
            rtt_c4::scatterv(
                (sendbuf.size()>0 ? &sendbuf[0]: NULL),
                &counts[0],
                &displs[0],
                (incoming_data.size()>0?&incoming_data[0]:NULL),
                count);
        }
        else
        {
            int count = incoming_data.size();
            rtt_c4::scatterv(static_cast<T*>(NULL),
                             static_cast<int*>(NULL),
                             static_cast<int*>(NULL),
                             &incoming_data[0],
                             count);
        }
    }
#else
    {
        // Only need to copy outgoing to incoming
        incoming_data = outgoing_data[0];
    }
#endif // C4_MPI

    return;
}

} // end namespace rtt_c4

#endif // c4_scatterv_t_hh

//---------------------------------------------------------------------------//
// end of c4/scatterv.t.hh
//---------------------------------------------------------------------------//
