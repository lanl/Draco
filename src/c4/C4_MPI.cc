//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_MPI.cc
 * \author Thomas M. Evans
 * \date   Thu Mar 21 16:56:17 2002
 * \brief  C4 MPI implementation.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <c4/config.h>
#ifdef C4_MPI

#include <unistd.h>
#include <sys/times.h>
#include <vector>
#include "C4_Functions.hh"
#include "C4_Req.hh"

#include "C4_MPI.hh"

namespace rtt_c4
{

//---------------------------------------------------------------------------//
// MPI COMMUNICATOR
//---------------------------------------------------------------------------//

MPI_Comm communicator = MPI_COMM_WORLD;
bool initialized(false);

//---------------------------------------------------------------------------//
// Null source/destination rank
//---------------------------------------------------------------------------//

const int proc_null = MPI_PROC_NULL;

//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//

void initialize(int &argc, char **&argv)
{
    int result = MPI_Init(&argc, &argv);
    initialized = (result == MPI_SUCCESS);
    Check( initialized );

    // Resync clocks for Darwin mpich
    double foo( MPI_Wtick() );
}

//---------------------------------------------------------------------------//

void finalize()
{
    MPI_Finalize();
}

//---------------------------------------------------------------------------//

void free_inherited_comm()
{
    if (communicator != MPI_COMM_WORLD)
    {
	MPI_Comm_free(&communicator);
	communicator = MPI_COMM_WORLD;
    }
}

//---------------------------------------------------------------------------//
// QUERY FUNCTIONS
//---------------------------------------------------------------------------//

int node()
{
    int node = 0; 
    MPI_Comm_rank(communicator, &node);
    Check (node >= 0);
    return node;
}

//---------------------------------------------------------------------------//

int nodes()
{
    int nodes = 0;
    MPI_Comm_size(communicator, &nodes);
    Check (nodes > 0);
    return nodes;
}

//---------------------------------------------------------------------------//
// BARRIER FUNCTIONS
//---------------------------------------------------------------------------//

void global_barrier()
{
    MPI_Barrier(communicator);
}

//---------------------------------------------------------------------------//
// TIMING FUNCTIONS
//---------------------------------------------------------------------------//

// overloaded function (no args)
double wall_clock_time()
{
    return MPI_Wtime();
}
// overloaded function (provide POSIX timer information).
double wall_clock_time( tms & now )
{
    // obtain posix timer information and return it to the user via the
    // reference value argument "now".
    times( &now );
    // This funtion will return the MPI wall-clock time.
    return MPI_Wtime();
}

//---------------------------------------------------------------------------//

double wall_clock_resolution()
{
    return MPI_Wtick();
}

//---------------------------------------------------------------------------//
// PROBE/WAIT FUNCTIONS
//---------------------------------------------------------------------------//

bool probe(int  source, 
	   int  tag,
	   int &message_size)
{
    Require(source>=0 && source<nodes());

    int flag;
    MPI_Status status;
    
    // post an MPI_Irecv (non-blocking receive)
    MPI_Iprobe(source, tag, communicator, &flag, &status);

    if (!flag) return false;

    MPI_Get_count(&status, MPI_CHAR, &message_size);
    
    return true;
}

void blocking_probe(int  source, 
                    int  tag,
                    int &message_size)
{
    Require(source>=0 && source<nodes());

    MPI_Status status;
    MPI_Probe(source, tag, communicator, &status);
    MPI_Get_count(&status, MPI_CHAR, &message_size);
}

void wait_all(int count,
              C4_Req *requests)
{
    using std::vector;
    
    vector<MPI_Request> array_of_requests(count);
    vector<MPI_Status> array_of_statuses(count);
    for (unsigned i=0; i<count; ++i)
    {
        array_of_requests[i] = requests[i].r();
    }
    MPI_Waitall(count, &array_of_requests[0], &array_of_statuses[0]);
}

//---------------------------------------------------------------------------//
// ABORT
//---------------------------------------------------------------------------//

int abort(int error)
{
    // This test is not recorded as tested by BullseyeCoverage because abort
    // terminates the execution and BullseyeCoverage only reports coverage for
    // function that return control to main().

    int rerror = MPI_Abort(communicator, error);
    return rerror;
}

//---------------------------------------------------------------------------//

bool isScalar()
{
    return ! initialized;
}

} // end namespace rtt_c4

#endif // C4_MPI

//---------------------------------------------------------------------------//
//                              end of C4_MPI.cc
//---------------------------------------------------------------------------//
