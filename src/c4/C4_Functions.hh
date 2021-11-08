//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/C4_Functions.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 11:42:03 2002
 * \brief  C4 Communication Functions.
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved.
 *
 * This file contains the declarations for communication functions provided by C4. This file allows
 * the client to include the message passing services provided by C4.  The function declarations and
 * class definitions are contained in the rtt_c4 namespace.  For backwards compatibility, the
 * old-style C4 functions and classes are declared in the C4 namespace.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_c4_C4_Functions_hh
#define rtt_c4_C4_Functions_hh

#include "C4_Datatype.hh"
#include "C4_Status.hh"
#include "C4_Traits.hh"
#include "C4_sys_times.h"

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
/*! Forward declarations
 *
 * We postpone including C4_Req.hh until C4_MPI.i.hh is loaded. This allows the 'friend'
 * declarations found in class C4_Req to be seen after the 'official' function declarations (with
 * optional default arguments) are loaded.
 */
class C4_Req;

//------------------------------------------------------------------------------------------------//
/*!
 * C4 unit tests.
 */
/*! \example c4/test/tstAbort.cc
 * Example of MPI abort functions.
 */
/*! \example c4/test/tstBroadcast.cc
 * Example of MPI broadcast-like functions
 */
/*! \example c4/test/tstComm_Dup.cc
 * Example
 */
/*! \example c4/test/tstPingPong.cc
 * Example of point-to-point communications
 */
/*! \example c4/test/tstReduction.cc
 * Example of many-to-one communications
 */
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
// GLOBAL CONSTANTS
//------------------------------------------------------------------------------------------------//

//! Any source rank
DLL_PUBLIC_c4 extern const int any_source;
//! Null source/destination rank
DLL_PUBLIC_c4 extern const int proc_null;

//------------------------------------------------------------------------------------------------//
// SETUP FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Initialize a parallel job.
 */
int initialize(int &argc, char **&argv, int required = DRACO_MPI_THREAD_SINGLE);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Finish a parallel job.
 */
void finalize();

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Inherit a communicator from another application.
 */
template <typename Comm> void inherit(const Comm &);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Free an inherited communicator from another application.
 */
void free_inherited_comm();

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Create up a new vector type.
 *
 * \param count Number of blocks in the data type
 * \param blocklength Length of each block (in units of base type)
 * \param stride Spacing between start of each block (in units of base type)
 * \param new_type On return, contains the new type descriptor.
 */
template <typename T>
int create_vector_type(unsigned count, unsigned blocklength, unsigned stride,
                       C4_Datatype &new_type);

//------------------------------------------------------------------------------------------------//
//! Free a user defined type, such as a vector type.

void type_free(C4_Datatype &old_type);

//------------------------------------------------------------------------------------------------//
// QUERY FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Get the node (rank) of the current processor.
 *
 * The rank is determined by the current communicator.
 */
int node();
uint32_t rank();

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Get the number of processors used for this job.
 *
 * The number of nodes is determined by the current communicator.
 */
int nodes();
uint32_t nranks();

//------------------------------------------------------------------------------------------------//
// BARRIER FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Set a global barrier for the communicator.
 */
void global_barrier();

//------------------------------------------------------------------------------------------------//
// BLOCKING SEND/RECEIVE OPERATIONS
//------------------------------------------------------------------------------------------------//
//! Do a point-to-point, blocking send.
template <typename T>
int send(const T *buffer, int size, int destination, int tag = C4_Traits<T *>::tag);

//------------------------------------------------------------------------------------------------//
//! Do a point-to-point, blocking send.
template <typename T> int send_custom(const T *buffer, int size, int destination, int tag);

//------------------------------------------------------------------------------------------------//
//! Do a point-to-point, blocking receive.
template <typename T> int receive(T *buffer, int size, int source, int tag = C4_Traits<T *>::tag);

//------------------------------------------------------------------------------------------------//
//! Do a point-to-point, blocking receive with a custom MPI type
template <typename T> int receive_custom(T *buffer, int size, int source, int tag);

//------------------------------------------------------------------------------------------------//
//! Do a point-to-point, blocking send of a user-defined type.
template <typename T>
int send_udt(const T *buffer, int size, int destination, C4_Datatype &,
             int tag = C4_Traits<T *>::tag);

//------------------------------------------------------------------------------------------------//
//! Do a point-to-point, blocking send-receive.
template <typename TS, typename TR>
int send_receive(TS *sendbuf, int sendcount, int destination, TR *recvbuf, int recvcount,
                 int source, int sendtag = C4_Traits<TS *>::tag,
                 int recvtag = C4_Traits<TR *>::tag);

//------------------------------------------------------------------------------------------------//
//! Do a point-to-point, blocking receive of a user-defined type.
template <typename T>
int receive_udt(T *buffer, int size, int source, C4_Datatype &, int tag = C4_Traits<T *>::tag);

//------------------------------------------------------------------------------------------------//
// NON-BLOCKING SEND/RECEIVE OPERATIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking send.
 *
 * \return C4_Req object to handle communication requests
 */
template <typename T>
C4_Req send_async(T const *buffer, int size, int destination, int tag = C4_Traits<T *>::tag);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking send.
 */
template <typename T>
void send_async(C4_Req &request, T const *buffer, int size, int destination,
                int tag = C4_Traits<T *>::tag);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking synchronous send.
 */
template <typename T>
void send_is(C4_Req &request, T const *buffer, int size, int destination, int tag);
template <typename T> void send_is(C4_Req &request, T const *buffer, int size, int destination) {
  int tag = C4_Traits<T *>::tag;
  send_is(request, buffer, size, destination, tag);
  return;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Get the size of a message with custom types
 *
 * \param[in] status C4_Status object that will hold MPI request status
 * \param[in] mpi_type The signature of the special type.
 * \return number of type T objects in the completed message
 */
template <typename T> int message_size_custom(C4_Status status, const T &mpi_type);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking send with a MPI custom type
 *
 * \param[in,out] request C4_Req object that will hold MPI request
 * \param[in,out] buffer array of data of type T that has an MPI type
 * \param[in] size buffer size
 * \param[in] destination rank that will receive this message
 * \param[in] tag message tag
 */
template <typename T>
void send_is_custom(C4_Req &request, const T *buffer, int size, int destination,
                    int tag = C4_Traits<T *>::tag);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking receive.
 *
 * \return C4_Req object to handle communication requests
 */
template <typename T>
C4_Req receive_async(T *buffer, int size, int source, int tag = C4_Traits<T *>::tag);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking receive.
 */
template <typename T>
void receive_async(C4_Req &request, T *buffer, int size, int source, int tag = C4_Traits<T *>::tag);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Post a non-blocking receive for a message of custom MPI type data
 *
 * \param[in,out] request C4_Req object that will hold MPI request
 * \param[in,out] buffer array of data of type T that has a registered MPI type
 * \param[in] size buffer size
 * \param[in] source remote rank sending message to this rank
 * \param[in] tag message tag
 */
template <typename T>
void receive_async_custom(C4_Req &request, T *buffer, int size, int source,
                          int tag = C4_Traits<T *>::tag);

//------------------------------------------------------------------------------------------------//
// BROADCAST
//------------------------------------------------------------------------------------------------//

/*------------------------------------------------------------------------------------------------*/
/*
 * \brief Send data from processor 0 to all other processors.
 *
 * These are declared and defined in C4_MPI.hh and in C4_Serial.hh.  KT is having trouble with
 * getting the DLL_PUBLIC_c4 to be correct, so delay declaration until the C4_MPI.hh or C4_Serial.hh
 * files are included.

template <typename T>
int broadcast(T *buffer, int size, int root);

template <typename ForwardIterator, typename OutputIterator>
void broadcast(ForwardIterator first, ForwardIterator last,
                             OutputIterator result);

template <typename ForwardIterator, typename OutputIterator>
void broadcast(ForwardIterator first, ForwardIterator last,
                             OutputIterator result, OutputIterator result_end);
*/

//------------------------------------------------------------------------------------------------//
// GATHER/SCATTER
//------------------------------------------------------------------------------------------------//

template <typename T> int gather(T *send_buffer, T *receive_buffer, int size);

template <typename T> int allgather(T *send_buffer, T *receive_buffer, int size);

template <typename T>
int gatherv(T *send_buffer, int send_size, T *receive_buffer, int *receive_sizes,
            int *receive_displs);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Gather from all ranks to a receiving array on each rank, using sizes and displacements
 *
 * \param[in] send_buffer array of data of type T that has a registered MPI type
 * \param[in] send_size size of send buffer
 * \param[in,out] receive_buffer array to gather at, according to sizes and displacements per rank
 * \param[in] receive_sizes anticipated data size per rank in the Communicator.
 * \param[in] receive_displs displacement into receive_buffer per rank in the Communicator.
 */
template <typename T>
int allgatherv(T *send_buffer, int send_size, T *receive_buffer, int *receive_sizes,
               int *receive_displs);

template <typename T> int scatter(T *send_buffer, T *receive_buffer, int size);

template <typename T>
int scatterv(T *send_buffer, int *send_sizes, int *send_displs, T *receive_buffer,
             int receive_size);

//------------------------------------------------------------------------------------------------//
// GLOBAL REDUCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do a global sum of a scalar variable.
 */
template <typename T> void global_sum(T &x);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do a non-blocking global sum of a scalar variable.
 *
 * \param[in,out] send_buffer scalar value on this processing element
 * \param[in,out] recv_buffer scalar value summed across all ranks
 * \param[in,out] request C4_Requst handle for testing completed message
 */
template <typename T> void global_isum(T &send_buffer, T &recv_buffer, C4_Req &request);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do a global product of a scalar variable.
 */
template <typename T> void global_prod(T &x);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do a global minimum of a scalar variable.
 */
template <typename T> void global_min(T &x);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do a global maximum of a scalar variable.
 */
template <typename T> void global_max(T &x);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do a global logical and (MPI_LAND) of a scalar variable
 */
template <typename T> void global_and(T &x);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global sum of an array.
 *
 * Only instantiate this function for types L that are integral. Details of this pattern are at
 * https://en.cppreference.com/w/cpp/types/enable_if.
 */
template <typename T, typename L,
          typename std::enable_if<std::is_integral<L>::value, bool>::type = true>
void global_sum(T *x, L n);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global product of an array.
 */
template <typename T> void global_prod(T *x, int n);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global minimum of an array.
 */
template <typename T> void global_min(T *x, int n);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global maximum of an array.
 */
template <typename T> void global_max(T *x, int n);

//------------------------------------------------------------------------------------------------//
// TIMING FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the wall-clock time in seconds.
 */
double wall_clock_time();
double wall_clock_time(DRACO_TIME_TYPE &now);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the resolution of wall_clock_time.
 */
double wall_clock_resolution();

//------------------------------------------------------------------------------------------------//
// PROBE/WAIT FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief See if a message is pending.
 *
 * \param source
 * Processor from which a message may be pending.
 * \param tag
 * Tag for pending message.
 * \param message_size
 * On return, size of the pending message in bytes.
 * \return \c true if a message from the specified processor with the specified tag is pending; \c
 * false otherwise.
 */
bool probe(int source, int tag, int &message_size);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Wait until a message (of unknown size) is pending.
 *
 * \param source
 * Processor from which a message of unknown size is expected.
 * \param tag
 * Tag for pending message.
 * \param message_size
 * On return, size of the pending message in bytes.
 */
void blocking_probe(int source, int tag, int &message_size);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Wait until every one of a set of posted sends/receives is complete.
 *
 * This version returns no status information.
 *
 * \param count
 * Size of the set of requests to wait on.
 * \param requests
 * Set of requests to wait on.
 */
void wait_all(const unsigned count, C4_Req *const requests);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Wait until every one of a set of posted sends/receives is complete.
 *
 * This version returns source rank information for RECEIVE REQUESTS ONLY WARNING: If you call this
 * function with send requests, the source rank ID might not be set in MPI_STATUS. REPEAT: THE
 * RETURN VALUE IS VALID FOR RECEIVE REQUESTS ONLY!!!!
 *
 * \pre All requests passed to this function must be receive requests
 *
 * \param count
 * Size of the set of requests to wait on.
 * \param requests
 * Set of requests to wait on.
 * \return A vector containing the source rank IDs for received messages.
 */
std::vector<int> wait_all_with_source(const unsigned count, C4_Req *const requests);

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Wait until one of a set of posted sends/receives is complete.
 *
 * \param count
 * Size of the set of requests to wait on.
 * \param requests
 * Set of requests to wait on.
 * \return The request that completed.
 */
unsigned wait_any(unsigned count, C4_Req *requests);

//------------------------------------------------------------------------------------------------//
// ABORT
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Abort across all processors.
 *
 * \param error suggested return error, defaults to 1
 */
int abort(int error = 1);

//------------------------------------------------------------------------------------------------//
// isScalar
//------------------------------------------------------------------------------------------------//
//! Is C4 executing in scalar-only mode?
bool isScalar();

//! Has MPI been initialized (always false for scalar-mode).
bool isMpiInit();

//------------------------------------------------------------------------------------------------//
// get_processor_name
//------------------------------------------------------------------------------------------------//
//! Return the processor name for each rank.
std::string get_processor_name();

//------------------------------------------------------------------------------------------------//
// prefix_sum
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the value of the prefix sum at this processor.
 *
 * \param node_value Current node's value of variable to be prefix summed
 * \return Sum of value over nodes up to and including this node.
 */
template <typename T> T prefix_sum(const T node_value);

/*!
 * \brief Return the value of the prefix sum at this processor.
 *
 * \param buffer Current node's starting buffer address to be prefix summed
 * \param n number of objects of type T in the buffer
 */
template <typename T> void prefix_sum(T *buffer, const int32_t n);

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// Include the appropriate header for an underlying message passing
// implementation.

#ifdef C4_SCALAR
#include "C4_Serial.hh"
#endif

#ifdef C4_MPI
#include "C4_MPI.hh"
#endif

#endif // rtt_c4_C4_Functions_hh

//------------------------------------------------------------------------------------------------//
// end of c4/C4_Functions.hh
//------------------------------------------------------------------------------------------------//
