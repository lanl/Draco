//----------------------------------*-C++-*----------------------------------//
// global_shmem.cc
// Geoffrey Furnish
// 24 May 1996
//---------------------------------------------------------------------------//
// @> Global C4 functions for Cray SHMEM.
//---------------------------------------------------------------------------//

#include "ds++/DynArray.hh"
#include "ds++/Assert.hh"

#include "c4/global_shmem.hh"

#include <mpp/shmem.h>

#define shmalloc malloc
#define shfree free

//---------------------------------------------------------------------------//
// Miscellaneous

#include <malloc.h>
#include <stdio.h>
#include <string.h>

C4_NAMESPACE_BEG

static void C4_shm_init_scalar_work_arrays();
static void C4_shm_init_pt2pt();

static int C4_shm_mynode;
static int C4_shm_nodes;

//---------------------------------------------------------------------------//
// Set up the work buffers we need in order to ensure that addresses we pass
// to shmem are in the address space.  Also, various state information to
// ensure that we have reliable, deterministic point to point communication
// semantics. 

static int *pe_msg_waiting, *pe_ready;

// const int C4_max_buf_sz = 1024;
const int C4_max_buf_sz = 4*1024;
const int C4_max_asyncs = 20;

struct pe_recv_buf_s {
    int length;
    int tag;
    char data[ C4_max_buf_sz ];
};

static pe_recv_buf_s *pe_recv_buf;

static int *pe_async_recv_reqs;
static Async_DB **pe_async_recv_req;

static int *pe_async_send_reqs;
static Async_DB **pe_async_send_req;

static int *pe_async_recvs_pending;
static int *pe_posted_recvs;
static Async_DB **pe_posted_recv;

static int *pe_async_sends_pending;
static int *pe_posted_sends;
static Async_DB **pe_posted_send;

static void C4_shm_init_pt2pt()
{
    int i, nodes = C4_shm_nodes;

    pe_msg_waiting = (int *) shmalloc( nodes * sizeof(int) );
    Assert( pe_msg_waiting );

    pe_ready = (int *) shmalloc( nodes * sizeof(int) );
    Assert( pe_ready );

    pe_recv_buf = (pe_recv_buf_s *)
	shmalloc( nodes * sizeof(struct pe_recv_buf_s) );
    Assert( pe_recv_buf );

// Setup the Async Recv structures.
    pe_async_recv_reqs = (int *) shmalloc( nodes * sizeof(int) );
    Assert( pe_async_recv_reqs );

    pe_async_recv_req = (Async_DB **) shmalloc( nodes * sizeof(Async_DB *) );
    for( i=0; i < nodes; i++ ) {
	pe_async_recv_reqs[i] = 0;

	pe_async_recv_req[i] = (Async_DB *) shmalloc( C4_max_asyncs *
						      sizeof(Async_DB) );
	Assert( pe_async_recv_req[i] );
	for( int j=0; j < C4_max_asyncs; j++ ) {
	// Initialize each Async_DB to a sane initial state.
	}
    }

// Setup the Async Send structures.

    pe_async_send_reqs = (int *) shmalloc( nodes * sizeof(int) );
    Assert( pe_async_send_reqs );

    pe_async_send_req = (Async_DB **) shmalloc( nodes * sizeof(Async_DB *) );
    for( i=0; i < nodes; i++ ) {
	pe_async_send_reqs[i] = 0;

	pe_async_send_req[i] = (Async_DB *) shmalloc( C4_max_asyncs *
						      sizeof(Async_DB) );
	Assert( pe_async_send_req[i] );
	for( int j=0; j < C4_max_asyncs; j++ ) {
	// Initialize each Async_DB to a sane initial state.
	}
    }

// Stuff for managing async postings.

    pe_async_recvs_pending = (int *) shmalloc( nodes * sizeof(int) );
    pe_posted_recvs        = (int *) shmalloc( nodes * sizeof(int) );
    for( i=0; i < nodes; i++ ) {
	pe_async_recvs_pending[i] = 0;
	pe_posted_recvs[i] = 0;
    }

    pe_posted_recv = (Async_DB **) shmalloc( nodes * sizeof(Async_DB *) );
    for( i=0; i < nodes; i++ ) {
	pe_posted_recv[i] = (Async_DB *) shmalloc( C4_max_asyncs *
						   sizeof(Async_DB) );
	Assert( pe_posted_recv[i] );
	for( int j=0; j < C4_max_asyncs; j++ ) {
	// Initialize each Async_DB to a sane initial state.
	}
    }

    pe_async_sends_pending = (int *) shmalloc( nodes * sizeof(int) );
    pe_posted_sends        = (int *) shmalloc( nodes * sizeof(int) );
    for( i=0; i < nodes; i++ ) {
	pe_async_sends_pending[i] = 0;
	pe_posted_sends[i] = 0;
    }

    pe_posted_send = (Async_DB **) shmalloc( nodes * sizeof(Async_DB *) );
    for( i=0; i < nodes; i++ ) {
	pe_posted_send[i] = (Async_DB *) shmalloc( C4_max_asyncs *
						   sizeof(Async_DB) );
	Assert( pe_posted_send[i] );
	for( int j=0; j < C4_max_asyncs; j++ ) {
	// Initialize each Async_DB to a sane initial state.
	}
    }

// Initialize the flag arrays.

    for( i=0; i < nodes; i++ ) {
	pe_msg_waiting[i] = 0;
	pe_ready[i] = 1;
    }
    gsync();
}

static Msg_Queue *msg_queue;

//---------------------------------------------------------------------------//
// Utility functions, not part of the public interface.

//---------------------------------------------------------------------------//
// Extract data from the inbox and store into a buffer.
//---------------------------------------------------------------------------//

static void pull_msg_from_inbox( int source, void *buf, int msglen )
{
    if ( msglen > C4_max_buf_sz ) {
	throw( "receive of packetized inbounds not implemented." );
    } else {
	memcpy( buf, pe_recv_buf[source].data, msglen );
    }
}

//---------------------------------------------------------------------------//
// Update book keeping data to reflect that the inbox is clear and source is
// free to send another message.
//---------------------------------------------------------------------------//

static void mark_inbox_as_clear( int source )
{
#if 0
// Clear my own msg_waiting flag, and tell soruce I'm ready to receive again.

// Must do it in this order to avoid a race.

    pe_msg_waiting[source] = 0;
    int one=1;
    shmem_put( (long *) &pe_ready[C4_shm_mynode],
	       (long *) &one, 1, source );
#endif
}

//---------------------------------------------------------------------------//
// Pull a message out of the inbox and stuff it onto the message queue.
//---------------------------------------------------------------------------//

static void empty_inbox_to_msg_queue( int source )
{
    if ( !pe_msg_waiting[source] ) return;

    int msgtag = pe_recv_buf[source].tag;
    int msglen = pe_recv_buf[source].length;

// Build the message to be put on our local queue.

    Msg_DB *pm = msg_queue[source].new_msg( msglen );
    pm->tag = msgtag;

    pull_msg_from_inbox( source, pm->buf, msglen );

// Okay, msg is now in pm, so clear the inbox.

    mark_inbox_as_clear( source );

// Now queue the message.

    msg_queue[source].enqueue(pm);
}

//---------------------------------------------------------------------------//
// Empty all inboxes except the one for the specified node.  This is used b/c
// the inbox for the specified node is being handled directly. 
//---------------------------------------------------------------------------//

static void empty_other_inboxes( int s )
{
    for( int i=0; i < C4_shm_nodes; i++ ) {
	if (i != s) {
	// queue messages from node i to msg_queue.
	}
    }
}

//---------------------------------------------------------------------------//
// This function searches the msg_queue for a given source node, to see if it
// contains a matching message.  Return value is 1 if a match is found, 0 if
// not.  If dequeue == 1 (which is the default), any matching message is
// copied to the user buffer and dequeued.
//---------------------------------------------------------------------------//

static int search_msg_queue( int source, void *buf, int size, int tag,
			     int dequeue =1 )
{
    int match = 0;

    for( int i=0; i < msg_queue[source].num_msgs(); i++ ) {

	Msg_DB *pm = msg_queue[source].msg(i);

	if ( pm->matches(tag) ) {
#ifdef SHMDBG
	    printf( "%d found matching msg on msg_queue\n", C4_shm_mynode );
#endif
	    match = 1;

	    if (dequeue) {

		Insist( size >= pm->msg_size,
			"recv buffer not large enough for message." );

	    // Pull data from message into user's output buffer.

		memcpy( buf, pm->buf, pm->msg_size );

		msg_queue[source].discard(i);
	    }
	}
    }

    return match;
}

//---------------------------------------------------------------------------//
// This function is used to clear the books for a receive request which had
// been posted to the source node, but which has by now been satisfied.
//---------------------------------------------------------------------------//

static void mark_recv_req_inactive( int source, int mid )
{
#if 0
    Async_DB& adb = pe_async_recv_req[source][mid];

    adb.state = Inactive;

    shmem_put( (long *) &pe_posted_recv[C4_shm_mynode][mid],
	       (long *) &adb, sizeof(Async_DB)/sizeof(double),
	       source );

    pe_posted_recvs[source]--;
    shmem_put( (long *) &pe_async_recvs_pending[C4_shm_mynode],
	       (long *) &pe_posted_recvs[source], 1, source );

    Assert( pe_posted_recvs[source] >= 0 );
    Assert( pe_posted_recvs[source] <= pe_async_recv_reqs[source]-1 );
#endif
}

//---------------------------------------------------------------------------//

void Init( int& argc, char **& argv )
{
    C4_shm_mynode = _my_pe();
    C4_shm_nodes = _num_pes();

    C4_shm_init_scalar_work_arrays();
    C4_shm_init_pt2pt();

    msg_queue = new Msg_Queue[ C4_shm_nodes ];

// Why is this not the default?
    shmem_set_cache_inv();
}

void Finalize()
{
// Whaddya think?  Should I shmem_free the symmetric data buffer?  heh heh.

// Check that the shmalloc books haven't been hosed.

//    Assert( !shmalloc_check(0) );

    shfree( pe_msg_waiting );
    shfree( pe_ready );
    shfree( pe_recv_buf );

    delete[] msg_queue;

// Check that the shmalloc books haven't been hosed.  Again, just to be sure.

//    Assert( !shmalloc_check(0) );
}

int node()
{
    return C4_shm_mynode;
}

int nodes()
{
    return C4_shm_nodes;
}

int group()
{
    int group = 0;
    return group;
}

void gsync()
{
    shmem_barrier_all();
}

//---------------------------------------------------------------------------//
// MPI send/receive calls (basic set)
//
// Synchronous:
//	MPI_Send(void* buf, int count, MPI_Datatype datatype,
//		 int dest, int tag, MPI_Comm comm);
//	MPI_Recv(void* buf, int count, MPI_Datatype datatype,
//		 int source, int tag, MPI_Comm comm, MPI_Status *status);
//
// Asynchronous:
//	MPI_Isend(void* buf, int count, MPI_Datatype datatype,
//		  int dest, int tag, MPI_Comm comm, MPI_Request *request);
//	MPI_Irecv(void* buf, int count, MPI_Datatype datatype,
//		  int source, int tag, MPI_Comm comm, MPI_Request *request);
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Perform a normal (blocking) send.
//---------------------------------------------------------------------------//

int SHM_Send( void *buf, int size, int dest, int tag, int group )
{
#if 0
// Check to see if we can short circuit send this to an already pending async
// receive posted by the dest node.

#ifdef SHMDBG
printf( "%d, pe_async_recvs_pending[%d]=%d\n",
	C4_shm_mynode, dest, pe_async_recvs_pending[dest] );
#endif

    if (pe_async_recvs_pending[dest]) {
	for( int i=0; i < C4_max_asyncs; i++ ) {
	    Async_DB& adb = pe_posted_recv[dest][i];
	    if (adb.state != Pending) continue;
	    if (adb.tag == C4_Any_Tag || adb.tag == tag) {
	    // Found a matching posted receive.  Send the data.

#ifdef SHMDBG
		printf( "%d, C4_Send, found matching posted async recv.\n",
			C4_shm_mynode );
#endif
	    // We should really do more checking than this.  We need not only
	    // to make sure that adb.size >= size, but really more precisely
	    // that dest buffer has enough words to accept the number of
	    // words we're sending.  For that matter, we should probably
	    // check alignment too.  But, ..., punt on that for now.

		Insist( size <= adb.size,
			"Receive buffer not large enough to send." );

		int words = size / 8 + ( size % 8 ? 1 : 0 );

		shmem_put( (long *) adb.data, (long *) buf, words, dest );

	    // Update the book keeping info.

		adb.state = No_Longer_Pending;

		shmem_put( (long *) &pe_async_recv_req[C4_shm_mynode][i],
			   (long *) &adb, sizeof(Async_DB)/sizeof(double),
			   dest );

	    // Note, to avoid races, we really can't decrement the pending
	    // receives count, but rather must let dest do it whenever it
	    // can. 

		return C4_SUCCESS;
	    }
	}
    }

// Okay, we couldn't short circuit send it to dest, so we'll have to wait for
// it to be ready to receive "the normal way".

    while( !pe_ready[dest] )
    // Do not have to flush cache, b/c C4_Init sets it for messages to
    // automatically invalidate the cache.
	;

    if (size > C4_max_buf_sz) {
    // packetize
	throw( "Sorry, packetized send not implemented yet." );
    } else {
    // put data to remote buffer, including tag and length

	struct {
	    int length;
	    int tag;
	} hdr;

	hdr.length = size; hdr.tag = tag;

	shmem_put( (long *) &pe_recv_buf[C4_shm_mynode],
 		   (long *) &hdr, 2, dest );

	int words = size / 8 + ( size % 8 ? 1 : 0 );
	Assert( words <= C4_max_buf_sz/8 );

	shmem_put( (long*) &pe_recv_buf[C4_shm_mynode].data,
 		   (long *) buf, words, dest );

    // indicate dest is not ready

	pe_ready[dest] = 0;

    // set msg_waiting flag in dest.

	int one=1;
	shmem_put( (long *) &pe_msg_waiting[C4_shm_mynode],
		   (long *) &one, 1, dest );
    }    
#endif
    return C4_SUCCESS;
}

//---------------------------------------------------------------------------//
// Perform a normal (blocking) receive.
//---------------------------------------------------------------------------//

int SHM_Recv( void *buf, int size, int source, int tag, int group )
{
    int msgtag, msglen = 0;

    if (source != C4_Any_Source) {
	int bytes_recvd = 0;

    // This is the "easy" case.  Look for a message from source which matches
    // tag. 

	Insist( source >= 0 && source < C4_shm_nodes,
		"C4_Recv: source node out of range." );

#ifdef SHMDBG
printf( "%d trying to receive a msg.\n", C4_shm_mynode );
#endif
    // First, check to see if the message queue is empty.  If it's not, the
    // message we want could be in there.

	if ( !msg_queue[source].is_empty() ) {

	// search msg queue. If find matching tag, copy to user buffer, and
	// remember to empty the inbox before returning.

	    int match = 0;
	    for( int i=0; i < msg_queue[source].num_msgs(); i++ ) {

		Msg_DB *pm = msg_queue[source].msg(i);

		if ( pm->matches(tag) ) {
		// Dequeue this message, copy to user buffer.

		    if (pm->msg_size > size) {
			throw( "recv buffer not large enough for message." );
		    }

		    bytes_recvd = pm->msg_size;
		    memcpy( buf, pm->buf, bytes_recvd );

		    msg_queue[source].discard(i);

		    match = 1;
		    break;
		}
	    }

	    if (match) {
	    // Queue any messages which are waiting.

	    // Mmm.  Failing to queue here should not change the delivery
	    // semantics.  So we will return to this speed optimization at a
	    // later time.

	    // And return.
		return bytes_recvd;
	    }
	}

    // It wasn't already in the queue.  So, ... let's just start pulling
    // messages in from the inbox.  If we find what we are looking for,
    // return it.  If not, queue it.

#ifdef SHMDBG
printf( "%d, it wasn't in the queue, looking in inbox.\n", C4_shm_mynode );
#endif
	int match = 0;
	do {
	// This is all bogus!  As written, it won't pull in an async send
	// which arrives at this point in time.  This will have to be
	// rewritten to handle that case correctly.

	// If there's nothing waiting for me right now, I'll presume I have
	// time to queue any inbound messages from other nodes.

	    if (!pe_msg_waiting[source])
		empty_other_inboxes(source);

	// Now busy wait till one shows up for me.

	    while( !pe_msg_waiting[source] )
		;

#ifdef SHMDBG
printf( "%d now there is a message waiting.\n", C4_shm_mynode );
#endif
	// Performance optimization opportunity:
	// empty_other_inboxes could be moved inside the busy wait.

	// Okay, something is finally here, now dispense.

	    msgtag = pe_recv_buf[source].tag;
	    msglen = pe_recv_buf[source].length;

#ifdef SHMDBG
printf( "%d, tag=%d, msgtag=%d, msglen=%d\n", C4_shm_mynode, tag, msgtag,
	msglen );
#endif
	    if (tag == C4_Any_Tag || tag == msgtag) {
#ifdef SHMDBG
printf( "%d, match!\n", C4_shm_mynode );
#endif
	    // Receive and exit this spinlock.
		if ( msglen > C4_max_buf_sz ) {
		    throw( "receive of packetized inbounds not implemented." );
		}
		else {
		    memcpy( buf, pe_recv_buf[source].data, msglen );
		}

		mark_inbox_as_clear( source );

		match = 1;
	    }
	    else {
		empty_inbox_to_msg_queue( source );

		msglen = 0;	// clear this, just to be safe.
	    }

	} while( !match );
    }
    else {
    // search for a message from any processor which matches tag

    // receive the message
	throw( "receive from any not implemented yet." );
    }

    return msglen;
}

//---------------------------------------------------------------------------//
// Perform a non blocking send.
//---------------------------------------------------------------------------//

C4_Req SHM_SendAsync( void *buf, int size, int dest, int tag, int group )
{
#if 0
// Actually, I think I should just go ahead and allocate an mid right now,
// b/c no matter what happens below, we still have to return the C4_Req in a
// valid state, since the user will have to wait() on it anyway.

    int mid = -1;

    if (pe_async_send_reqs[dest] == 0)
	mid = 0;
    else {
	for( int i=0; i < C4_max_asyncs; i++ )
	    if (pe_async_send_req[dest][i].state == Inactive) {
		mid = i;
		break;
	    }
	Insist( mid > -1, "Too many async messages." );
    }

    C4_Req r;
    r.set();

// C4_req's mid must encode which async message for which destination pe.
    r.mid = dest << 32 + mid;
    r.type = Async_Send;

    Async_DB& adb = pe_async_send_req[dest][mid];

    adb.state = Setting_Up;

    pe_async_send_reqs[dest]++;

#ifdef SHMDBG
    printf( "C4_SendAsync: mid=%d, r.mid=%x\n", mid, r.mid );
#endif
// Check to see if we can short circuit send this to an already pending async
// receive posted by the dest node.

    if (pe_async_recvs_pending[dest]) {
	for( int i=0; i < C4_max_asyncs; i++ ) {
	    Async_DB& pdb = pe_posted_recv[dest][i];
	    if (pdb.state != Pending) continue;
	    if (pdb.tag == C4_Any_Tag || pdb.tag == tag) {
	    // Found a matching posted receive.  Send the data.

#ifdef SHMDBG
		printf( "%d, C4_SendAsync, found matching posted async recv.\n",
			C4_shm_mynode );
#endif
	    // Should really do a better check on message buffer sizing
	    // stuff.  Later...

		Insist( size <= pdb.size,
			"Receive buffer not large enough to send." );

		int words = size / 8 + ( size % 8 ? 1 : 0 );

		shmem_put( (long *) pdb.data, (long *) buf, words, dest );

	    // Update the book keeping info.

		pdb.state = No_Longer_Pending;

		shmem_put( (long *) &pe_async_recv_req[C4_shm_mynode][i],
			   (long *) &pdb, sizeof(Async_DB)/sizeof(double),
			   dest );

	    // Note, to avoid races, we really can't decrement the pending
	    // receives count, but rather must let dest do it whenever it
	    // can. 

		adb.state = Complete;
		return r;
	    }
	}
    }

// Check if can transfer straight to recv buffer on remote node.  If so, then
// do so, and return.

    if (pe_ready[dest] && size <= C4_max_buf_sz) {
    // put data to remote buffer, including tag and length

	struct {
	    int length;
	    int tag;
	} hdr;

	hdr.length = size; hdr.tag = tag;

	shmem_put( (long *) &pe_recv_buf[C4_shm_mynode],
		   (long *) &hdr, 2, dest );

	int words = size/8;
	if (size % 8) 
	    words++;

#ifdef SHMDBG
	printf( "%d sending %d words to %d\n",
		C4_shm_mynode, words, dest );
	{
	    for( int i=0; i < words; i++ )
		printf( "word %d = %x\n", i, ((int *)buf)[i] );
	}
#endif

	shmem_put( (long*) &pe_recv_buf[C4_shm_mynode].data,
		   (long *) buf, words, dest );

    // indicate dest is not ready

	pe_ready[dest] = 0;

	adb.state = Complete;

    // set msg_waiting flag in dest.

	int one=1;
	shmem_put( (long *) &pe_msg_waiting[C4_shm_mynode],
		   (long *) &one, 1, dest );

#ifdef SHMDBG
	printf( "%d sent message to %d\n", C4_shm_mynode, dest );
#endif
	return r;
    }

// Post async send request with necessary data to remote node.

    throw( "Not able to post async send requests yet." );
#endif
    C4_Req r;
    return r;
}

//---------------------------------------------------------------------------//
// Perform a non blocking receive.
//---------------------------------------------------------------------------//

C4_Req SHM_RecvAsync( void *buf, int size, int source, int tag, int group )
{
#if 0
    Assert( source == C4_Any_Source || 
	    (source >= 0 && source < C4_shm_nodes) );
    Assert( size >= 0 );

#ifdef SHMDBG
printf( "%d entering C4_RecvAsync\n", C4_shm_mynode );
#endif

// Hmmm.  Contrarty to the above assertion, I don't yet have a decent plan
// for handling receive from any.  Punt on that for now...

    Insist( source != C4_Any_Source,
	    "Not able to handle receive-from-any yet." );

// Allocate an mid for this request.

    int mid = -1;

#ifdef SHMDBG
    printf( "%d searching for free Async_DB, pe_a_r_r[%d]=%d\n",
	    C4_shm_mynode, source, pe_async_recv_reqs[source] );
#endif

    if (pe_async_recv_reqs[source] == 0)
	mid = 0;
    else {
	for( int i=0; i < C4_max_asyncs; i++ ) {
#ifdef SHMDBG
	    printf( "%d checking mid %d, state=%d\n",
		    C4_shm_mynode, i, pe_async_recv_req[source][i].state );
#endif
	    if (pe_async_recv_req[source][i].state == Inactive) {
		mid = i;
		break;
	    }
	}
	Insist( mid > -1, "Too many async messages." );
    }

#ifdef SHMDBG
printf( "%d, mid=%d\n", C4_shm_mynode, mid );
#endif

    C4_Req r;
    r.set();

// C4_req's mid must encode which async message for which destination pe.
    r.mid = source << 32 + mid;
    r.type = Async_Recv;

    Async_DB& adb = pe_async_recv_req[source][mid];

    adb.state = Setting_Up;

    pe_async_recv_reqs[source]++;

// Check the msg queue to see if there is a match waiting for us already.

    if ( !msg_queue[source].is_empty() ) {

	int match = search_msg_queue( source, buf, size, tag );

	if (match) {
	    adb.state = Complete;

	// Think about proactively emptying the inbox.  This will induce a
	// memcpy which will slow down this node, but will free up the inbox
	// which might let a different node proceed faster.  

	// Have to think about this.  Note that it would be a bummer if the
	// msg in the inbox was a packetized one, since I don't think we'd
	// want to stall completion of this async request while doing a hairy
	// packet handshake.  Hmmmm.  May have to come back to this.

	    empty_inbox_to_msg_queue( source );

	    return r;
	}
    }

// Okay, there was no matching message in the msg_queue.  Now let's check the
// inbox.

    if (pe_msg_waiting[source]) {
#ifdef SHMDBG
printf( "%d msg waiting in inbox.\n", C4_shm_mynode );
#endif
    // If it matches, do the short circuit thing, otherwise queue it.

	int msgtag = pe_recv_buf[source].tag;
	int msglen = pe_recv_buf[source].length;

	if (tag == C4_Any_Tag || tag == msgtag) {
#ifdef SHMDBG
printf( "%d inbox was a match\n", C4_shm_mynode );
#endif
	// We have a match!

	    Insist( msglen <= size,
		    "Inbound message too large for data buffer." );

	    pull_msg_from_inbox( source, buf, msglen );

	    mark_inbox_as_clear( source );

	// Now update state of r, and go home.

	    adb.state = Complete;
	    return r;
	}
	else {
	// Queue this message.
	    empty_inbox_to_msg_queue( source );
	}
    }

// Check to see if any async sends have been posted to this node.  If so,
// check them for matches, and if a match, fetch data and return.

    if (pe_async_sends_pending[source]) {

	throw( "Don't know how to check for posted async sends to this node." );
    }

// Okay, nothings was here, and nothing was coming.  So, post an async
// receive request, and go about our business.

    adb.data = buf;
    adb.tag = tag;
    adb.size = size;
    adb.state = Pending;

//  shmem_put( target, source, words, dest );

    shmem_put( (long *) &pe_posted_recv[C4_shm_mynode][mid],
	       (long *) &adb, sizeof(Async_DB)/sizeof(double),
	       source );

    pe_posted_recvs[source]++;
    shmem_put( (long *) &pe_async_recvs_pending[C4_shm_mynode],
	       (long *) &pe_posted_recvs[source], 1, source );

#ifdef SHMDBG
printf( "%d posted async recv to pe %d\n", C4_shm_mynode, source );
printf( "%d pe_posted_recvs[%d]=%d\n", C4_shm_mynode, source,
	pe_posted_recvs[source] );
printf( "%d state=%d\n", C4_shm_mynode, adb.state );
#endif
#endif
C4_Req r; 
    return r;
}

//---------------------------------------------------------------------------//
// "Optimized" forms of the above.  These are supposed to avoid spurios
// object creation and copy.  However, for the purpose of getting going,
// these are implemented by just calling the ones above.  Later we can switch
// the implementations so that these go faster.
//---------------------------------------------------------------------------//

void SHM_SendAsync( C4_Req& r, void *buf, int size, int dest, int tag, 
		    int group /*=0*/ )
{
// Not checking that r is not in use, which is of course a concern...
    r = SHM_SendAsync( buf, size, dest, tag, group );
}

void SHM_RecvAsync( C4_Req& r, void *buf, int size, int source, int tag, 
		    int group /*=0*/ )
{
// Not checking that r is not in use, which is of course a concern...
    r = SHM_RecvAsync( buf, size, source, tag, group );
}

//---------------------------------------------------------------------------//
// Wait on an asynchronous message to complete.
//---------------------------------------------------------------------------//

void C4_Wait( int m, int t )
{
#ifdef SHMDBG
    printf( "%d waiting, C4_Req.mid = %x\n", C4_shm_mynode, m );
#endif
    int node = m >> 32 ;
    int mid = m & 0xFFFFFFFF;

#ifdef SHMDBG
    printf( "Trying to wait on mid %d for node %d, type=%s\n",
	    mid, node, 
	    (t == Async_Send ? "send" : "recv" ) );
#endif
    switch( t ) {

    case Async_Send:
	while( pe_async_send_req[node][mid].state != Complete )
	    ;
	pe_async_send_req[node][mid].state = Inactive;
	pe_async_send_reqs[node]--;
	Assert( pe_async_send_reqs[node] >= 0 );
	break;

    case Async_Recv: {
	int source = node;
	Async_DB& adb = pe_async_recv_req[node][mid];
	while( adb.state != Complete && adb.state != Inactive ) {
	// Process inbound messages.  sheesh, what a nightmare.

	// Okay, let's think about how horrid this could be.  The sender
	// could've sent the msg via short circuit, resulting in a
	// No_Longer_Pending state.  Or, it could've just missed the posting
	// of the recv req, and wound up sending the msg to the inbox (in
	// which case it could either be in the inbox, /or/ on the msg queue
	// by now).  Finally (I think), it could have been an async send
	// which missed the posting of the recv, so that it wound up here in
	// the posted inbound list (posting sends isn't implemented yet).

	// Do gobs of work to exhaust these possibilities.

	    if (adb.state == No_Longer_Pending) {
#ifdef SHMDBG
		printf( "%d, async recv msg mid=%d, No_Longer_Pending\n",
			C4_shm_mynode, mid );
#endif
		mark_recv_req_inactive( source, mid );
		continue;
	    }

	// Uhh, now come the hard parts...
#ifdef SHMDBG
printf( "%d, async recv msg mid=%d, state=%d\n",
	C4_shm_mynode, mid, adb.state );
#endif
	    if (!msg_queue[source].is_empty()) {
	    // Check the msg queue.
 		int match = search_msg_queue( source, adb.data,
					      adb.size, adb.tag );

		if (match) {
		    mark_recv_req_inactive( source, mid );
		    continue;
		}
	    }

	// Check the inbox.
	    if (pe_msg_waiting[source]) {
		int msgtag = pe_recv_buf[source].tag;
		int msglen = pe_recv_buf[source].length;

		if (adb.tag == C4_Any_Tag || adb.tag == msgtag) {
#ifdef SHMDBG
		    printf( "%d inbox was a match\n", C4_shm_mynode );
#endif
		// We have a match!

		    Insist( msglen <= adb.size,
			    "Inbound message too large for data buffer." );

		    pull_msg_from_inbox( source, adb.data, msglen );

		    mark_inbox_as_clear( source );

		// The only way to get here is with bizarre timing
		// situations, and I don't frankly know how to write a test
		// case which will exercise this code block.  Grrrrrrrrrr.

		// Ture hackers know how to use logic analyzers and how to
		// induce wait states.

		    mark_recv_req_inactive( source, mid );
		    continue;
		}
		else {
		// Queue this message.
		    empty_inbox_to_msg_queue( source );
		}
	    }
	    {
	    // Check to see if a posted send satisfies this request.  If so,
	    // pull it in, and clear all the books.  Remember that this recv
	    // req may have been posted to source node.  grrr.
		if (pe_async_sends_pending[source]) {

		    throw( "Don't know how to check for posted async sends to this node." );
		}
	    }

	// Okay, well, looks like we're gonna busy wait, so might as well
	// take this chance to do whatever else can be done to expedite
	// message delivery from other nodes.

	// NOTE: That is a speed optimization, not necessary for maintaining
	// proper semantics.  So, defer it for now.

	}
	adb.state = Inactive;
	pe_async_recv_reqs[node]--;
#ifdef SHMDBG
	printf( "%d, after waiting, pe_a_r_r[%d]=%d\n", C4_shm_mynode,
		node, pe_async_recv_reqs[node] );
	printf( "%d, after waiting, pe_posted_recvs[%d]=%d\n",
		C4_shm_mynode, node, pe_posted_recvs[node] );
#endif
	Assert( pe_async_recv_reqs[node] >= 0 );
	break;
    }

    default:
	throw( "Unknown message type, something is /really/ hosed." );
    }
}
#if 0
//---------------------------------------------------------------------------//
// Send a buffer of integers.
//---------------------------------------------------------------------------//

int C4_Send( int *buf, int nels, int dest, int group /*=0*/ )
{
    return C4_Send( (void *)buf, nels * sizeof(int),
		    dest, C4_int_ptr_Tag, group );
}

//---------------------------------------------------------------------------//
// Receive a buffer of integers.  nels is how many we can receive, the return
// value is how many we did receive.
//---------------------------------------------------------------------------//

int C4_Recv( int *buf, int nels, int source, int group /*=0*/ )
{
    int cnt = C4_Recv( (void *) buf, nels * sizeof(int),
		       source, C4_int_ptr_Tag, group );
    return cnt/sizeof(int);
}

//---------------------------------------------------------------------------//
// Send a buffer of floats.
//---------------------------------------------------------------------------//

int C4_Send( float *buf, int nels, int dest, int group /*=0*/ )
{
    return C4_Send( (void *)buf, nels * sizeof(float),
		    dest, C4_float_ptr_Tag, group );
}

//---------------------------------------------------------------------------//
// Receive a buffer of floats.  nels is how many we can receive, the return
// value is how many we did receive.
//---------------------------------------------------------------------------//

int C4_Recv( float *buf, int nels, int source, int group /*=0*/ )
{
    int cnt = C4_Recv( (void *) buf, nels * sizeof(float),
		       source, C4_float_ptr_Tag, group );
    return cnt/sizeof(float);
}

//---------------------------------------------------------------------------//
// Send a buffer of doubles.
//---------------------------------------------------------------------------//

int C4_Send( double *buf, int nels, int dest, int group /*=0*/ )
{
    return C4_Send( (void *)buf, nels * sizeof(double),
		    dest, C4_double_ptr_Tag, group );
}

//---------------------------------------------------------------------------//
// Receive a buffer of floats.  nels is how many we can receive, the return
// value is how many we did receive.
//---------------------------------------------------------------------------//

int C4_Recv( double *buf, int nels, int source, int group /*=0*/ )
{
    int cnt = C4_Recv( (void *) buf, nels * sizeof(double),
		       source, C4_double_ptr_Tag, group );
    return cnt/sizeof(double);
}
#endif
//---------------------------------------------------------------------------//
// Global reduction operations.
//
// The call most like that in NX is:
//
//	MPI_Allreduce(void* sendbuf, void* recvbuf, int count,
//		      MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
// 
// which returns the result to all processes in the group.  "op" determines
// the type of reduction performed.
//
// Available reduction operators, description, and allowed types: 
//
// MPI_MAX, MPI_MIN		min, max			I, F
// MPI_SUM, MPI_PROD		sum, product			I, F
// MPI_BAND, MPI_BOR, MPI_BXOR	bitwise and, or, xor		I, B
// MPI_LAND, MPI_LOR, MPI_LXOR	logical and, or, xor		I
// MPI_MAXLOC, MPI_MINLOC	min, max value and location
//
// where types are:
//
// I:	MPI_INT, MPI_LONG, MPI_SHORT, MPI_UNSIGNED_SHORT, 
//	MPI_UNSIGNED, MPI_UNSIGNED_LONG
// F:	MPI_FLOAT, MPI_DOUBLE, MPI_LONG_DOUBLE
// B:	MPI_BYTE
//---------------------------------------------------------------------------//

// static DynArray<int>    ibuf(10);
// static DynArray<long>   lbuf(10);
// static DynArray<float>  fbuf(10);
// static DynArray<double> dbuf(10);

// "sr" == "scalar reduce"

static int sri_pWrk[ _SHMEM_REDUCE_MIN_WRKDATA_SIZE ];
static float srf_pWrk[ _SHMEM_REDUCE_MIN_WRKDATA_SIZE ];
static double srd_pWrk[ _SHMEM_REDUCE_MIN_WRKDATA_SIZE ];
static long sr_pSync[ _SHMEM_COLLECT_SYNC_SIZE ];

static void C4_shm_init_scalar_work_arrays()
{
    for( int i=0; i < _SHMEM_COLLECT_SYNC_SIZE; i++ ) {
	sr_pSync[i] = _SHMEM_SYNC_VALUE;
    }

    gsync();
}

//---------------------------------------------------------------------------//
// Sum, scalar

static int sr_i;

void gsum( int& x )
{
    sr_i = x;

    shmem_int_sum_to_all( &sr_i, &sr_i, 1, 0, 0, C4_shm_nodes,
			  sri_pWrk, sr_pSync );
    x = sr_i;

    gsync();
}

void gsum( long& x )
{
    sr_i = (int) x;		// Should be valid on an Alpha.

    shmem_int_sum_to_all( &sr_i, &sr_i, 1, 0, 0, C4_shm_nodes,
			  sri_pWrk, sr_pSync );

    x = sr_i;

    gsync();
}

static float sr_f;

void gsum( float& x )
{
    sr_f = x;

    shmem_float_sum_to_all( &sr_f, &sr_f, 1, 0, 0, C4_shm_nodes,
			    srf_pWrk, sr_pSync );

    x = sr_f;

    gsync();
}

static double sr_d;

void gsum( double& x )
{
    sr_d = x;

    shmem_double_sum_to_all( &sr_d, &sr_d, 1, 0, 0, C4_shm_nodes,
			     srd_pWrk, sr_pSync );

    x = sr_d;

    gsync();
}

#ifdef max
#undef max
#endif
#define max(a,b) ( (a) > (b) ? (a) : (b) )

//---------------------------------------------------------------------------//
// Sum, array

static int *ar_i = NULL;
static int ari_sz = 0;
static int *ari_pWrk = NULL;
static int ari_Wsz = 0;

void gsum( int *px, int n )
{
    int i;
    if (n < 1) return;

// Make sure all the work arrays are big enough.
    if (ari_sz < n) {
    // Make the source/target array larger.
	if (ar_i)
	    shfree( ar_i );
	ar_i = (int *) shmalloc( n * sizeof(int) );
	ari_sz = n;

    // See if we have to lengthen pWrk.
	int nWsz = max( n/2 + 1, _SHMEM_REDUCE_MIN_WRKDATA_SIZE );

	if (nWsz > ari_Wsz) {
	    if (ari_pWrk)
		shfree( ari_pWrk );

	    ari_pWrk = (int *) shmalloc( nWsz * sizeof(int) );
	    ari_Wsz = nWsz;
	}
    }

// Copy data to symmetric work array.
    for( i=0; i < n; i++ )
	ar_i[i] = px[i];

// Perform reduction.
    shmem_int_sum_to_all( ar_i, ar_i, n, 0, 0, C4_shm_nodes,
			  ari_pWrk, sr_pSync );

// Copy data back to user buffer.
    for( i=0; i < n; i++ )
	px[i] = ar_i[i];

    gsync();
}

void gsum( long *px, int n )
{
    gsum( (int *) px, n );
}

static float *ar_f = NULL;
static int arf_sz = 0;
static float *arf_pWrk = NULL;
static int arf_Wsz = 0;

void gsum( float *px, int n )
{
    int i;
    if (n < 1) return;

// Make sure all the work arrays are big enough.
    if (arf_sz < n) {
    // Make the source/target array larger.
	if (ar_f)
	    shfree( ar_f );
	ar_f = (float *) shmalloc( n * sizeof(float) );
	arf_sz = n;

    // See if we have to lengthen pWrk.
	int nWsz = max( n/2 + 1, _SHMEM_REDUCE_MIN_WRKDATA_SIZE );

	if (nWsz > arf_Wsz) {
	    if (arf_pWrk)
		shfree( arf_pWrk );

	    arf_pWrk = (float *) shmalloc( nWsz * sizeof(float) );
	    arf_Wsz = nWsz;
	}
    }

// Copy data to symmetric work array.
    for( i=0; i < n; i++ )
	ar_f[i] = px[i];

// Perform reduction.
    shmem_float_sum_to_all( ar_f, ar_f, n, 0, 0, C4_shm_nodes,
			    arf_pWrk, sr_pSync );

// Copy data back to user buffer.
    for( i=0; i < n; i++ )
	px[i] = ar_f[i];

    gsync();
}

static double *ar_d = NULL;
static int ard_sz = 0;
static double *ard_pWrk = NULL;
static int ard_Wsz = 0;

void gsum( double *px, int n )
{
    int i;
    if (n < 1) return;

// Make sure all the work arrays are big enough.
    if (ard_sz < n) {
    // Make the source/target array larger.
	if (ar_d)
	    shfree( ar_d );
	ar_d = (double *) shmalloc( n * sizeof(double) );
	ard_sz = n;

    // See if we have to lengthen pWrk.
	int nWsz = max( n/2 + 1, _SHMEM_REDUCE_MIN_WRKDATA_SIZE );

	if (nWsz > ard_Wsz) {
	    if (ard_pWrk)
		shfree( ard_pWrk );

	    ard_pWrk = (double *) shmalloc( nWsz * sizeof(double) );
	    ard_Wsz = nWsz;
	}
    }

// Copy data to symmetric work array.
    for( i=0; i < n; i++ )
	ar_d[i] = px[i];

// Perform reduction.
    shmem_double_sum_to_all( ar_d, ar_d, n, 0, 0, C4_shm_nodes,
			     ard_pWrk, sr_pSync );

// Copy data back to user buffer.
    for( i=0; i < n; i++ )
	px[i] = ar_d[i];

    gsync();
}

//---------------------------------------------------------------------------//
// Min, scalar

void gmin( int& x )
{
    sr_i = x;

    shmem_int_min_to_all( &sr_i, &sr_i, 1, 0, 0, C4_shm_nodes,
			  sri_pWrk, sr_pSync );
    x = sr_i;

    gsync();
}

void gmin( long& x )
{
    sr_i = (int) x;		// Should be valid on an Alpha.

    shmem_int_min_to_all( &sr_i, &sr_i, 1, 0, 0, C4_shm_nodes,
			  sri_pWrk, sr_pSync );

    x = sr_i;

    gsync();
}

void gmin( float& x )
{
    sr_f = x;

    shmem_float_min_to_all( &sr_f, &sr_f, 1, 0, 0, C4_shm_nodes,
			    srf_pWrk, sr_pSync );

    x = sr_f;

    gsync();
}

void gmin( double& x )
{
    sr_d = x;

    shmem_double_min_to_all( &sr_d, &sr_d, 1, 0, 0, C4_shm_nodes,
			     srd_pWrk, sr_pSync );

    x = sr_d;

    gsync();
}

//---------------------------------------------------------------------------//
// Min, array.

void gmin( int *px, int n )
{
    int i;
    if (n < 1) return;

// Make sure all the work arrays are big enough.
    if (ari_sz < n) {
    // Make the source/target array larger.
	if (ar_i)
	    shfree( ar_i );
	ar_i = (int *) shmalloc( n * sizeof(int) );
	ari_sz = n;

    // See if we have to lengthen pWrk.
	int nWsz = max( n/2 + 1, _SHMEM_REDUCE_MIN_WRKDATA_SIZE );

	if (nWsz > ari_Wsz) {
	    if (ari_pWrk)
		shfree( ari_pWrk );

	    ari_pWrk = (int *) shmalloc( nWsz * sizeof(int) );
	    ari_Wsz = nWsz;
	}
    }

// Copy data to symmetric work array.
    for( i=0; i < n; i++ )
	ar_i[i] = px[i];

// Perform reduction.
    shmem_int_min_to_all( ar_i, ar_i, n, 0, 0, C4_shm_nodes,
			  ari_pWrk, sr_pSync );

// Copy data back to user buffer.
    for( i=0; i < n; i++ )
	px[i] = ar_i[i];

    gsync();
}

void gmin( long *px, int n )
{
    gmin( (int *) px, n );
}

void gmin( float *px, int n )
{
    int i;
    if (n < 1) return;

// Make sure all the work arrays are big enough.
    if (arf_sz < n) {
    // Make the source/target array larger.
	if (ar_f)
	    shfree( ar_f );
	ar_f = (float *) shmalloc( n * sizeof(float) );
	arf_sz = n;

    // See if we have to lengthen pWrk.
	int nWsz = max( n/2 + 1, _SHMEM_REDUCE_MIN_WRKDATA_SIZE );

	if (nWsz > arf_Wsz) {
	    if (arf_pWrk)
		shfree( arf_pWrk );

	    arf_pWrk = (float *) shmalloc( nWsz * sizeof(float) );
	    arf_Wsz = nWsz;
	}
    }

// Copy data to symmetric work array.
    for( i=0; i < n; i++ )
	ar_f[i] = px[i];

// Perform reduction.
    shmem_float_min_to_all( ar_f, ar_f, n, 0, 0, C4_shm_nodes,
			    arf_pWrk, sr_pSync );

// Copy data back to user buffer.
    for( i=0; i < n; i++ )
	px[i] = ar_f[i];

    gsync();
}

void gmin( double *px, int n )
{
    int i;
    if (n < 1) return;

// Make sure all the work arrays are big enough.
    if (ard_sz < n) {
    // Make the source/target array larger.
	if (ar_d)
	    shfree( ar_d );
	ar_d = (double *) shmalloc( n * sizeof(double) );
	ard_sz = n;

    // See if we have to lengthen pWrk.
	int nWsz = max( n/2 + 1, _SHMEM_REDUCE_MIN_WRKDATA_SIZE );

	if (nWsz > ard_Wsz) {
	    if (ard_pWrk)
		shfree( ard_pWrk );

	    ard_pWrk = (double *) shmalloc( nWsz * sizeof(double) );
	    ard_Wsz = nWsz;
	}
    }

// Copy data to symmetric work array.
    for( i=0; i < n; i++ )
	ar_d[i] = px[i];

// Perform reduction.
    shmem_double_min_to_all( ar_d, ar_d, n, 0, 0, C4_shm_nodes,
			     ard_pWrk, sr_pSync );

// Copy data back to user buffer.
    for( i=0; i < n; i++ )
	px[i] = ar_d[i];

    gsync();
}

//---------------------------------------------------------------------------//
// Max, scalar

void gmax( int& x )
{
    sr_i = x;

    shmem_int_max_to_all( &sr_i, &sr_i, 1, 0, 0, C4_shm_nodes,
			  sri_pWrk, sr_pSync );
    x = sr_i;

    gsync();
}

void gmax( long& x )
{
    sr_i = (int) x;		// Should be valid on an Alpha.

    shmem_int_max_to_all( &sr_i, &sr_i, 1, 0, 0, C4_shm_nodes,
			  sri_pWrk, sr_pSync );

    x = sr_i;

    gsync();
}

void gmax( float& x )
{
    sr_f = x;

    shmem_float_max_to_all( &sr_f, &sr_f, 1, 0, 0, C4_shm_nodes,
			    srf_pWrk, sr_pSync );

    x = sr_f;

    gsync();
}

void gmax( double& x )
{
    sr_d = x;

    shmem_double_max_to_all( &sr_d, &sr_d, 1, 0, 0, C4_shm_nodes,
			     srd_pWrk, sr_pSync );

    x = sr_d;

    gsync();
}

//---------------------------------------------------------------------------//
// Max, array.

void gmax( int *px, int n )
{
    int i;
    if (n < 1) return;

// Make sure all the work arrays are big enough.
    if (ari_sz < n) {
    // Make the source/target array larger.
	if (ar_i)
	    shfree( ar_i );
	ar_i = (int *) shmalloc( n * sizeof(int) );
	ari_sz = n;

    // See if we have to lengthen pWrk.
	int nWsz = max( n/2 + 1, _SHMEM_REDUCE_MIN_WRKDATA_SIZE );

	if (nWsz > ari_Wsz) {
	    if (ari_pWrk)
		shfree( ari_pWrk );

	    ari_pWrk = (int *) shmalloc( nWsz * sizeof(int) );
	    ari_Wsz = nWsz;
	}
    }

// Copy data to symmetric work array.
    for( i=0; i < n; i++ )
	ar_i[i] = px[i];

// Perform reduction.
    shmem_int_max_to_all( ar_i, ar_i, n, 0, 0, C4_shm_nodes,
			  ari_pWrk, sr_pSync );

// Copy data back to user buffer.
    for( i=0; i < n; i++ )
	px[i] = ar_i[i];

    gsync();
}

void gmax( long *px, int n )
{
    gmax( (int *) px, n );
}


void gmax( float *px, int n )
{
    int i;
    if (n < 1) return;

// Make sure all the work arrays are big enough.
    if (arf_sz < n) {
    // Make the source/target array larger.
	if (ar_f)
	    shfree( ar_f );
	ar_f = (float *) shmalloc( n * sizeof(float) );
	arf_sz = n;

    // See if we have to lengthen pWrk.
	int nWsz = max( n/2 + 1, _SHMEM_REDUCE_MIN_WRKDATA_SIZE );

	if (nWsz > arf_Wsz) {
	    if (arf_pWrk)
		shfree( arf_pWrk );

	    arf_pWrk = (float *) shmalloc( nWsz * sizeof(float) );
	    arf_Wsz = nWsz;
	}
    }

// Copy data to symmetric work array.
    for( i=0; i < n; i++ )
	ar_f[i] = px[i];

// Perform reduction.
    shmem_float_max_to_all( ar_f, ar_f, n, 0, 0, C4_shm_nodes,
			    arf_pWrk, sr_pSync );

// Copy data back to user buffer.
    for( i=0; i < n; i++ )
	px[i] = ar_f[i];

    gsync();
}

void gmax( double *px, int n )
{
    int i;
    if (n < 1) return;

// Make sure all the work arrays are big enough.
    if (ard_sz < n) {
    // Make the source/target array larger.
	if (ar_d)
	    shfree( ar_d );
	ar_d = (double *) shmalloc( n * sizeof(double) );
	ard_sz = n;

    // See if we have to lengthen pWrk.
	int nWsz = max( n/2 + 1, _SHMEM_REDUCE_MIN_WRKDATA_SIZE );

	if (nWsz > ard_Wsz) {
	    if (ard_pWrk)
		shfree( ard_pWrk );

	    ard_pWrk = (double *) shmalloc( nWsz * sizeof(double) );
	    ard_Wsz = nWsz;
	}
    }

// Copy data to symmetric work array.
    for( i=0; i < n; i++ )
	ar_d[i] = px[i];

// Perform reduction.
    shmem_double_max_to_all( ar_d, ar_d, n, 0, 0, C4_shm_nodes,
			     ard_pWrk, sr_pSync );

// Copy data back to user buffer.
    for( i=0; i < n; i++ )
	px[i] = ar_d[i];

    gsync();
}

//---------------------------------------------------------------------------//
// Now some debugging functions, just in case we ever have any bugs :-).

void C4_shm_dbg_1()
{
#if 0
    static int seq;
    seq = 0;
    gsync();

    if (C4_shm_mynode == 0) seq = 1;

    while(!seq) ;

    cout << "C4_shm_dbg_1:\n" << flush;
    {
	for( int i=0; i < C4_shm_nodes; i++ ) {
	    printf( "node %d: pe_ready[%d]=%d\n",
		    C4_shm_mynode, i, pe_ready[i] );
	    printf( "node %d: pe_msg_waiting[%d]=%d\n",
		    C4_shm_mynode, i, pe_msg_waiting[i] );
	    printf( "node %d: pe_async_recv_reqs[%d]=%d\n",
		    C4_shm_mynode, i, pe_async_recv_reqs[i] );
	    printf( "node %d: msg_queue[%d].is_empty()=%d\n",
		    C4_shm_mynode, i, msg_queue[i].is_empty() );
	}
    }

    if (C4_shm_mynode < C4_shm_nodes-1) {
	shmem_put( (long *) &seq, (long *) &seq, 1, C4_shm_mynode+1 );
    }

    gsync();
#endif
}

C4_NAMESPACE_END

#include "ds++/DynArray.cc"

template class DynArray<C4::Msg_DB *>;

//---------------------------------------------------------------------------//
//                              end of global_shmem.cc
//---------------------------------------------------------------------------//
