//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/DD_Transporter.t.hh
 * \author Todd J. Urbatsch and Thomas M. Evans
 * \date   Wed Apr 19 15:37:14 2000
 * \brief  DD_Transporter template definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_DD_Transporter_t_hh
#define rtt_imc_DD_Transporter_t_hh

#include "DD_Transporter.hh"
#include "mc/Particle_Stack.hh"
#include "ds++/Assert.hh"
#include <iostream>
#include <iomanip>

namespace rtt_imc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief DD_Transporter constructor.
 *
 * Builds an IMC Transporter for full and general DD topologies.  The
 * constructor checks to make sure that the proper Topology exists for this
 * construction.  
 */
template<class MT, class FT, class PT>
DD_Transporter<MT,FT,PT>::DD_Transporter(SP_Topology top)
    : Transporter<MT,FT,PT>(),
    topology(top),
    num_run(0),
    nsrc_run(0),
    num_done(0),
    num_to_run(0),
    delta_t(0),
    cycle(0),
    print_f(0)
{
    Require (!mesh);
    Require (!opacity);
    Require (!mat_state);
    Require (!source);
    Require (!tally);
    Require (!communicator);
    Require (topology);

    // check the topology
    Insist (topology->get_parallel_scheme() == "DD" ||
	    topology->get_parallel_scheme() == "DD/replication", 
	    "Invalid topology assigned to Transporter!");
}

//---------------------------------------------------------------------------//
// IMC DD TRANSPORT (PUBLIC INTERFACE)
//---------------------------------------------------------------------------//
/*!
 * \brief Do IMC, DD-based replication transport.
 
 * This does F&C IMC transport for one timestep.  The new census created
 * during the timestep is returned.  All fundamental objects are unset at the
 * end of the transport step.

 * All particles placed into the census have \b global cell numbers.
 * Components in rtt_imc expect to operate on census particles with global
 * cell numbers.

 * \param dt timestep size
 * \param cycle_in current cycle
 * \param print_f_in particle notification frequency

 * \param num_to_run_in total number of particles to run across all
 * processors

 * \param verbose do verbose messaging

 * \return census for this timestep

 */
template<class MT, class FT, class PT>
typename DD_Transporter<MT,FT,PT>::SP_Census 
DD_Transporter<MT,FT,PT>::transport(double dt,
				    int    cycle_in, 
				    int    print_f_in, 
				    int    num_to_run_in, 
				    bool   verbose)
{
    using std::cerr;
    using std::cout;
    using std::endl;
    using std::setw;

    Require (ready()); 

    // data initialized through function arguments
    num_to_run = num_to_run_in;
    delta_t    = dt;
    print_f    = print_f_in;
    cycle      = cycle_in;

    // initialize `finished' status to false
    finished = 0;

    // initialize total and source particle counters
    num_run  = 0;
    nsrc_run = 0;
    num_done = 0;

    // start transport
    cerr << ">> Doing transport for cycle " << cycle
	 << " on proc " << C4::node() << " using " 
	 << topology->get_parallel_scheme() << "." << endl;

    // particle history diagnostic
    rtt_dsxx::SP<typename PT::Diagnostic> check;
    if (verbose)
        check = new typename PT::Diagnostic(cout, true); 

    // make a census and communication bank on this node
    SP_Census new_census_bank(new Census());
    Bank bank;

    // post arecvs to communicating processors
    communicator->post();

    // post arecvs from each node to master and from master to all nodes
    post_step_arecvs();

    // begin timing the transport on this processor
    double trans_begin = C4::Wtime();

    // transport particles
    while (!finished)
    {
	// transport a source particle and any incoming particles
	if (*source)
	    trans_src_async(check, bank, new_census_bank);
	
	// transport an incoming particle from another domain (processor)
	else if (bank.size())
	    trans_domain_async(check, bank, new_census_bank);

	// if send-buffers are not empty, flush them
	else if (communicator->get_send_size())
	    communicator->flush();

	// receive particles as second-to-last option
	else if (communicator->arecv_post(bank))
	{
	    int bsize = bank.size();
	    Check (bank.size() > 0);
	}

	// report num_done back to master; see if we are finished
	else
	    update();
    }

    // stop timing the transport on this processor
    double trans_end = C4::Wtime();

    // synchronize before wrapping this cycle up
    C4::gsync();
    
    // complete the asynchronous recvs (last comm w/out posting new recvs)
    complete_step_arecvs();
    
    // clean up receive buffers 
    Check (!communicator->get_send_size());
    Check (!bank.size());
    communicator->asend_end();
    communicator->arecv_end();
    Check (!communicator->arecv_status());
    Check (!communicator->asend_status());

    // finished with this timestep
    cerr << ">> Finished particle transport for cycle " 
	 << cycle << " on proc " << C4::node() << " in " 
	 << trans_end - trans_begin << " seconds." << endl;

    // unset the transport stuff
    unset();

    Ensure (!ready());
    Ensure (new_census_bank);
    
    return new_census_bank;
}

//---------------------------------------------------------------------------//
// PRIVATE IMC TRANSPORT FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Run a source particle and, possibly, incoming particles.
 */
template<class MT, class FT, class PT>
void DD_Transporter<MT,FT,PT>::trans_src_async(
    SP_PT_Diagnostic check,
    Bank            &bank,
    SP_Census        new_census_bank)
{
    using rtt_dsxx::SP;
    using std::cerr;
    using std::endl;
    using std::setw;

    // get a source particle
    SP<PT> particle = source->get_Source_Particle(delta_t); 
    Check (particle->status());

    // transport the particle; update counters
    particle->transport(*mesh, *opacity, *tally, check);
    num_run++;
    nsrc_run++;
    
    // particle is no longer active
    Check (!particle->status());

    // increment number of completed particles, if this one is done
    if (particle->get_descriptor() == PT::CENSUS)
	num_done++;
    else if (particle->get_descriptor() == PT::ESCAPE)
	num_done++;
    else if (particle->get_descriptor() == PT::KILLED)
	num_done++;
   
    // if particle goes to census, write to new census bank
    if (particle->get_descriptor() == PT::CENSUS)
    {
	// convert the particle cell index back to a global cell index
	int local_cell  = particle->get_cell();
	int global_cell = topology->global_cell(local_cell);
	particle->set_cell(global_cell);
	
	// push particle to census
	new_census_bank->push(particle);
    }
    
    // if particle crosses a domain boundary, communicate the particle
    if (particle->get_descriptor() == PT::CROSS_BOUNDARY)
        communicator->communicate(particle);

    // The conditional has been commented out because too many MPI buffers
    // were building up.  To check the incoming buffer for every source
    // particle is probably too inefficient but more robust.  Maybe it should
    // be checked every max(1,0.01*buffer_size/nprocs) or something more
    // optimal between 1 and buffer_size.  For now, we'll do it after every
    // particle. 
    // 
    // >>> clip <<< 
    // if we have transported the same number of particles as
    // >>> the buffer // size; then... if (!(nsrc_run %
    // >>> Particle_Buffer<PT>::get_maximum_num_particles())) 
    // >>> clip <<<
    //
    // check on, receive, and transport the incoming buffer
    if (communicator->arecv_post(bank))
	while (bank.size())
	    trans_domain_async(check, bank, new_census_bank);
    
    // message particle counter
    if (!(nsrc_run % print_f)) 
	cerr << "Ran " << setw(10) << nsrc_run << " source (and " << setw(10) 
	     << (num_run - nsrc_run) << " incoming) particles on proc " 
	     << C4::node() << endl;
}

//---------------------------------------------------------------------------//
/*!  
 * \brief Run an incoming particle and, possibly, check for more incoming
 * particles. 
 */
template<class MT, class FT, class PT>
void DD_Transporter<MT,FT,PT>::trans_domain_async(
    SP_PT_Diagnostic check, 
    Bank            &bank,
    SP_Census        new_census_bank)
{
    using rtt_dsxx::SP;
    using std::cerr;
    using std::endl;
    using std::setw;

    // get a particle from the bank and activate it
    SP<PT> particle = bank.top();
    bank.pop();
    particle->reset_status();
    
    // transport the particle
    particle->transport(*mesh, *opacity, *tally, check);
    num_run++;
    
    // particle is no longer active; take further action accordingly
    Check (!particle->status());

    // increment number of completed particles, if this one is done
    if (particle->get_descriptor() == PT::CENSUS)
	num_done++;
    else if (particle->get_descriptor() == PT::ESCAPE)
	num_done++;
    else if (particle->get_descriptor() == PT::KILLED)
	num_done++;

    if (particle->get_descriptor() == PT::CENSUS)
    {
	// convert the particle cell index back to a global cell index
	int local_cell  = particle->get_cell();
	int global_cell = topology->global_cell(local_cell);
	particle->set_cell(global_cell);
	
	// push particle to census
	new_census_bank->push(particle);
    }
    
    if (particle->get_descriptor() == PT::CROSS_BOUNDARY)
	communicator->communicate(particle);
  
    // during the source block, this statement could fill the bank such that
    // all incoming particles are run before another source particle is run.
    // As in trans_src_async, checking for incoming particles after running
    // num_run particles was building up too many buffers.  Thus, we have
    // commented out the conditional such that we check for incoming
    // particles after every incoming particle is run.  Checking so
    // frequently is more inefficient, but seems to be more robust.  Ideally, 
    // we would want to check at some optimal frequency, somewhere between 1
    // and buffer_size.  
    // if (!(num_run % Particle_Buffer<PT>::get_maximum_num_particles()))
    communicator->arecv_post(bank);

    // message particle counter
    if (!(num_run % print_f)) 
	cerr << "Total of " << setw(10) << num_run 
	     << " particles run on proc " << C4::node() << endl;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Communicate number of particles transported from IMC nodes to
 * master node and check if finished.
 */
template<class MT, class FT, class PT>
void DD_Transporter<MT,FT,PT>::update()
{
    using C4::C4_Req;

    if (C4::node())
    {
	// if IMC node has completed any particles, inform master
	if (num_done > 0)
	{
	    C4_Req sender = C4::SendAsync(&num_done, 1, 0, 500);
	    sender.wait();
	    num_done = 0;
	}

	// Check on "finished" status from master
	if (rcv501_fin.complete())
	{
	    Check (recv_finished == 1);
	    finished = recv_finished;
	    // won't post new receive; should be finished
	}
    }
    else if (!C4::node())
    {
	// make sure (probably too late) num_to_run is reasonable
        Require (num_to_run > 0);

	// see if the master has completed the last particle,
	if (num_done == num_to_run)
	    finished = 1;

	// check if any IMC nodes have finished more particles
	int i = 0;
	while (!finished && ++i < C4::nodes())
	    if (rcv500_ndone[i-1].complete())
	    {
		Check (recv_num_done[i-1] > 0);
		num_done += recv_num_done[i-1];
		Check (num_done <= num_to_run);
		rcv500_ndone[i-1] = C4::RecvAsync(&recv_num_done[i-1], 1, i,
						  500);  
		if (num_done == num_to_run)
		    finished = 1;
	    }
	
	// master sends 'finished' to all other nodes.
	if (finished)
	    for (int i = 1; i < C4::nodes(); i++)
	    {
		C4_Req sender = C4::SendAsync(&finished, 1, i, 501);
		sender.wait();
	    }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Post asynchronous receives.
 */
template<class MT, class FT, class PT>
void DD_Transporter<MT,FT,PT>::post_step_arecvs()
{
    using C4::nodes;
    using C4::RecvAsync;

    // master
    if (!C4::node())
    {
	// resize C4 request tags to number of other nodes
	rcv500_ndone.resize(nodes()-1);
	
	// resize flags from other nodes
	recv_num_done.resize(nodes()-1);

	// post first receives to all other nodes
	for (int i = 1; i < nodes(); i++)
	{
	    rcv500_ndone[i-1] = RecvAsync(&recv_num_done[i-1], 1, i, 500);
	}
    }

    // IMC nodes
    else if (C4::node())
    {
	// post first (and only) `receive' to the master
	rcv501_fin = RecvAsync(&recv_finished, 1, 0, 501);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Complete asynchronous receives.
 */
template<class MT, class FT, class PT>
void DD_Transporter<MT,FT,PT>::complete_step_arecvs()
{
    using C4::C4_Req;

    // IMC nodes
    if (C4::node())
    {
	// IMC nodes sending final "num_done" to master   
	num_done = 0;
	C4_Req sender = C4::SendAsync(&num_done, 1, 0, 500);
	sender.wait();
    }

    // master
    else if (!C4::node())
    {
        int r500count = 0;

	// receive final "num_done" from IMC nodes
	for (int i = 1; i < C4::nodes(); i++)
	    while (!rcv500_ndone[i-1].complete())
		r500count++;
    }
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Set fundamental objects for a transport step.
 
 * This function sets the fundamental transport objects: MT, Opacity,
 * Mat_State, Source, Tally, and Communicator for a current timestep.  After
 * the transport function is executed these objects are automatically unset
 * using the unset function.  Thus, the set is valid for \b one execution of
 * transport.

 * The Communicator must be assigned for DD_Transport topologies.
 
 * \param mesh_in rtt_dsxx::SP to a fully replicated mesh.
 * \param mat_state_in rtt_dsxx::SP to a valid Mat_State object
 * \param opacity_in rtt_dsxx::SP to a valid Opacity object
 * \param source_in rtt_dsxx::SP to a valid Source object
 * \param tally_in rtt_dsxx::SP to a valid Tally object
 * \param random_walk_in rtt_dsxx::SP to a random walk object (can be null)
 * \param communicator_in rtt_dsxx::SP to a valid Communicator object

 */
template<class MT, class FT, class PT>
void DD_Transporter<MT,FT,PT>::set(SP_Mesh         mesh_in,
				   SP_Mat_State    mat_state_in,
				   SP_Opacity      opacity_in,
				   SP_Source       source_in,
				   SP_Tally        tally_in,
				   SP_Random_Walk  random_walk_in,
				   SP_Communicator communicator_in)
{
    Require (mesh_in);
    Require (opacity_in);
    Require (source_in);
    Require (mat_state_in);
    Require (tally_in);
    Require (communicator_in);

    // assign objects
    mesh         = mesh_in;
    opacity      = opacity_in;
    source       = source_in;
    mat_state    = mat_state_in;
    tally        = tally_in;
    random_walk  = random_walk_in;
    communicator = communicator_in;
    
    // number of global cells is the same number of cells on processor
    int num_cells = topology->num_cells(C4::node());

    Ensure (num_cells == mesh->num_cells());
    Ensure (num_cells == opacity->num_cells());
    Ensure (num_cells == source->num_cells());
    Ensure (num_cells == mat_state->num_cells());
    Ensure (num_cells == tally->num_cells());
    Ensure (random_walk ? num_cells == random_walk->num_cells() : true);
    Ensure (communicator);
    Ensure (topology);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unset transport objects.

 * This function is used to unset the fundamental transport objects assigned
 * in set.  It is provided in the public interface as a memory conservation
 * feature.  It is automatically called at the end of transport.  We include
 * it in the public interface so that users may optionally unset object.

 * The function does \b not check to see if the fundamental objects are
 * assigned prior to unsetting them.

 */
template<class MT, class FT, class PT>
void DD_Transporter<MT,FT,PT>::unset()
{
    Require (topology);

    // assign the fundamental objects to null pointers
    mesh         = SP_Mesh();
    opacity      = SP_Opacity();
    mat_state    = SP_Mat_State();
    tally        = SP_Tally();
    source       = SP_Source();
    communicator = SP_Communicator();
    random_walk  = SP_Random_Walk();

    Ensure (!mesh);
    Ensure (!opacity);
    Ensure (!mat_state);
    Ensure (!tally);
    Ensure (!source);
    Ensure (!random_walk);
    Ensure (!communicator);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Query to see if all objects are ready for transport.
 
 * This function checks to make sure that all fundamental IMC transport
 * objects are set in the transporter.  If everything is ready a value of
 * true is returned; otherwise, ready returns false.
 
 */
template<class MT, class FT, class PT>
bool DD_Transporter<MT,FT,PT>::ready() const
{
    Require (topology);

    bool indicator = true;

    // if any of the fundamental objects does not exist we are not ready
    if (!mesh) 
	indicator = false;
    else if (!opacity)
	indicator = false;
    else if (!mat_state)
	indicator = false;
    else if (!source) 
	indicator = false;
    else if (!tally)
	indicator = false;
    else if (!communicator)
	indicator = false;

    return indicator;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Query to see number of particles run to completion by the
 * transporter. 
 *
 * \return the total, global number of particles run in the transporter
 */
template<class MT, class FT, class PT>
int DD_Transporter<MT,FT,PT>::get_num_run() const
{
    // local number of particles finished
    int num_particles_run = 0;

    // on the host node this number will be the total number run
    if (rtt_c4::node() == 0)
	num_particles_run = num_done;
	
    // sum up and return
    rtt_c4::global_sum(num_particles_run);

    // the num_particles_run should be the number asked for or zero if they
    // haven't been run yet
    Ensure (num_particles_run == num_to_run || num_particles_run == 0);
    return num_particles_run;
}

} // end namespace rtt_imc

#endif                          // rtt_imc_DD_Transporter_t_hh

//---------------------------------------------------------------------------//
//                        end of imc/DD_Transporter.t.hh
//---------------------------------------------------------------------------//
