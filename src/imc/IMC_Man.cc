//----------------------------------*-C++-*----------------------------------//
// IMC_Man.cc
// Thomas M. Evans
// Wed Jun  3 10:36:11 1998
//---------------------------------------------------------------------------//
// @> IMC_Man class implementation file.
//---------------------------------------------------------------------------//

#include "imc/IMC_Man.hh"
#include "ds++/Assert.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

IMCSPACE

// draco components
using C4::node;
using C4::nodes;
using C4::HTSyncSpinLock;

// stl components
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::ios;
using std::setiosflags;
using std::setw;
using std::fill;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// default constructor

template<class MT, class BT, class IT, class PT>
IMC_Man<MT,BT,IT,PT>::IMC_Man(bool verbose_)
    : delta_t(0), cycle(0), max_cycle(1), print_f(0), verbose(verbose_)
{
  // all the SPs should be null defined

  // objects used by all processors
    Check (!mesh);
    Check (!opacity);
    Check (!mat_state);
    Check (!rnd_con);
    Check (!parallel_builder);
    Check (!buffer);
    Check (!source);
    Check (!tally);

  // objects used only by the host
    Check (!source_init);
    Check (!global_state);

  // objects used to control the census
    Check (!new_census);
}

//---------------------------------------------------------------------------//
// execute IMC
//---------------------------------------------------------------------------//
// run IMC from beginning to end

template<class MT, class BT, class IT, class PT>
void IMC_Man<MT,BT,IT,PT>::execute_IMC(char *argv)
{
  // run through IMC timesteps until we have reached the last cycle
    while (cycle < max_cycle)
    {
      // initialize IMC on the host processor
	host_init(argv);

      // initialize the IMC processors
	IMC_init();

      // do a time_step
	step_IMC();

      // regroup all the results on the host
	regroup();
    }
}

//---------------------------------------------------------------------------//
// host processor initialization
//---------------------------------------------------------------------------//
// here we build the Mesh, Opacity, Mat_State, and Source_Init on the host
// processor 

template<class MT, class BT, class IT, class PT>
void IMC_Man<MT,BT,IT,PT>::host_init(char *argv)
{
  // update the cycle
    cycle++;

  // run only on the host node
    if (node())
	return;

  // build the objects on the host processor
    if (verbose)
	cout << ">> Running through problem builders on proc " << node()
	     << endl;
    
  // run the interface (parser to input)
    string infile = argv;
    SP<IT> interface = new IT(infile);
    interface->parser();
    if (verbose)
	cout << " ** Read input file " << argv << " on node " << node()
	     << endl;

  // make the rnd_control object, set the buffer size, and assign delta_t on
  // the first cycle
    if (!rnd_con)
    {
	delta_t   = interface->get_delta_t();
	max_cycle = interface->get_max_cycle();
	print_f   = interface->get_printf();
	rnd_con   = new Rnd_Control(9836592);
	Particle_Buffer<PT>::set_buffer_size(interface->get_buffer());
    }

  // initialize the mesh builder and build mesh
    {
	BT mt_builder(interface);
	mesh = mt_builder.build_Mesh();
	if (verbose)
	{
	    cout << " ** Built mesh on node " << node() << endl;
	    cout << *mesh << endl;
	}
    }

  // initialize the opacity builder and build the state
    {
	Opacity_Builder<MT>  opacity_builder(interface);

      // build the mat_state
	mat_state = opacity_builder.build_Mat(mesh);

      // update the Mat_State from our last cycle
	if (global_state) 
	    global_state->update_Mat(*mat_state);

      // now build the opacity
	opacity = opacity_builder.build_Opacity(mesh, mat_state);	
	if (verbose)
	{
	    cout << " ** Built Mat_State and Opacity on node " << node() 
		 << endl;
	    cout << *mat_state << endl;
	    cout << *opacity << endl;
	}
    }

  // prepare to do the source initialization
    source_init = new Source_Init<MT>(interface, mesh);

  // update the Source with census info from the last cycle
    if (global_state)
	global_state->update_Source_Init(*source_init);

  // do source initialization
    source_init->initialize(mesh, opacity, mat_state, rnd_con, cycle); 
    if (verbose)
	cout << " ** Did the source initialization on node " << node()
	     << endl;
    
  // make a Global_Buffer to hold T, Cv, evol_net, etc, for all time
    if (!global_state)
	global_state = new Global_Buffer<MT>(*mesh, *mat_state,
					     *source_init); 
}

//---------------------------------------------------------------------------//
// IMC processor initialization
//---------------------------------------------------------------------------//
// here we make communication objects and send stuff out to the IMC problem
// processors

template<class MT, class BT, class IT, class PT>
void IMC_Man<MT,BT,IT,PT>::IMC_init()
{
  // first lets do stuff on the host node
    if (!node())
    {
      // make sure we have performed the host_init
	Require (mesh);
	Require (opacity);
	Require (mat_state);
	Require (rnd_con);
	Require (source_init);
	Require (global_state);

      // make the Parallel Builder and Particle Buffer on the first cycle
	if (!parallel_builder)
	{
	    parallel_builder = new Parallel_Builder<MT>(*mesh, *source_init);
	    buffer           = new Particle_Buffer<PT>(*mesh, *rnd_con);

	  // make sure the global mesh is perserved here
	    Check (parallel_builder->num_cells() ==
		   global_state->num_cells() && mesh->num_cells() ==
		   parallel_builder->num_cells());  

	  // send the buffer sizes, Rnd seed, and delta_t
	    for (int i = 1; i < nodes(); i++)
	    {
		Send (rnd_con->get_seed(), i, 200);
		Send (Particle_Buffer<PT>::get_buffer_s(), i, 201);
		Send (buffer->get_dsize(), i, 202);	
		Send (buffer->get_isize(), i, 203);
		Send (buffer->get_csize(), i, 204);
		Send (delta_t, i, 205);
		Send (max_cycle, i, 206);
		Send (print_f, i, 207);
	    }
	}

      // now send out the mesh, Opacity (full rep right now)
	parallel_builder->send_Mesh(*mesh);
	parallel_builder->send_Opacity(*opacity);

      // send out the Mat_State
	mat_state = parallel_builder->send_Mat(mesh, *mat_state);

      // send out the Source
	source = parallel_builder->send_Source(mesh, mat_state, rnd_con,
					       *source_init, *buffer);

      // make a tally
	tally = new Tally<MT>(mesh);

      // kill objects we no longer need
	kill(source_init);

      // make sure the Global_State census is empty now that it has been sent 
      // out to the IMC processors
	Ensure (global_state->get_census()->size() == 0);
    }

  // lets receive stuff on the IMC nodes
    if (node())
    {
      // make the Rnd_Control, Parallel_Builder, and Particle_Buffer on this
      // processor if we are on the first cycle
	if (!parallel_builder)
	{
	  // receive the data from the host proc the first cycle
	    int seed;
	    int s, d, i, c;
	    Recv (seed, 0, 200);
	    Recv (s, 0, 201);
	    Recv (d, 0, 202);
	    Recv (i, 0, 203);
	    Recv (c, 0, 204);
	    Recv (delta_t, 0, 205);
	    Recv (max_cycle, 0, 206);
	    Recv (print_f, 0, 207);
	    
	  // now make the objects
	    rnd_con = new Rnd_Control(seed);
	    Particle_Buffer<PT>::set_buffer_size(s);
	    parallel_builder = new Parallel_Builder<MT>();
	    buffer           = new Particle_Buffer<PT>(d, i, c);
	}
	
      // receive the objects necessary for transport
	mesh      = parallel_builder->recv_Mesh();
	opacity   = parallel_builder->recv_Opacity(mesh);
	mat_state = parallel_builder->recv_Mat(mesh);
	source    = parallel_builder->recv_Source(mesh, mat_state, rnd_con, 
						  *buffer);
	tally     = new Tally<MT>(mesh);
    }
    
    if (verbose)
	cout << *source << endl;

  // make sure each processor has the requisite objects to do transport
    Ensure (mesh);
    Ensure (opacity);
    Ensure (mat_state);
    Ensure (source);
    Ensure (buffer);
    Ensure (parallel_builder);
    Ensure (rnd_con);
    Ensure (tally);

  // objects we shouldn't have
    Ensure (!source_init);
    Ensure (!new_census);

  // ensure we have times straight
    Ensure (delta_t   > 0);
    Ensure (max_cycle > 0);

  // print a message indicating success
    if (verbose)
	cout << " ** Node " << node() << " has all objects necessary" 
	     << " for IMC transport; initialization complete." << endl;
}

//---------------------------------------------------------------------------//
// IMC transport functions
//---------------------------------------------------------------------------//
// run particles on the problem geometry for one time step

template<class MT, class BT, class IT, class PT>
void IMC_Man<MT,BT,IT,PT>::step_IMC()
{
    cout << ">> Doing particle transport for cycle " << cycle
	 << " on proc " << node() << endl;

  // make a new census Comm_Buffer on this node
    new_census = new Particle_Buffer<PT>::Comm_Buffer();

  // diagnostic
    SP<typename PT::Diagnostic> check = new PT::Diagnostic(cout, true);

  // get source particles and run them to completion
    int counter = 0;
    while (*source)
    {
      // get a particle from the source
	SP<PT> particle = source->get_Source_Particle(delta_t);
	Check (particle->status());

      // transport the particle
	particle->transport(*mesh, *opacity, *tally);
	counter++;

      // after the particle is no longer active take appropriate action
	Check (!particle->status());
	
      // if census write to file
	if (particle->desc() == "census")
	    buffer->buffer_particle(*new_census, *particle);
      // no other choices right now

      // message particle counter
	if (!(counter % print_f)) 
	    cout << setw(10) << counter << " particles run on proc " 
		 << node()   << endl;
    }

  // cout << *tally << endl;
    tally->cycle_print(cout);

  // finished with this timestep
    cout << ">> Finished particle transport for cycle " << cycle
	 << " on proc " << node() << endl;

  // now remove objects we no longer need
    kill (source);
    kill (mat_state);
    kill (opacity);
    kill (mesh);

  // object inventory
    Ensure (rnd_con);
    Ensure (buffer);
    Ensure (parallel_builder);
    Ensure (tally);	
    Ensure (new_census);
    Ensure (!source);
    Ensure (!mat_state);
    Ensure (!opacity);
    Ensure (!source_init);
    Ensure (!mesh);
}

//---------------------------------------------------------------------------//
// concatentate results on the master processor and update the Global-mesh

template<class MT, class BT, class IT, class PT>
void IMC_Man<MT,BT,IT,PT>::regroup()
{
  // send the tallies to the master node
    if (node())
    {
      // calculate the number of cells on this node
	int num_cells = tally->num_cells();

      // make allocatable arrays for tally data
	double *edep_send = new double[num_cells];
	int *ncen_send    = new int[num_cells];

      // assign tally data to arrays
	int ncentot = 0;
	for (int i = 0; i < num_cells; i++)
	{
	    edep_send[i] = tally->get_energy_dep(i+1);
	    ncen_send[i] = tally->get_new_ncen(i+1);
	    ncentot     += ncen_send[i];
	}
	Check (ncentot == tally->get_new_ncen_tot());

      // send tally to master
	Send (num_cells, 0, 300);
	Send (edep_send, num_cells, 0, 301);
	Send (ncen_send, num_cells, 0, 302);

      // send back the Census buffers created on the processor
	buffer->send_buffer(*new_census, 0);

      // release memory
	delete [] edep_send;
	delete [] ncen_send;
	kill (tally);
	kill (new_census);
    }

  // receive the tallies on the master and update the global buffer
    if (!node())
    {
      // accumulated tally data from all processors
	vector<double> accumulate_edep(global_state->num_cells(), 0.0);
	vector<int> accumulate_ncen(global_state->num_cells());
	fill (accumulate_ncen.begin(), accumulate_ncen.end(), 0);

      // census container that is presently in the Global_State
	SP<Particle_Buffer<PT>::Census> census = global_state->get_census();
	Check (census->size() == 0);

      // loop through processors and get stuff
	int ncentot = 0;
	for (int i = 1; i < nodes(); i++)
	{
	    int num_cells;
	    Recv (num_cells, i, 300);
	    Check (num_cells == parallel_builder->num_cells(i));
	    double *edep_recv = new double[num_cells];
	    int *ncen_recv    = new int[num_cells];
	    Recv (edep_recv, num_cells, i, 301);
	    Recv (ncen_recv, num_cells, i, 302);

	  // assign data to accumulate_tally by looping over IMC cells
	    for (int cell = 1; cell <= num_cells; cell++)
	    {
		int global_cell = parallel_builder->master_cell(cell, i);
		accumulate_edep[global_cell-1] += edep_recv[cell-1];
		accumulate_ncen[global_cell-1] += ncen_recv[cell-1];
	    }

	  // reclaim storage
	    delete [] edep_recv;
	    delete [] ncen_recv;
	    
	  // receive the census particles from the remote processors
	    SP<Particle_Buffer<PT>::Comm_Buffer> temp_buffer = 
		buffer->recv_buffer(i);
	    ncentot += temp_buffer->n_part;

	  // write the buffers to the census container
	    buffer->add_to_bank(*temp_buffer, *census);
	}
     
      // add the results from the master processor
	int ncenmaster = 0;
	for (int cell = 1; cell <= tally->num_cells(); cell++)
	{
	    int global_cell = parallel_builder->master_cell(cell, 0);
	    accumulate_edep[global_cell-1] += tally->get_energy_dep(cell);
	    accumulate_ncen[global_cell-1] += tally->get_new_ncen(cell);
	    ncenmaster += tally->get_new_ncen(cell);
	}
	Check (ncenmaster == tally->get_new_ncen_tot());

      // write the master processor census buffer to the census container
	ncentot += new_census->n_part;
	buffer->add_to_bank(*new_census, *census);
	Check (ncentot == census->size());

      // update the global state
	global_state->update_T(accumulate_edep);
	global_state->update_cen(accumulate_ncen);

      // dump this timestep to the screen
	cycle_dump();

      // reclaim objects on the master processor
	kill (tally);
	kill (new_census);
    }

  // do our inventory
    Ensure (!mesh);
    Ensure (!opacity);
    Ensure (!mat_state);
    Ensure (!source);
    Ensure (!tally);
    Ensure (!new_census);
    Ensure (!source_init);
    Ensure (parallel_builder);
    Ensure (buffer);
    Ensure (rnd_con);
}

//---------------------------------------------------------------------------//
// print diagnostics
//---------------------------------------------------------------------------//
// do a screen dump per cycle

template<class MT, class BT, class IT, class PT>
void IMC_Man<MT,BT,IT,PT>::cycle_dump() const
{
    cout << endl;
    cout << ">>> RESULTS FOR CYCLE " << setw(4) << cycle << " <<<" << endl;
    cout << "==============================" << endl;
    cout.precision(4);
    cout << " ** Total time : " << delta_t * cycle << endl;
    cout << *global_state << endl;
}
		
CSPACE

//---------------------------------------------------------------------------//
//                              end of IMC_Man.cc
//---------------------------------------------------------------------------//
