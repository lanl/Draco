//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Gray_Particle.t.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 29 13:30:27 2002
 * \brief  Gray_Particle class implementation.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Gray_Particle.hh"

namespace rtt_imc
{

//---------------------------------------------------------------------------//
// STATIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate and return the size of the packed particle.
 * 
 * \param dimension spatial dimension of the problem
 * \param control rtt_rng::Rnd_Control object
 */
template<class MT>
int Gray_Particle<MT>::get_packed_particle_size(
    int                         dimension,
    const rtt_rng::Rnd_Control &control)
{
    Require (dimension > 0 && dimension <= 3);

    // determine the size of the packed random number
    int size_rn = control.get_size();
    Check (size_rn > 0);
    
    // calculate size: int for dimension, int for cell, int for size_rn;
    // dimension doubles for position, 3 doubles for omega, 1 double for ew,
    // 1 double for time_left, 1 double for fraction; size_rn characters for
    // random number state
    int size = 3 * sizeof(int) + (dimension + 6) * sizeof(double) + size_rn;
    
    return size;
}

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 *
 * This constructor is used to unpack a particle that has been packed with
 * the Particle::pack() function.
 */
template<class MT>
Gray_Particle<MT>::Gray_Particle(const std::vector<char> &packed)
{
    Require (packed.size() >= 3 * sizeof(int) + 6 * sizeof(double));

    // make an unpacker
    rtt_dsxx::Unpacker u;
    
    // set it
    u.set_buffer(packed.size(), &packed[0]);

    // unpack the spatial dimension of the particle
    int dimension = 0;
    u >> dimension;
    Check (dimension > 0 && dimension <= 3);

    // size the dimension and direction 
    Base::r.resize(dimension);
    Base::omega.resize(3);

    // unpack the position
    for (int i = 0; i < dimension; i++)
	u >> Base::r[i];

    // unpack the rest of the data
    u >> Base::omega[0] >> Base::omega[1] >> Base::omega[2] 
      >> Base::cell >> Base::ew >> Base::time_left
      >> Base::fraction;
    Check (Base::time_left >= 0.0);
    Check (Base::fraction  >= 0.0);
    Check (Base::cell      >  0);
    Check (Base::ew        >= 0.0);

    // get the size of the RN state
    int size_rn = 0;
    u >> size_rn;
    Check (size_rn > 0);

    // make a packed rn vector
    std::vector<char> prn(size_rn);

    // unpack the rn state
    for (int i = 0; i < size_rn; i++)
	u >> prn[i];

    // rebuild the rn state
    Base::random = new Rnd_Type(prn);
    Check (Base::random->get_num() >= 0);
    Check (Base::random->get_id());

    // assign the descriptor and status
    Base::descriptor = Base::UNPACKED;
    Base::alive      = true;
    
    Ensure (Base::status());
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Do a transport step for a Gray_Particle.
 *
 * This function transports a particle until it:
 * - leaves the mesh
 * - goes to census
 * - gets killed
 * .
 * The optional arguments specify transport options.  If no optional
 * arguments are given then simple transport is performed as described in
 * \ref trans_method.  If random_walk is active then Random Walk is turned on
 * in thick cells as described in \ref trans_method_walk.
 *
 * \subsection trans_method Straight Transport Method
 *
 * Particles undergo implicit absorption until their energy weight drops
 * below 0.01 of their original value. During this time their energy weight
 * decreases exponentialy. Once they attain the cutoff energy weight they are
 * explicitly absorbed and their energy weight remains constant.
 *
 * We unify some of the treatment of the two modes (analog and implicit
 * absorption) by defining sigma_analog_abs to be the actual effective
 * absorption for analog (light) particles and 0.0 for implicit (heavy)
 * particles. The total collision cross section is always
 * sigma_scatter+sigma_eff_abs.
 *
 * To prevent the accumulation of particles with energy weights below the
 * cutoff point, particles stream no further than required to reach cutoff if
 * nothing else happens first. The distance to the cutoff is computed from
 * the current fractional energy weight and the effective absorption in the
 * current cell.
 *
 * \subsection trans_method_walk Transport Method with Random Walk
 *
 * Particles undergo implicit absorption until their energy weight drops
 * below 0.01 of their original value. During this time their energy weight
 * decreases exponentialy. Once they attain the cutoff energy weight they are
 * explicitly absorbed and their energy weight remains constant.
 *
 * We unify some of the treatment of the two modes (analog and implicit
 * absorption) by defining sigma_analog_abs to be the actual effective
 * absorption for analog (light) particles and 0.0 for implicit (heavy)
 * particles. The total collision cross section is always
 * sigma_scatter+sigma_eff_abs.
 *
 * To prevent the accumulation of particles with energy weights below the
 * cutoff point, particles stream no further than required to reach cutoff if
 * nothing else happens first. The distance to the cutoff is computed from
 * the current fractional energy weight and the effective absorption in the
 * current cell.
 *
 * Random walk is used as a hybrid diffusion scheme to speed up the
 * transport.  The conditions for Random Walk are:
 * \f[
 * R_{0} > \sigma_{R}
 * \f]
 * \f[
 * R_{0} > d_{\mathrm{collision}} 
 * \f]
 * and 
 * \f[
 * R_{0} < d_{\mathrm{census}}
 * \f]
 * 
 * The Random Walk algorithm is described in Fleck and Canfield, \e J. \e
 * Comp. \e Phys., \b 54, 1984.  
 *
 * \param mesh the mesh object (determined by MT)
 * \param xs the opacity object
 * \param tally the tally object
 * \param random_walk (optional) rtt_dsxx::SP to a Random_Walk object; if the
 * SP is active then random walk is on
 * \param diagnostic (optional) rtt_dsxx::SP to a diagnostic object; if the
 * SP is active then particle tracks will get dumped to standard out.
 */
template<class MT>
void Gray_Particle<MT>::transport(
    const MT                         &mesh, 
    const Opacity<MT,Gray_Frequency> &xs, 
    Tally<MT>                        &tally, 
    SP_Random_Walk                    random_walk,
    SP_Diagnostic                     diagnostic)
{
    Require (Base::cell > 0);
    Require (Base::cell <= mesh.num_cells());
    Require (tally.num_cells() == mesh.num_cells());

    // >>> RUN THE TRANSPORT

    // run transport with random walk
    if (random_walk)
	rw_transport(mesh, xs, tally, random_walk, diagnostic);

    // otherwise run simple transport
    else
	straight_transport(mesh, xs, tally, diagnostic);

    Ensure (!Base::alive);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pack a particle into a char stream for communication and
 * persistence. 
 */
template<class MT>
std::vector<char> Gray_Particle<MT>::pack() const
{
    Require (Base::omega.size() == 3);

    // make a packer
    rtt_dsxx::Packer p;

    // first pack the random number state
    std::vector<char> prn = Base::random->pack();

    // determine the size of the packed particle: 1 int for cell, + 1 int for
    // size of packed RN state + 1 int for dimension of space; dimension +
    // 6 doubles; + size of RN state chars
    int size = 3 * sizeof(int) + (Base::r.size() + 6) * sizeof(double) +
	prn.size();

    // set the packed buffer
    std::vector<char> packed(size);
    p.set_buffer(size, &packed[0]);

    // pack the spatial dimension
    p << static_cast<int>(Base::r.size());
    
    // pack the dimension
    for (int i = 0; i < Base::r.size(); i++)
	p << Base::r[i];
    
    // pack the rest of the data
    p << Base::omega[0] << Base::omega[1] << Base::omega[2]
      << Base::cell << Base::ew << Base::time_left
      << Base::fraction;

    // pack the RN state
    p << static_cast<int>(prn.size());
    for (int i = 0; i < prn.size(); i++)
	p << prn[i];

    Ensure (p.get_ptr() == &packed[0] + size);
    return packed;
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Do a transport step for the Gray_Particle.
 *
 * This transports a particle until it leaves the mesh, goes to census, or
 * gets killed.
 */
template<class MT>
void Gray_Particle<MT>::straight_transport(
    const MT                         &mesh, 
    const Opacity<MT,Gray_Frequency> &xs, 
    Tally<MT>                        &tally, 
    SP_Diagnostic                     diagnostic)
{
    Require (Base::alive);

    // initialize diagnostics
    if (diagnostic)
    {
	diagnostic->header();
	diagnostic->print(*this);
    }

    // distance to collision, boundary, census and cutoff definitions
    double d_collide, d_boundary, d_census, d_cutoff;
    
    // streaming distance definition
    double d_stream;

    // cell face definition
    int face = 0;

    // event probability definitions
    double prob_thomson_scatter, prob_scatter, prob_abs;

    // intermediate cross section definitions:
    double sigma_thomson_scatter, sigma_eff_scatter, sigma_eff_abs, 
	sigma_scatter, sigma_analog_abs, sigma_collide;
  
    // !!! BEGIN TRANSPORT LOOP !!!

    // transport loop, ended when alive = false
    while (Base::alive)
    {
	// initialize cell face
	face = 0; 

	// get cross sections
	sigma_thomson_scatter = xs.get_sigma_thomson(Base::cell);
	sigma_eff_scatter     = xs.get_sigeffscat(Base::cell);
	sigma_eff_abs         = xs.get_sigeffabs(Base::cell);

	// calculate total scattering cross section
	sigma_scatter = sigma_thomson_scatter + sigma_eff_scatter;

	// determine analog absorption cross section
	if (Base::use_analog_absorption())
 	    sigma_analog_abs = sigma_eff_abs;
	else
	    sigma_analog_abs = 0.0;

	// determine total collision cross section
	sigma_collide = sigma_scatter + sigma_analog_abs;
	Check(sigma_collide>=0);

	// accumulate momentum deposition from volume emission particles
	if (Base::descriptor == Base::VOL_EMISSION)
	    tally.accumulate_momentum(Base::cell, -Base::ew, Base::omega);
        
	// sample distance-to-scatter/absorption (effective scatter or hardball)
	if (sigma_collide == 0 ) 
	{
	    d_collide = rtt_mc::global::huge;    

	    prob_thomson_scatter = 0.0;
	    prob_scatter         = 0.0;
	    prob_abs             = 0.0;
	}
	else 
	{
	    d_collide = -std::log(Base::random->ran()) / sigma_collide;

	    prob_thomson_scatter = sigma_thomson_scatter / sigma_collide;
	    prob_scatter         = sigma_scatter         / sigma_collide;
	    prob_abs             = sigma_analog_abs      / sigma_collide;
	}

	Check(d_collide > 0);

	// get distance-to-boundary and cell face
	d_boundary = mesh.get_db(Base::r, Base::omega, Base::cell, face);  
	Check(d_boundary>=0);

	// distance to census (end of time step)
	d_census = rtt_mc::global::c * Base::time_left;   
	Check(d_census);

	// distance until cutoff weight is reached:
	if (sigma_eff_abs == 0 || Base::use_analog_absorption() )
	{ 
	    d_cutoff = rtt_mc::global::huge;
	}
	else
	{
	    d_cutoff = std::log(Base::fraction / Base::minwt_frac) / 
		sigma_eff_abs;
	}

	Check(d_cutoff > 0);

	// detailed diagnostics
	if (diagnostic)
	    if (diagnostic->detail_status())
	    {
		diagnostic->print_dist(d_collide, d_boundary, d_census, 
				       Base::cell); 
		diagnostic->print_xs(xs, Base:: cell);
	    }

	// determine limiting event
	if      (d_collide < d_boundary  &&  d_collide < d_census   &&
		 d_collide < d_cutoff  )
	{
	    Base::descriptor = Base::COLLISION;
	    d_stream         = d_collide;
	}	
	else if (d_boundary < d_collide  &&  d_boundary < d_census  &&
		 d_boundary < d_cutoff )
	{
	    Base::descriptor = Base::BOUNDARY;
	    d_stream         = d_boundary;
	}
	else if (d_census < d_collide    &&  d_census < d_boundary  &&  
		 d_census < d_cutoff   )
	{
	    Base::descriptor = Base::CENSUS;
	    d_stream         = d_census;
	}
	else if (d_cutoff < d_collide    &&  d_cutoff < d_boundary  &&
		 d_cutoff < d_census   )
	{
	    Base::descriptor = Base::CUTOFF;
	    d_stream         = d_cutoff;
	}
	else
	{
	    throw rtt_dsxx::assertion(
		"Transport could not decide limiting event!");
	}


	// Stream the particle, according to its status:
	if (Base::use_analog_absorption())
	{
	    // Light particle (analog) streaming.
	    Base::stream_analog_capture(tally, d_stream);          
	}
	else
	{
	    // Heavy particle (implicit) streaming
	    stream_implicit_capture(xs, tally, d_stream);    
	}

	// Adjust the time remaining till the end of the time step
	Base::time_left -= d_stream / rtt_mc::global::c;

	// Process collisions, boundary crossings, going to census or
	// reaching cutoff events.
	switch (Base::descriptor) 
	{

	case Base::COLLISION:

	    // process the collision
	    collision_event(mesh, tally, prob_scatter, prob_thomson_scatter, 
			    prob_abs);
	    break;

	case Base::CUTOFF:

	    Check(rtt_mc::global::soft_equiv(Base::fraction, 
					     Base::minwt_frac));

	    // Ensure light weight from now on:
	    Base::fraction = Base::minwt_frac * 0.5;  
	    break;

	case Base::BOUNDARY:

	    // process a boundary event
	    Base::boundary_event(mesh, tally, face);
	    break;

	case Base::CENSUS:
	    
	    // process a census event
	    Base::census_event(tally);
	    break;

	default:

	    // throw an assertion
	    throw rtt_dsxx::assertion("Undefined event in Gray_Particle.");
	}
    } 

    // !!! END OF TRANSPORT LOOP !!!
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do a transport step for the Gray_Particle with Random Walk.
 *
 * This transports a particle until it leaves the mesh, goes to census, or
 * gets killed. Random Walk is used in thick, diffusive cells.
 */
template<class MT>
void Gray_Particle<MT>::rw_transport(
    const MT                         &mesh, 
    const Opacity<MT,Gray_Frequency> &xs, 
    Tally<MT>                        &tally, 
    SP_Random_Walk                    random_walk,
    SP_Diagnostic                     diagnostic)
{
    Require (Base::alive);

    // initialize diagnostics
    if (diagnostic)
    {
	diagnostic->header();
	diagnostic->print(*this);
    }

    // distance to collision, boundary, census and cutoff definitions
    double d_collide, d_boundary, d_census, d_cutoff;
    
    // streaming distance definition
    double d_stream;

    // random walk sphere radius
    double rw_radius;

    // cell face definition
    int face = 0;

    // event probability definitions
    double prob_thomson_scatter, prob_scatter, prob_abs;

    // intermediate cross section definitions:
    double sigma_thomson_scatter, sigma_eff_scatter, sigma_eff_abs, 
	sigma_scatter, sigma_rosseland, sigma_analog_abs, sigma_collide;

    // are we doing random walk on this step
    bool do_a_random_walk;

    // time taken during a random walk
    double rw_time;
  
    // !!! BEGIN TRANSPORT LOOP !!!

    // transport loop, ended when alive = false
    while (Base::alive)
    {
	// initialize cell face
	face = 0; 

	// get cross sections
	sigma_thomson_scatter = xs.get_sigma_thomson(Base::cell);
	sigma_eff_scatter     = xs.get_sigeffscat(Base::cell);
	sigma_eff_abs         = xs.get_sigeffabs(Base::cell);

	// calculate total scattering cross section
	sigma_scatter = sigma_thomson_scatter + sigma_eff_scatter;

	// determine analog absorption cross section
	if (Base::use_analog_absorption())
 	    sigma_analog_abs = sigma_eff_abs;
	else
	    sigma_analog_abs = 0.0;

	// determine total collision cross section
	sigma_collide = sigma_scatter + sigma_analog_abs;
	Check(sigma_collide>=0);

	// accumulate momentum deposition from volume emission particles
	if (Base::descriptor == Base::VOL_EMISSION)
	    tally.accumulate_momentum(Base::cell, -Base::ew, Base::omega);
        
	// sample distance-to-scatter/absorption (effective scatter or hardball)
	if (sigma_collide == 0 ) 
	{
	    d_collide = rtt_mc::global::huge;    

	    prob_thomson_scatter = 0.0;
	    prob_scatter         = 0.0;
	    prob_abs             = 0.0;
	}
	else 
	{
	    d_collide = -std::log(Base::random->ran()) / sigma_collide;

	    prob_thomson_scatter = sigma_thomson_scatter / sigma_collide;
	    prob_scatter         = sigma_scatter         / sigma_collide;
	    prob_abs             = sigma_analog_abs      / sigma_collide;
	}
	Check(d_collide > 0);

	// distance to census (end of time step)
	d_census = rtt_mc::global::c * Base::time_left;   
	Check(d_census);

	// distance until cutoff weight is reached:
	if (sigma_eff_abs == 0 || Base::use_analog_absorption() )
	{ 
	    d_cutoff = rtt_mc::global::huge;
	}
	else
	{
	    d_cutoff = std::log(Base::fraction / Base::minwt_frac) / 
		sigma_eff_abs;
	}
	Check(d_cutoff > 0);

	// check to see if we should do a random walk; check on initial
	// conditions first
	if (Base::descriptor == Base::VOL_EMISSION || 
	    Base::descriptor == Base::EFF_SCATTER)
	{
	    // get the size of the random walk sphere
	    rw_radius = mesh.get_orthogonal_dist_to_bnd(r, cell);
	    Check (rw_radius >= 0.0);

	    // check to see if the random walk conditions are valid
	    do_a_random_walk = random_walk->do_a_random_walk(
		Base::cell, rw_radius, d_collide, d_census);
	}
	else
	{
	    // initial conditions not met, don't do a random walk
	    do_a_random_walk = false;
	}
	 
	// calculate distance to boundary if we are not doing random walk
	if (!do_a_random_walk)
	{
	    // get distance-to-boundary and cell face
	    d_boundary = mesh.get_db(Base::r, Base::omega, Base::cell, face);  
	    Check (d_boundary >= 0);
	}
	else   
	{
	    // we don't need distance to boundary so set it to zero
	    d_boundary = 0.0;
	}

	// detailed diagnostics
	if (diagnostic)
	    if (diagnostic->detail_status())
	    {
		diagnostic->print_dist(d_collide, d_boundary, d_census, 
				       Base::cell); 
		diagnostic->print_xs(xs, Base:: cell);
	    }

	// determine limiting event
	if      (do_a_random_walk)
	{
	    Base::descriptor = Base::RANDOM_WALK;
	    d_stream         = 0.0;
	}
	else if (d_collide < d_boundary  &&  d_collide < d_census   &&
		 d_collide < d_cutoff  )
	{
	    Base::descriptor = Base::COLLISION;
	    d_stream         = d_collide;
	}	
	else if (d_boundary < d_collide  &&  d_boundary < d_census  &&
		 d_boundary < d_cutoff )
	{
	    Base::descriptor = Base::BOUNDARY;
	    d_stream         = d_boundary;
	}
	else if (d_census < d_collide    &&  d_census < d_boundary  &&  
		 d_census < d_cutoff   )
	{
	    Base::descriptor = Base::CENSUS;
	    d_stream         = d_census;
	}
	else if (d_cutoff < d_collide    &&  d_cutoff < d_boundary  &&
		 d_cutoff < d_census   )
	{
	    Base::descriptor = Base::CUTOFF;
	    d_stream         = d_cutoff;
	}
	else
	{
	    throw rtt_dsxx::assertion(
		"Transport could not decide limiting event!");
	}


	// stream the particle: transport or random walk
	if (do_a_random_walk)
	{
	    // boolean to see if particle goes to census
	    bool to_census;

	    // do a random walk
	    rw_time = random_walk->random_walk(Base::r, Base::omega, 
					       Base::time_left, Base::cell, 
					       *Base::random, to_census);

	    // set descriptor if particle goes to census
	    if (to_census)
		Base::descriptor = Base::CENSUS;
	}
	else
	{
	    // Stream the particle, according to its status:
	    if (Base::use_analog_absorption())
	    {
		// Light particle (analog) streaming.
		Base::stream_analog_capture(tally, d_stream);          
	    }
	    else
	    {
		// Heavy particle (implicit) streaming
		stream_implicit_capture(xs, tally, d_stream);    
	    }

	    // Adjust the time remaining till the end of the time step
	    Base::time_left -= d_stream / rtt_mc::global::c;
	}

	// Process collisions, boundary crossings, going to census or
	// reaching cutoff events.
	switch (Base::descriptor) 
	{

	case Base::RANDOM_WALK:

	    // process the random walk absorption
	    random_walk_event(rw_time, tally, xs);
	    break;
	
	case Base::COLLISION:

	    // process the collision
	    collision_event(mesh, tally, prob_scatter, prob_thomson_scatter, 
			    prob_abs);
	    break;

	case Base::CUTOFF:

	    Check(rtt_mc::global::soft_equiv(Base::fraction,
					     Base::minwt_frac));

	    // Ensure light weight from now on:
	    Base::fraction = Base::minwt_frac * 0.5;  
	    break;

	case Base::BOUNDARY:

	    // process a boundary event
	    Base::boundary_event(mesh, tally, face);
	    break;

	case Base::CENSUS:
	    
	    // process a census event
	    Base::census_event(tally);
	    break;

	default:

	    // throw an assertion
	    throw rtt_dsxx::assertion("Undefined event in Gray_Particle.");
	}
    } 

    // !!! END OF TRANSPORT LOOP !!!
}

//---------------------------------------------------------------------------//
/*!
 * \brief Stream for implicit capture.
 */
template<class MT>
void Gray_Particle<MT>::stream_implicit_capture(
    const Opacity<MT,Gray_Frequency> &xs, 
    Tally<MT>                        &tally,
    double                            distance)
{
    Check(distance >= 0.0);
    
    // exponential argument
    double argument = -xs.get_sigeffabs(Base::cell) * distance;

    // calculate multiplicative reduction in energy-weight
    double factor = std::exp(argument);

    // calculate new energy weight; change in energy-weight
    double new_ew = Base::ew * factor;
    double del_ew = Base::ew - new_ew;

    // tally deposited energy and momentum
    tally.deposit_energy(Base::cell, del_ew);
    tally.accumulate_momentum(Base::cell, del_ew, Base::omega);

    // accumulate tallies for energy-weighted path length 
    //     ewpl == energy-weighted-path-length = 
    //             int_0^d e^(-sig*x) dx =
    //             (1/sig)ew(1-e^(-sig*x)), 
    //             or if sig=0, 
    //             ewpl = ew*d.
    if (xs.get_sigeffabs(Base::cell) > 0)
    {
	// integrate exponential from 0 to distance
	tally.accumulate_ewpl(Base::cell, del_ew / 
			      xs.get_sigeffabs(Base::cell));
    }
    else if (xs.get_sigeffabs(Base::cell) == 0)
    {
	// integrate constant
	tally.accumulate_ewpl(Base::cell, distance * Base::ew);
    }
    else if (xs.get_sigeffabs(Base::cell) < 0)
    {
	Insist (0, "Effective absorption is negative!");
    }

    // update the fraction of the particle's original weight
    Base::fraction *= factor;
    Check (Base::fraction > Base::minwt_frac || 
	   rtt_dsxx::soft_equiv(Base::fraction, Base::minwt_frac));

    // update particle energy-weight
    Base::ew = new_ew;

    Check(Base::ew > 0.0);

    // Physically transport the particle
    Base::stream(distance); 
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Process a particle that has undergone random walk.
 */
template<class MT>
void Gray_Particle<MT>::random_walk_event(
    double                            rw_time,
    Tally<MT>                        &tally,
    const Opacity<MT,Gray_Frequency> &opacity)
{
    // adjust weight of random walk particle
    double exponent = -rtt_mc::global::c * 
	opacity.get_sigeffscat(Base::cell) * rw_time;

    // calculate weight factor
    double weight_factor = 0.0;

    // check the exponent against the minimum allowed
    if (exponent > std::numeric_limits<double>::min_exponent)
	weight_factor = std::exp(exponent);

    // adjust weight
    double new_ew   = weight_factor * Base::ew;
    double delta_ew = Base::ew - new_ew;
    Check (delta_ew >= 0.0);

    // do momentum deposition
    tally.accumulate_momentum(Base::cell, delta_ew, Base::omega);
	    
    // do energy deposition
    tally.deposit_energy(Base::cell, delta_ew);

    // tally random walk event
    tally.accum_n_random_walks();

    // update particle energy weight
    Base::ew = new_ew;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process a collision.
 */
template<class MT> 
void Gray_Particle<MT>::collision_event(
    const MT  &mesh, 
    Tally<MT> &tally, 
    double     prob_scatter, 
    double     prob_thomson_scatter, 
    double     prob_abs)
{
    Check (mesh.num_cells() == tally.num_cells());

    // get a random number
    double rand_selector = Base::random->ran();
    
    if (rand_selector < prob_scatter) 
    { 
	// Scatter
	if (rand_selector < prob_thomson_scatter)
	{ 
	    // Thomson scatter
	    Base::descriptor = Base::THOM_SCATTER;
	    tally.accum_n_thomscat();
	}
	else
	{ 
	    // Effective scatter
	    Base::descriptor = Base::EFF_SCATTER;
	    tally.accum_n_effscat();
	}
	
	// accumulate momentum from before the scatter
	tally.accumulate_momentum(Base::cell, Base::ew, Base::omega);
	
	// scatter the particle -- update direction cosines (we really should
	// call effective_scatter for effective scatters because the particle
	// is re-emitted so the direction should be sampled from rest, since
	// everything is isotropic we will wait to do this later (so as not
	// to hose all our regression tests); it has no effect on the
	// practical outcome of the transport
	Base::scatter(mesh);
	
	// accumulate momentum from after the scatter
	tally.accumulate_momentum(Base::cell, -Base::ew, Base::omega);
	
    }
    else if (rand_selector < prob_scatter + prob_abs)
    { 
	// Absorption
	
	// tally absorption data
	tally.deposit_energy(Base::cell, Base::ew);
	tally.accum_n_killed();
	tally.accum_ew_killed(Base::ew);
	tally.accumulate_momentum(Base::cell, Base::ew, Base::omega);

	// set the descriptor and particle status
	Base::descriptor = Base::KILLED; 
	Base::alive      = false;
    }
    else
    {
	Insist(0,"D'oh! Transport could not pick a random event!");
    }   
}

//---------------------------------------------------------------------------//
// DIAGNOSTIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Test equality of two Gray_Particle objects.
 *
 * Test Gray_Particle equality. With regard to the random number state of the
 * particle, this simply checks to see if the Sprng random number object has
 * the same streamnum, not ID pointer.
 */
template<class MT>
bool Gray_Particle<MT>::operator==(const Gray_Particle<MT> &rhs) const
{
    // check particle data
    if (Base::ew != rhs.Base::ew)
	return false;
    else if (Base::r != rhs.Base::r)
	return false;
    else if (Base::omega != rhs.Base::omega)
	return false;
    else if (Base::cell != rhs.Base::cell)
	return false;
    else if (Base::time_left != rhs.Base::time_left)
	return false;
    else if (Base::fraction != rhs.Base::fraction)
	return false;
    else if (Base::alive != rhs.Base::alive)
	return false;
    else if (Base::descriptor != rhs.Base::descriptor)
	return false;

    if (Base::random->get_num() != rhs.Base::random->get_num())
	return false;

    // if all these things check out then the particles are equal
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Test inequality of two Gray_Particle objects.
 */
template<class MT>
bool Gray_Particle<MT>::operator!=(const Gray_Particle<MT> &rhs) const
{
    return !(*this == rhs);
}

//===========================================================================//
// GRAY_PARTICLE::DIAGNOSTIC FUNCTIONS
//===========================================================================//

template<class MT>
void Gray_Particle<MT>::Diagnostic::print_xs(
    const Opacity<MT,Gray_Frequency> &xs,
    int                               cell_in) const
{
    using std::setw;
    using std::endl;
    using std::ios;
    using std::setiosflags;

    // do detailed diagnostic print of particle event cross sections
    Base::Diagnostic::output << setw(20) << setiosflags(ios::right) 
			     << "Opacity: " << setw(12) 
			     << xs.get_sigma_abs(cell_in) << endl;
    Base::Diagnostic::output << setw(20) << setiosflags(ios::right) 
			     << "Eff. scatter: " << setw(12) 
			     << xs.get_sigeffscat(cell_in) << endl; 
    Base::Diagnostic::output << setw(20) << setiosflags(ios::right) 
			     << "Eff. absorption: " << setw(12) 
			     << xs.get_sigeffabs(cell_in) << endl; 
    Base::Diagnostic::output << setw(20) << setiosflags(ios::right) 
			     << "Thomson scatter: " << setw(12) 
			     << xs.get_sigma_thomson(cell_in) << endl; 
}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                        end of imc/Gray_Particle.t.hh
//---------------------------------------------------------------------------//
