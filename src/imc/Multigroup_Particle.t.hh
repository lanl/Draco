//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Multigroup_Particle.t.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 29 13:31:05 2002
 * \brief  Multigroup_Particle implementation.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Multigroup_Particle_t_hh__
#define __imc_Multigroup_Particle_t_hh__

#include "Multigroup_Particle.hh"

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
int Multigroup_Particle<MT>::get_packed_particle_size(
    int                         dimension,
    const rtt_rng::Rnd_Control &control)
{
    Require (dimension > 0 && dimension <= 3);

    // determine the size of the packed random number
    int size_rn = control.get_size();
    Check (size_rn > 0);
    
    // calculate size: int for dimension, int for cell, int for size_rn, int
    // for group_index; dimension doubles for position, 3 doubles for omega,
    // 1 double for ew, 1 double for time_left, 1 double for fraction;
    // size_rn characters for random number state
    int size = 4 * sizeof(int) + (dimension + 6) * sizeof(double) + size_rn;
    
    return size;
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Do a transport step for the Mulitgroup_Particle.
 *
 * This transports a particle until it leaves the mesh, goes to census, or
 * gets killed.
 *
 * \p Transport Method
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
 */
template<class MT>
void Multigroup_Particle<MT>::transport(const MT   &mesh, 
				  const MG_Opacity &xs, 
				  Tally<MT>        &tally, 
				  SP_Diagnostic     diagnostic)
{
    Require (alive);

    // initialize diagnostics
    if (diagnostic)
    {
	diagnostic->header();
	diagnostic->print(*this);
    }
  
    // !!! BEGIN TRANSPORT LOOP !!!

    // transport loop, ended when alive = false
    while (alive)
    {
	// distance to collision, boundary, census and cutoff defnintions
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

	sigma_thomson_scatter = xs.get_sigma_thomson(cell, group_index);
	sigma_eff_scatter     = xs.get_sigeffscat(cell, group_index);
	sigma_eff_abs         = xs.get_sigeffabs(cell, group_index);

	sigma_scatter         = sigma_thomson_scatter + sigma_eff_scatter;

	if (use_analog_absorption())
 	    sigma_analog_abs  = sigma_eff_abs;
	else
	    sigma_analog_abs  = 0.0;

	sigma_collide = sigma_scatter + sigma_analog_abs;

	Check(sigma_collide >= 0);

	// accumulate momentum deposition from volume emission particles
	if (descriptor == VOL_EMISSION)
	    tally.accumulate_momentum(cell, -ew, omega);
        
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
	    d_collide = -log(random->ran()) / sigma_collide;

	    prob_thomson_scatter = sigma_thomson_scatter / sigma_collide;
	    prob_scatter         = sigma_scatter         / sigma_collide;
	    prob_abs             = sigma_analog_abs      / sigma_collide;
	}

	Check (d_collide > 0);

	// get distance-to-boundary and cell face
	d_boundary = mesh.get_db(r, omega, cell, face); 
	Check (d_boundary >= 0);

	// distance to census (end of time step)
	d_census = rtt_mc::global::c * time_left;  
	Check (d_census);

	// distance until cutoff weight is reached:
	if (sigma_eff_abs == 0 || use_analog_absorption() )
	{ 
	    d_cutoff = rtt_mc::global::huge;
	}
	else
	{
	    d_cutoff = log(fraction/minwt_frac) / sigma_eff_abs;
	}

	Check (d_cutoff > 0);


	// detailed diagnostics
	if (diagnostic)
	    if (diagnostic->detail_status())
	    {
		diagnostic->print_dist(d_collide, d_boundary, d_census, cell); 
		diagnostic->print_xs(xs, cell, group_index);
	    }


	// determine limiting event
	if      (d_collide < d_boundary  &&  d_collide < d_census   &&
		 d_collide < d_cutoff  )
	{
	    descriptor = COLLISION;
	    d_stream   = d_collide;
	}	
	else if (d_boundary < d_collide  &&  d_boundary < d_census  &&
		 d_boundary < d_cutoff )
	{
	    descriptor = BOUNDARY;
	    d_stream   = d_boundary;
	}
	else if (d_census < d_collide    &&  d_census < d_boundary  &&  
		 d_census < d_cutoff   )
	{
	    descriptor = CENSUS;
	    d_stream   = d_census;
	}
	else if (d_cutoff < d_collide    &&  d_cutoff < d_boundary  &&
		 d_cutoff < d_census   )
	{
	    descriptor = CUTOFF;
	    d_stream   = d_cutoff;
	}  
	else
	{
	    Insist(0,"D'oh! Transport could not decide limiting event!");
	}


	// Stream the particle, according to its status:
	if (use_analog_absorption())
	{
	    // Light particle (analog) streaming.
	    stream_analog_capture(tally, d_stream);          
	}
	else
	{
	    // Heavy particle (implicit) streaming
	    stream_implicit_capture(xs, tally, d_stream);    
	}


	// Adjust the time remaining till the end of the time step
	time_left -= d_stream / rtt_mc::global::c;


	// Process collisions, boundary crossings, going to census or
	// reaching cutoff events.
	switch (descriptor) 
	{

	case COLLISION:

	    // process the collision event
	    collision_event(mesh, tally, xs, prob_scatter, 
			    prob_thomson_scatter, prob_abs);
	    break;

	case CUTOFF:

	    Check(rtt_mc::global::soft_equiv(fraction, minwt_frac));

	    // Ensure light weight from now on:
	    fraction = minwt_frac * 0.5;  
	    break;

	case BOUNDARY:

	    // process a boundary event
	    boundary_event(mesh, tally, face);
	    break;

	case CENSUS:
	    
	    // process a census event
	    census_event(tally);
	    break;

	default:

	    // throw an assertion
	    throw 
		rtt_dsxx::assertion("Undefined event in Multigroup_Particle.");
	}

    } 

    // !!! END OF TRANSPORT LOOP !!!
}

//---------------------------------------------------------------------------//
// DIAGNOSTIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Print out a Multigroup_Particle to some stream.
 */
template<class MT>
void Multigroup_Particle<MT>::print(std::ostream &output) const
{
    using std::ios;
    using std::setiosflags;
    using std::endl;
    using std::setw;

    // call the base class print
    Particle<MT>::print(output);

    // add the group index to the list
    output << setw(20) << setiosflags(ios::right) << "Group Index: " 
           << setw(12) << group_index << endl;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Test equality of two Multigroup_Particle objects.
 *
 * Test Multigroup_Particle equality.  With regard to the random number state
 * of the particle, this simply checks to see if the Sprng random number
 * object has the same streamnum, not ID pointer
 */
template<class MT>
bool Multigroup_Particle<MT>::operator==(const Multigroup_Particle<MT> &rhs)
    const
{
    // check particle data
    if (ew != rhs.ew)
	return false;
    else if (r != rhs.r)
	return false;
    else if (omega != rhs.omega)
	return false;
    else if (cell != rhs.cell)
	return false;
    else if (time_left != rhs.time_left)
	return false;
    else if (fraction != rhs.fraction)
	return false;
    else if (alive != rhs.alive)
	return false;
    else if (descriptor != rhs.descriptor)
	return false;
    else if (group_index != rhs.group_index)
	return false;

    if (random->get_num() != rhs.random->get_num())
	return false;

    // if all these things check out then the particles are equal
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Test inequality of two Gray_Particle objects.
 */
template<class MT>
bool Multigroup_Particle<MT>::operator!=(const Multigroup_Particle<MT> &rhs) 
    const
{
    return !(*this == rhs);
}

//===========================================================================//
// MULTIGROUP_PARTICLE::DIAGNOSTIC FUNCTIONS
//===========================================================================//

template<class MT>
void Multigroup_Particle<MT>::Diagnostic::print_xs(
    const MG_Opacity &xs,
    int               cell,
    int               group_index) const
{
    using std::setw;
    using std::endl;
    using std::ios;
    using std::setiosflags;

    // do detailed diagnostic print of particle event cross sections
    output << setw(20) << setiosflags(ios::right) << "Opacity: " 
	   << setw(12) << xs.get_sigma_abs(cell, group_index)  << endl;
    output << setw(20) << setiosflags(ios::right) << "Eff. scatter: "
	   << setw(12) << xs.get_sigeffscat(cell, group_index) << endl; 
    output << setw(20) << setiosflags(ios::right) << "Eff. absorption: " 
	   << setw(12) << xs.get_sigeffabs(cell, group_index)  << endl; 
    output << setw(20) << setiosflags(ios::right) << "Thomson scatter: " 
	   << setw(12) << xs.get_sigma_thomson(cell, group_index)  << endl; 
}

} // end namespace rtt_imc

#endif                          // __imc_Multigroup_Particle_t_hh__

//---------------------------------------------------------------------------//
//                        end of imc/Multigroup_Particle.t.hh
//---------------------------------------------------------------------------//
