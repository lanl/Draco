//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Gray_Particle.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 29 13:30:27 2002
 * \brief  Gray_Particle class definition.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef RTT_imc_Gray_Particle_HH
#define RTT_imc_Gray_Particle_HH

#include "Particle.hh"
#include "Random_Walk.hh"
#include "ds++/Soft_Equivalence.hh"
#include <cmath>

namespace rtt_imc
{

// Forward declarations
class Gray_Frequency;
 
//===========================================================================//
/*!
 * \class Gray_Particle
 *
 * \brief Particle for doing gray IMC transport.
 */
// revision history:
// -----------------
// 0) original
// 1) 10-FEB-03 : removed #define's for Base class scoping; added real
//                scoping
// 2) 18-FEB-03 : added random walk version of transport function
// 
//===========================================================================//

template<class MT>
class Gray_Particle : public Particle<MT>
{
  public:

    // >>> NESTED TYPES

    /*!
     * \brief Diagnostic class for tracking gray particle histories.
     */
    class Diagnostic : public Particle<MT>::Diagnostic
    {
      public:
	//! Constructor.
	Diagnostic(std::ostream &output_, bool detail_ = false) 
	    : Particle<MT>::Diagnostic(output_, detail_) {}
	
	// Diagnostic print functions.
	void print_xs(const Opacity<MT,Gray_Frequency> &, int) const;
    };
    
    // Friend declarations.
    friend class Diagnostic;

    // Useful typedefs.
    typedef std::vector<double>            sf_double;
    typedef rtt_rng::Sprng                 Rnd_Type;
    typedef rtt_dsxx::SP<Rnd_Type>         SP_Rnd_Type;
    typedef std::string                    std_string;
    typedef rtt_dsxx::SP<Diagnostic>       SP_Diagnostic;
    typedef rtt_dsxx::SP<Random_Walk<MT> > SP_Random_Walk;

  private:
    // Typedef for base class scoping.
    typedef Particle<MT> Base;

  private:
    // >>> IMPLEMENTATION

    // Stream for implicit capture.
    inline void stream_implicit_capture(const Opacity<MT,Gray_Frequency> &,
					Tally<MT> &, double);

    // Process a collision event.
    inline void collision_event(const MT &, Tally<MT> &, double, double, 
				double);

  public:
    // Particle constructor.
    inline Gray_Particle(const sf_double &, const sf_double &, double, int,
			 Rnd_Type, double = 1, double = 1, int = Base::BORN);

    // Unpacking constructor.
    inline Gray_Particle(const std::vector<char> &);

    // >>> TRANSPORT INTERFACE

    // IMC transport step.
    void transport(const MT &, const Opacity<MT,Gray_Frequency> &, 
		   Tally<MT> &, SP_Diagnostic = SP_Diagnostic()); 

    // IMC transport step with Random Walk
    void transport(const MT &, const Opacity<MT,Gray_Frequency> &, 
		   Tally<MT> &, SP_Random_Walk,
		   SP_Diagnostic = SP_Diagnostic()); 

    // >>> DIAGNOSTIC FUNCTIONS

    // Overloaded equals operator.
    bool operator==(const Gray_Particle<MT> &) const;

    //! Overloaded not equals operator.
    bool operator!=(const Gray_Particle<MT> &) const;
 
    // >>> PACKING FUNCTIONALITY
    
    // Pack function
    inline std::vector<char> pack() const;

    // Get the size of the packed particle.
    static int get_packed_particle_size(int, const rtt_rng::Rnd_Control &);
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Gray_Particle constructor.
 *
 * The constructor is declared inline for optimization purposes.
 */
template<class MT>
Gray_Particle<MT>::Gray_Particle(const sf_double& r_, 
				 const sf_double& omega_, 
				 double           ew_,
				 int              cell_, 
				 Rnd_Type         random_, 
				 double           frac,
				 double           tleft, 
				 int              desc)
    : Particle<MT>(r_, omega_, ew_, cell_, random_, frac, tleft, desc)
{
    // non-default particle constructor
    Check (Base::r.size() < 4);
    Check (Base::omega.size() == 3);
    Check (Base::cell > 0); 
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 *
 * This constructor is used to unpack a particle that has been packed with
 * the Particle::pack() function.
 *
 * It is declared inline for optimization.
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

} // end namespace rtt_imc

#endif                          // RTT_imc_Gray_Particle_HH

//---------------------------------------------------------------------------//
//                              end of imc/Gray_Particle.hh
//---------------------------------------------------------------------------//
