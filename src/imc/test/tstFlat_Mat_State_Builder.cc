//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstFlat_Mat_State_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed Mar  5 15:53:10 2003
 * \brief  test Flat_Mat_State_Builder class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "IMC_Test.hh"
#include "../Release.hh"
#include "../Flat_Mat_State_Builder.hh"
#include "../Opacity.hh"
#include "../Mat_State.hh"
#include "../Diffusion_Opacity.hh"
#include "../Frequency.hh"
#include "../Global.hh"
#include "../Particle.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "cdi/CDI.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

using rtt_imc_test::Parser;
using rtt_imc_test::IMC_Flat_Interface;
using rtt_imc::Flat_Mat_State_Builder;
using rtt_imc::Mat_State_Builder;
using rtt_imc::Mat_State;
using rtt_imc::Diffusion_Opacity;
using rtt_imc::Opacity;
using rtt_mc::OS_Builder;
using rtt_mc::OS_Mesh;
using rtt_cdi::CDI;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

typedef rtt_mc::OS_Mesh               MT;
typedef rtt_imc::Gray_Frequency       G;
typedef rtt_imc::Multigroup_Frequency MG;
typedef rtt_imc::Particle<MT>         PT;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void gray_flat_mat_state_test()
{
    // build a parser, mesh, and interface
    SP<Parser>                  parser(new Parser("OS_Input"));
    SP<OS_Builder>              mb(new OS_Builder(parser));
    SP<OS_Mesh>                 mesh = mb->build_Mesh();
    SP<IMC_Flat_Interface<PT> > interface(new IMC_Flat_Interface<PT>(mb));
    
    // pointer to a mat state builder
    SP<Mat_State_Builder<MT,G> > builder;

    // make a flat mat state builder
    builder = new Flat_Mat_State_Builder<MT,G>(interface);

    // check typeinfo
    if (typeid(*builder) != typeid(Flat_Mat_State_Builder<MT,G>))  ITFAILS;
    if (typeid(builder.bp()) != typeid(Mat_State_Builder<MT,G> *)) ITFAILS;
    if (typeid(*builder) == typeid(Flat_Mat_State_Builder<MT,MG>)) ITFAILS;

    // build objects
    builder->build_mat_classes(mesh);

    // build the frequency
    SP<G> frequency = builder->get_Frequency();
    
    // brief check
    if (!frequency->is_gray()) ITFAILS;

    // make a Mat_State
    SP<Mat_State<MT> > mat_state = builder->get_Mat_State();

    // check the mat_state
    if (mat_state->num_cells() != 6) ITFAILS;

    for (int cell = 1; cell <= 3; cell++)
    {
	double d        = mat_state->get_rho(cell);
	double T        = mat_state->get_T(cell);
	double Cv       = mat_state->get_spec_heat(cell);
	double dedT     = mat_state->get_dedt(cell);
	double dedT_ref = 0.1 * 1.0 * mesh->volume(cell);

	if (!soft_equiv(d, 1.0))         ITFAILS;
	if (!soft_equiv(T, 10.0))        ITFAILS;
	if (!soft_equiv(Cv, 0.1))        ITFAILS;
	if (!soft_equiv(dedT, dedT_ref)) ITFAILS;
    }
    for (int cell = 4; cell <= 6; cell++)
    {
	double d        = mat_state->get_rho(cell);
	double T        = mat_state->get_T(cell);
	double Cv       = mat_state->get_spec_heat(cell);
	double dedT     = mat_state->get_dedt(cell);
	double dedT_ref = 0.2 * 2.0 * mesh->volume(cell);

	if (!soft_equiv(d, 2.0))         ITFAILS;
	if (!soft_equiv(T, 20.0))        ITFAILS;
	if (!soft_equiv(Cv, 0.2))        ITFAILS;
	if (!soft_equiv(dedT, dedT_ref)) ITFAILS;
    }

    if (rtt_imc_test::passed)
	PASSMSG("Flat Mat_State passes all tests in Gray test.");

    // make an Opacity
    SP<Opacity<MT,G> > opacity = builder->get_Opacity();

    // check the opacity
    if (opacity->num_cells() != 6) ITFAILS;

    for (int cell = 1; cell <= 3; cell++)
    {
	double sig_abs = opacity->get_sigma_abs(cell);
	double sig_sc  = opacity->get_sigma_thomson(cell);
	double fleck   = opacity->get_fleck(cell);
	double T       = mat_state->get_T(cell);
	double dedT    = mat_state->get_dedt(cell);
	double effabs  = opacity->get_sigeffabs(cell);
	double effsc   = opacity->get_sigeffscat(cell);
	
	double fleck_ref   = 1.0 / 
	    (1.0 + 4.0*rtt_mc::global::a*rtt_mc::global::c*1000.0*
	     mesh->volume(cell)*0.001/dedT*0.1);

	double effabs_ref  = fleck_ref * 0.1;
	double effsc_ref   = (1.0 - fleck_ref) * 0.1;

	if (!soft_equiv(sig_abs, .1))          ITFAILS;
	if (!soft_equiv(sig_sc, .5))           ITFAILS;
	if (!soft_equiv(fleck, fleck_ref))     ITFAILS;
	if (!soft_equiv(effabs, effabs_ref))   ITFAILS;
	if (!soft_equiv(effsc, effsc_ref))     ITFAILS;
    }

    for (int cell = 4; cell <= 6; cell++)
    {
	double sig_abs = opacity->get_sigma_abs(cell);
	double sig_sc  = opacity->get_sigma_thomson(cell);
	double fleck   = opacity->get_fleck(cell);
	double T       = mat_state->get_T(cell);
	double dedT    = mat_state->get_dedt(cell);
	double effabs  = opacity->get_sigeffabs(cell);
	double effsc   = opacity->get_sigeffscat(cell);
	
	double fleck_ref   = 1.0 / 
	    (1.0 + 4.0*rtt_mc::global::a*rtt_mc::global::c*8000.0*
	     mesh->volume(cell)*0.001/dedT*0.02);

	double effabs_ref  = fleck_ref * 0.02;
	double effsc_ref   = (1.0 - fleck_ref) * 0.02;

	if (!soft_equiv(sig_abs, .02))         ITFAILS;
	if (!soft_equiv(sig_sc, 0.0))          ITFAILS;
	if (!soft_equiv(fleck, fleck_ref))     ITFAILS;
	if (!soft_equiv(effabs, effabs_ref))   ITFAILS;
	if (!soft_equiv(effsc, effsc_ref))     ITFAILS;
    }
    
    // we shouldn't have a diffusion opacity
    SP<Diffusion_Opacity<MT> > diff = builder->get_Diffusion_Opacity();
    if (diff) ITFAILS;

    if (rtt_imc_test::passed)
	PASSMSG("Flat Opacity passes all tests in Gray test.");
}

//---------------------------------------------------------------------------//

void gray_flat_diffusion_mat_state_test()
{
    // build a parser, mesh, and interface
    SP<Parser>                  parser(new Parser("OS_Input"));
    SP<OS_Builder>              mb(new OS_Builder(parser));
    SP<OS_Mesh>                 mesh = mb->build_Mesh();
    SP<IMC_Flat_Interface<PT> > interface(new IMC_Flat_Interface<PT>(mb, 0, 1));
    
    // pointer to a mat state builder
    SP<Mat_State_Builder<MT,G> > builder;

    // make a flat mat state builder
    builder = new Flat_Mat_State_Builder<MT,G>(interface);

    // build objects
    builder->build_mat_classes(mesh);

    // get the opacity and mat states
    SP<G>                      frequency = builder->get_Frequency();
    SP<Mat_State<MT> >         mat_state = builder->get_Mat_State();
    SP<Opacity<MT,G> >         opacity   = builder->get_Opacity();
    SP<Diffusion_Opacity<MT> > diff      = builder->get_Diffusion_Opacity();
    if (!diff)                  ITFAILS;
    if (diff->num_cells() != 6) ITFAILS;

    for (int cell = 1; cell <= 3; cell++)
    {
	double sig_abs = opacity->get_sigma_abs(cell);
	double sig_sc  = opacity->get_sigma_thomson(cell);
	double rosref  = sig_abs + sig_sc;
	double dedT    = mat_state->get_dedt(cell);

	double fleck = diff->get_fleck(cell);
	double ros   = diff->get_Rosseland_opacity(cell);

	double fleck_ref   = 1.0 / 
	    (1.0 + 4.0*rtt_mc::global::a*rtt_mc::global::c*1000.0*
	     mesh->volume(cell)*0.001/dedT*0.1);

	if (!soft_equiv(sig_abs, .1))      ITFAILS;
	if (!soft_equiv(sig_sc, .5))       ITFAILS;
	if (!soft_equiv(fleck, fleck_ref)) ITFAILS;
	if (!soft_equiv(ros, rosref))      ITFAILS;
    }

    for (int cell = 4; cell <= 6; cell++)
    {
	double sig_abs = opacity->get_sigma_abs(cell);
	double sig_sc  = opacity->get_sigma_thomson(cell);
	double rosref  = sig_abs + sig_sc;
	double dedT    = mat_state->get_dedt(cell);

	double fleck = diff->get_fleck(cell);
	double ros   = diff->get_Rosseland_opacity(cell);
	
	double fleck_ref   = 1.0 / 
	    (1.0 + 4.0*rtt_mc::global::a*rtt_mc::global::c*8000.0*
	     mesh->volume(cell)*0.001/dedT*0.02);

	if (!soft_equiv(sig_abs, .02))     ITFAILS;
	if (!soft_equiv(sig_sc, 0.0))      ITFAILS;
	if (!soft_equiv(fleck, fleck_ref)) ITFAILS;
	if (!soft_equiv(ros, rosref))      ITFAILS;
    }
    
    if (rtt_imc_test::passed)
	PASSMSG("Flat Diffusion_Opacity passes all tests in Gray test.");
}

//---------------------------------------------------------------------------//

void mg_flat_mat_state_test()
{
    // build a parser, mesh, and interface
    SP<Parser>                    parser(new Parser("OS_Input"));
    SP<OS_Builder>                mb(new OS_Builder(parser));
    SP<OS_Mesh>                   mesh = mb->build_Mesh();
    SP<IMC_Flat_Interface<PT> >   interface(new IMC_Flat_Interface<PT>(mb));
    
    // pointer to a mat state builder
    SP<Mat_State_Builder<MT,MG> > builder;

    // make a CDI mat state builder
    builder = new Flat_Mat_State_Builder<MT,MG>(interface);

    // check typeinfo
    if (typeid(*builder) != typeid(Flat_Mat_State_Builder<MT,MG>))  ITFAILS;
    if (typeid(builder.bp()) != typeid(Mat_State_Builder<MT,MG> *)) ITFAILS;
    if (typeid(*builder) == typeid(Flat_Mat_State_Builder<MT,G>))   ITFAILS;

    // build objects
    builder->build_mat_classes(mesh);

    // build the frequency
    SP<MG> frequency = builder->get_Frequency();
    
    // brief check
    if (!frequency->is_multigroup()) ITFAILS;

    // make a Mat_State
    SP<Mat_State<MT> > mat_state = builder->get_Mat_State();

    // check the mat state
    {
	if (mat_state->num_cells() != 6) ITFAILS;

	if (!soft_equiv(mat_state->get_rho(1), 1.0))       ITFAILS;
	if (!soft_equiv(mat_state->get_rho(2), 1.0))       ITFAILS;
	if (!soft_equiv(mat_state->get_rho(3), 1.0))       ITFAILS;
	if (!soft_equiv(mat_state->get_rho(4), 2.0))       ITFAILS;
	if (!soft_equiv(mat_state->get_rho(5), 2.0))       ITFAILS;
	if (!soft_equiv(mat_state->get_rho(6), 2.0))       ITFAILS;

	if (!soft_equiv(mat_state->get_T(1), 10.0))        ITFAILS;
	if (!soft_equiv(mat_state->get_T(2), 10.0))        ITFAILS;
	if (!soft_equiv(mat_state->get_T(3), 10.0))        ITFAILS;
	if (!soft_equiv(mat_state->get_T(4), 20.0))        ITFAILS;
	if (!soft_equiv(mat_state->get_T(5), 20.0))        ITFAILS;
	if (!soft_equiv(mat_state->get_T(6), 20.0))        ITFAILS;

	if (!soft_equiv(mat_state->get_spec_heat(1), 0.1)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(2), 0.1)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(3), 0.1)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(4), 0.2)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(5), 0.2)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(6), 0.2)) ITFAILS;
    }

    if (rtt_imc_test::passed)
	PASSMSG("Flat Mat_State passes all tests in MG test.");

    // make an Opacity
    SP<Opacity<MT,MG> > opacity = builder->get_Opacity();

    // check the groups
    vector<double> ref_groups(4);
    ref_groups[0] = 0.01;
    ref_groups[1] = 0.1;
    ref_groups[2] = 15.0;
    ref_groups[3] = 100.0;
    {
	vector<double> groups = frequency->get_group_boundaries();
	if (!soft_equiv(groups.begin(), groups.end(), 
			ref_groups.begin(), ref_groups.end())) ITFAILS;
    }

    // check the opacity
    {
	if (!soft_equiv(opacity->get_sigma_abs(1, 1), 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(1, 2), 0.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(1, 3), 0.1)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2, 1), 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2, 2), 0.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2, 3), 0.1)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 1), 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 2), 0.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 3), 0.1)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 1), 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 2), 1.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 3), 1.1)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 1), 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 2), 1.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 3), 1.1)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6, 1), 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6, 2), 1.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6, 3), 1.1)) ITFAILS;

	for (int c = 1; c <= 6; c++)
	    for (int g = 1; g <= frequency->get_num_groups(); g++)
		if (!soft_equiv(opacity->get_sigma_thomson(c,g), 0.0)) 
		    ITFAILS;
	
	// check integrated norm Planck
	double int_Planck;

	int_Planck = CDI::integratePlanckSpectrum(0.01, 100.0, 10.0);
	for (int c = 1; c <= 3; c++)
	    if (!soft_equiv(opacity->get_integrated_norm_Planck(c),
			    int_Planck)) ITFAILS;

	int_Planck = CDI::integratePlanckSpectrum(0.01, 100.0, 20.0);
	for (int c = 4; c <= 6; c++)
	    if (!soft_equiv(opacity->get_integrated_norm_Planck(c),
			    int_Planck)) ITFAILS;

	// check effective cross sections
	for (int c = 1; c <= 6; c++)
	{
	    pair<double,double> b;
	    double              sum = 0.0;
	    vector<double>      ref_emission(3);
	    double              T = mat_state->get_T(c);
	    for (int g = 1; g <= frequency->get_num_groups(); g++)
	    {
		b                 = frequency->get_group_boundaries(g);
		sum              += opacity->get_sigma_abs(c,g) *
		    CDI::integratePlanckSpectrum(b.first, b.second, T);
		ref_emission[g-1] = sum;
	    }
	    
	    vector<double> emission = opacity->get_emission_group_cdf(c);

	    if (!soft_equiv(emission.begin(), emission.end(),
			    ref_emission.begin(), ref_emission.end())) ITFAILS;

	    double planck = sum / opacity->get_integrated_norm_Planck(c);

	    double beta   = 4.0 * rtt_mc::global::a * T*T*T * mesh->volume(c) /
		mat_state->get_dedt(c);

	    double fleck  = 1.0 / 
		(1.0 + .001 * rtt_mc::global::c * beta * planck);  

	    if (!soft_equiv(opacity->get_fleck(c), fleck)) ITFAILS;

	    for (int g = 1; g <= frequency->get_num_groups(); g++)
	    {
		double effabs = opacity->get_sigma_abs(c,g) * fleck;
		double effsct = opacity->get_sigma_abs(c,g) * (1.0-fleck);

		if (!soft_equiv(opacity->get_sigeffscat(c,g),effsct)) ITFAILS;
		if (!soft_equiv(opacity->get_sigeffabs(c,g),effabs))  ITFAILS;
	    }
	}
    }

    if (rtt_imc_test::passed)
	PASSMSG("Flat MG Frequency and MG Opacity passes all tests.");
}

//---------------------------------------------------------------------------//

void mg_flat_diffusion_mat_state_test()
{
    // build a parser, mesh, and interface
    SP<Parser>                  parser(new Parser("OS_Input"));
    SP<OS_Builder>              mb(new OS_Builder(parser));
    SP<OS_Mesh>                 mesh = mb->build_Mesh();

    // build dummy interface with hybrid diffusion on and common mg opacities
    // off
    SP<IMC_Flat_Interface<PT> > interface(new IMC_Flat_Interface<PT>(mb, 0, 1));
    
    // pointer to a mat state builder
    SP<Mat_State_Builder<MT,MG> > builder;

    // make a flat mat state builder
    builder = new Flat_Mat_State_Builder<MT,MG>(interface);

    // build objects
    builder->build_mat_classes(mesh);

    // get the opacity and mat states
    SP<MG>                     frequency = builder->get_Frequency();
    SP<Mat_State<MT> >         mat_state = builder->get_Mat_State();
    SP<Opacity<MT,MG> >        opacity   = builder->get_Opacity();
    SP<Diffusion_Opacity<MT> > diff      = builder->get_Diffusion_Opacity();
    if (!diff)                  ITFAILS;
    if (diff->num_cells() != 6) ITFAILS;

    // check the groups
    vector<double> ref_groups(4);
    ref_groups[0] = 0.01;
    ref_groups[1] = 0.1;
    ref_groups[2] = 15.0;
    ref_groups[3] = 100.0;
    {
	vector<double> groups = frequency->get_group_boundaries();
	if (!soft_equiv(groups.begin(), groups.end(), 
			ref_groups.begin(), ref_groups.end())) ITFAILS;
    }

    // check the opacities
    double int_Planck;
    {
	if (!soft_equiv(opacity->get_sigma_abs(1, 1), 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(1, 2), 0.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(1, 3), 0.1)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2, 1), 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2, 2), 0.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2, 3), 0.1)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 1), 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 2), 0.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 3), 0.1)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 1), 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 2), 1.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 3), 1.1)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 1), 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 2), 1.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 3), 1.1)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6, 1), 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6, 2), 1.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6, 3), 1.1)) ITFAILS;

	for (int c = 1; c <= 6; c++)
	    for (int g = 1; g <= frequency->get_num_groups(); g++)
		if (!soft_equiv(opacity->get_sigma_thomson(c,g), 0.0)) 
		    ITFAILS;

	int_Planck = CDI::integratePlanckSpectrum(0.01, 100.0, 10.0);
	for (int c = 1; c <= 3; c++)
	    if (!soft_equiv(opacity->get_integrated_norm_Planck(c),
			    int_Planck)) ITFAILS;

	int_Planck = CDI::integratePlanckSpectrum(0.01, 100.0, 20.0);
	for (int c = 4; c <= 6; c++)
	    if (!soft_equiv(opacity->get_integrated_norm_Planck(c),
			    int_Planck)) ITFAILS;

	// check effective cross sections in opacity and diffusion opacity
	for (int c = 1; c <= 6; c++)
	{
	    pair<double,double> b;
	    double              sum = 0.0;
	    vector<double>      ref_emission(3);
	    double              T = mat_state->get_T(c);
	    for (int g = 1; g <= frequency->get_num_groups(); g++)
	    {
		b                 = frequency->get_group_boundaries(g);
		sum              += opacity->get_sigma_abs(c,g) *
		    CDI::integratePlanckSpectrum(b.first, b.second, T);
		ref_emission[g-1] = sum;
	    }
	    
	    vector<double> emission = opacity->get_emission_group_cdf(c);

	    if (!soft_equiv(emission.begin(), emission.end(),
			    ref_emission.begin(), ref_emission.end())) ITFAILS;

	    double planck = sum / opacity->get_integrated_norm_Planck(c);

	    double beta   = 4.0 * rtt_mc::global::a * T*T*T * mesh->volume(c) /
		mat_state->get_dedt(c);

	    double fleck  = 1.0 / 
		(1.0 + .001 * rtt_mc::global::c * beta * planck);  

	    if (!soft_equiv(diff->get_fleck(c), fleck)) ITFAILS;

	    for (int g = 1; g <= frequency->get_num_groups(); g++)
	    {
		double effabs = opacity->get_sigma_abs(c,g) * fleck;
		double effsct = opacity->get_sigma_abs(c,g) * (1.0-fleck);

		if (!soft_equiv(opacity->get_sigeffscat(c,g),effsct)) ITFAILS;
		if (!soft_equiv(opacity->get_sigeffabs(c,g),effabs))  ITFAILS;
	    }

	    // check Rosseland opacities in diffusion opacity
	    double ros_sum         = 0.0;
	    double inv_sig_ros_sum = 0.0;
	    double r_g             = 0.0;
	    double ros_ref         = 0.0;
	    for (int g = 1; g <= frequency->get_num_groups(); g++)
	    {
		b   = frequency->get_group_boundaries(g);
		r_g = CDI::integrateRosselandSpectrum(b.first, b.second, T);

		// sums (scattering is zero)
		ros_sum         += r_g;
		inv_sig_ros_sum += r_g / opacity->get_sigma_abs(c, g);
	    }

	    // check ros sum
	    if (!soft_equiv(
		    CDI::integrateRosselandSpectrum(0.01, 100.0, T),
		    ros_sum)) ITFAILS;

	    // calculate reference rosseland and compare
	    ros_ref = ros_sum / inv_sig_ros_sum;

	    // check it
	    if (!soft_equiv(diff->get_Rosseland_opacity(c), ros_ref)) ITFAILS;
	}
    }

    if (rtt_imc_test::passed)
	PASSMSG("Flat Diffusion_Opacity passes all tests in MG test.");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::initialize(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (rtt_c4::node() == 0)
		cout << argv[0] << ": version " << rtt_imc::release() 
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS

	// gray tests
	gray_flat_mat_state_test();
	gray_flat_diffusion_mat_state_test();

	// mg tests
	mg_flat_mat_state_test();
	mg_flat_diffusion_mat_state_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstFlat_Mat_State_Builder, " << ass.what()
	     << endl;
	rtt_c4::finalize();
	return 1;
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_imc_test::passed) 
	{
	    cout << "**** tstFlat_Mat_State_Builder Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    rtt_c4::global_barrier();

    cout << "Done testing tstFlat_Mat_State_Builder on " << rtt_c4::node() 
	 << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstFlat_Mat_State_Builder.cc
//---------------------------------------------------------------------------//
