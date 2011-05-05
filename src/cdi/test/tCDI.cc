//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/tCDI.cc
 * \author Thomas M. Evans
 * \date   Tue Oct  9 15:52:01 2001
 * \brief  CDI test executable.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_test.hh"
#include "DummyGrayOpacity.hh"
#include "DummyMultigroupOpacity.hh"
#include "DummyOdfmgOpacity.hh"
#include "DummyEoS.hh"
#include "ds++/Release.hh"
#include "../CDI.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <typeinfo>
#include <iomanip>

using namespace std;

using rtt_cdi_test::match;
using rtt_cdi_test::DummyGrayOpacity;
using rtt_cdi_test::DummyMultigroupOpacity;
using rtt_cdi_test::DummyOdfmgOpacity;
using rtt_cdi_test::DummyEoS;
using rtt_cdi::CDI;
using rtt_cdi::GrayOpacity;
using rtt_cdi::MultigroupOpacity;
using rtt_cdi::OdfmgOpacity;
using rtt_cdi::EoS;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void check_CDI(const CDI &cdi)
{
    // check cdi, note that the different combinations of rtt_cdi::Model and
    // rtt_cdi::Reaction will yield the same results because DummyOpacity,
    // DummyMultigroupOpacity, and DummyEoS all yield the same stuff.  These
    // have all been tested in tDummyOpacity and tDummyEoS.  Here we just
    // check the types

    // check for gray
    if (typeid(*cdi.gray(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION)) == 
        typeid(DummyGrayOpacity))
    {
        PASSMSG("CDI gray() returned the correct type!");
    }
    else
    {
        FAILMSG("CDI gray() did not return the correct type!");
    }

    if (typeid(*cdi.gray(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)) == 
        typeid(DummyGrayOpacity))
    {
        PASSMSG("CDI gray() returned the correct type!");
    }
    else
    {
        FAILMSG("CDI gray() did not return the correct type!");
    }

    // check for multigroup
    if (typeid(*cdi.mg(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION)) == 
        typeid(DummyMultigroupOpacity))
    {
        PASSMSG("CDI mg() returned the correct type!");
    }
    else
    {
        FAILMSG("CDI mg() did not return the correct type!");
    }


    if (typeid(*cdi.mg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)) == 
        typeid(DummyMultigroupOpacity))
    {
        PASSMSG("CDI mg() returned the correct type!");
    }
    else
    {
        FAILMSG("CDI mg() did not return the correct type!");
    }

    //          check for odfmg
    if (typeid(*cdi.odfmg(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION)) == 
        typeid(DummyOdfmgOpacity))
    {
        PASSMSG("CDI odfmg() returned the correct type!");
    }
    else
    {
        FAILMSG("CDI odfmg() did not return the correct type!");
    }

    if (typeid(*cdi.odfmg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)) == 
        typeid(DummyOdfmgOpacity))
    {
        PASSMSG("CDI odfmg() returned the correct type!");
    }
    else
    {
        FAILMSG("CDI odfmg() did not return the correct type!");
    }


    if (typeid(*cdi.eos()) == typeid(DummyEoS))
    {
        PASSMSG("CDI eos() returned the correct type!");
    }
    else
    {
        FAILMSG("CDI eos() did not return the correct type!");
    }

    // gray test case: Find the value of opacity at T=0.35 keV and rho = 27.2
    // g/cm^3.  For DummyGrayOpacity the value should be .35272 cm^2/g.

    double temp       = 0.35;               // keV
    double dens       = 27.2;               // g/cm^3
    double refOpacity = temp + dens/1000.0; // cm^2/g

    double opacity = cdi.gray(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION)->
                     getOpacity(temp, dens);

    if (match(opacity, refOpacity))
    {
        PASSMSG("CDI.gray()->getOpacity is ok.");
    }
    else
    {
        FAILMSG("CDI.gray()->getOpacity is not ok.");
    }

    // mg test case: Find the mg opacities at T=0.35 keV and rho = 27.2
    // g/cm^3.  For DummyMultigroupOpacity the values should be { }.  Three
    // groups    
    size_t ng = 3;

    // The energy groups in DummyMultigroupOpacity are hardwired to
    // be { 0.05, 0.5, 5.0, 50.0 } keV.
    std::vector< double > energyBoundary(ng+1);
    energyBoundary[0] = 0.05;
    energyBoundary[1] = 0.5;
    energyBoundary[2] = 5.0;
    energyBoundary[3] = 50.0;

    if (cdi.mg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)->getNumGroups() == ng)
    {
        PASSMSG("CDI.mg()->getNumGroups() is ok.");
    }
    else
    {
        FAILMSG("CDI.mg()->getNumGroups() is not ok.");
    }

    std::vector< double > vRefOpacity( 3 );
    for ( size_t group=0; group<ng; ++group )
        vRefOpacity[group] = 2.0*(temp+dens/1000.0)
                             /(energyBoundary[group]+energyBoundary[group+1]);

    std::vector< double > vOpacity( ng );
    vOpacity = cdi.mg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)->
               getOpacity( temp, dens );

    if ( match( vOpacity, vRefOpacity ) )
    {
        PASSMSG("CDI.mg()->getOpacity(T,rho) is ok.");
    }
    else
    {
        FAILMSG("CDI.mg()->getOpacity(T,rho) is not ok.");
    }
	
    // odfmg test case: Find the mg opacities at T=0.35 keV and rho = 27.2
    // g/cm^3.  For DummyOdfmgOpacity the values should be { }.  Three
    // groups, four bands
    size_t numBands( 4);
    size_t numGroups(3);


    // The energy groups in DummyOdfmgOpacity are the same as in multigroup
    // be { 0.05, 0.5, 5.0, 50.0 } keV.
    energyBoundary.resize(numGroups + 1);
    energyBoundary[0] = 0.05;
    energyBoundary[1] = 0.5;
    energyBoundary[2] = 5.0;
    energyBoundary[3] = 50.0;

    if (cdi.odfmg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)->getNumGroups()
        == numGroups )
    {
        PASSMSG("CDI.odfmg()->getNumGroups() is ok.");
    }
    else
    {
        FAILMSG("CDI.odfmg()->getNumGroups() is not ok.");
    }

    std::vector< double > bandBoundary(numBands + 1);
    bandBoundary[0] = 0.0;
    bandBoundary[1] = 0.125;
    bandBoundary[2] = 0.25;
    bandBoundary[3] = 0.5;
    bandBoundary[4] = 1.0;

    if (cdi.odfmg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)->getNumBands()
        == numBands )
    {
        PASSMSG("CDI.odfmg()->getNumBands() is ok.");
    }
    else
    {
        FAILMSG("CDI.odfmg()->getNumBands() is not ok.");
    }

    std::vector< std::vector< double > > odfmgRefOpacity( numGroups );
    for ( size_t group=0; group<numGroups; ++group )
    {
        odfmgRefOpacity[group].resize( numBands );
        for ( size_t band=0; band<numBands; ++band )
        {
            odfmgRefOpacity[group][band]
                = 2.0 * (temp + dens/1000.0)
                / (energyBoundary[group] + energyBoundary[group+1])
                * pow(10.0, static_cast<int>(band) - 2);
        }
    }

    std::vector< std::vector< double > > odfmgOpacity( numGroups );

    odfmgOpacity = cdi.odfmg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)->
                   getOpacity( temp, dens );

    if ( match( odfmgOpacity, odfmgRefOpacity ) )
    {
        PASSMSG("CDI.odfmg()->getOpacity(T,rho) is ok.");
    }
    else
    {
        FAILMSG("CDI.odfmg()->getOpacity(T,rho) is not ok.");
    }


    // Test the EoS plug-in component of this CDI.

    double refCve = temp + dens/1000.0;

    double Cve = cdi.eos()->getElectronHeatCapacity( temp, dens );

    if ( match( Cve, refCve ) )
    {
        ostringstream message;
        message << "CDI.eos()->getElectronHeatCapacity( dbl, dbl )\n\t"
                << "returned a value that matched the reference value.";
        PASSMSG(message.str());
    }
    else
    {
        ostringstream message;
        message << "CDI.eos()->getElectronHeatCapacity( dbl, dbl )\n\t"
                << "returned a value that was out of spec.";
        FAILMSG(message.str());
    }

    // Print the material ID string

    // Test Rosseland integration with MG opacities.
    {
        // integrate on the interval [0.5, 5.0] keV for T=1.0 keV.
        double int_total1 = cdi.integrateRosselandSpectrum(0.5, 5.0, 1.0);
        // group 1 has the same energy range.
        double int_total2 = cdi.integrateRosselandSpectrum(2, 1.0);
        if (soft_equiv(int_total1, int_total2 , 1.0e-7))
        {
            ostringstream message;
            message.precision(10);
            message << "Calculated a total normalized Rosseland integral of "
                    << setw(12) << setiosflags(ios::fixed) << int_total2;
            PASSMSG(message.str());
        }
        else
        {
            ostringstream message;
            message.precision(10);
            message << "Calculated a total normalized Rosseland integral of "
                    << setw(12) << setiosflags(ios::fixed) << int_total2 
                    << " instead of " << setw(12) << setiosflags(ios::fixed)
                    << int_total1 << "." << endl;
            FAILMSG(message.str());
        }
    }

    return;
}

//---------------------------------------------------------------------------//

void test_CDI()
{
    // make SPs to opacity and EoS objects
    SP<const GrayOpacity>       gray_planck_abs;
    SP<const GrayOpacity>       gray_iso_scatter;
    SP<const MultigroupOpacity> mg_planck_abs;
    SP<const MultigroupOpacity> mg_iso_scatter;
    SP<const MultigroupOpacity> mg_diff_bound;
    SP<const OdfmgOpacity>      odfmg_planck_abs;
    SP<const OdfmgOpacity>      odfmg_iso_scatter;
    SP<const OdfmgOpacity>      odfmg_diff_bound;
    SP<const OdfmgOpacity>      odfmg_diff_bound2;
    SP<const EoS>               eos;

    // assign to dummy state objects
    gray_planck_abs  = new DummyGrayOpacity(rtt_cdi::ABSORPTION,
                                            rtt_cdi::PLANCK);
    gray_iso_scatter = new DummyGrayOpacity(rtt_cdi::SCATTERING,
                                            rtt_cdi::ISOTROPIC);

    mg_planck_abs    = new DummyMultigroupOpacity(rtt_cdi::ABSORPTION,
                                                  rtt_cdi::PLANCK);
    mg_iso_scatter   = new DummyMultigroupOpacity(rtt_cdi::SCATTERING,
                                                  rtt_cdi::ISOTROPIC);
    //multigroup with different boundaries
    mg_diff_bound    = new DummyMultigroupOpacity(rtt_cdi::SCATTERING,
                                                  rtt_cdi::THOMSON,
                                                  6);

    odfmg_planck_abs    = new DummyOdfmgOpacity(rtt_cdi::ABSORPTION,
                                                rtt_cdi::PLANCK);
    odfmg_iso_scatter   = new DummyOdfmgOpacity(rtt_cdi::SCATTERING,
                                                rtt_cdi::ISOTROPIC);
    //odfmg with different group boundaries
    odfmg_diff_bound    = new DummyOdfmgOpacity(rtt_cdi::SCATTERING,
                                                rtt_cdi::THOMSON,
                                                6, 5);
    //odfmg with different band boundaries
    odfmg_diff_bound2   = new DummyOdfmgOpacity(rtt_cdi::SCATTERING,
                                                rtt_cdi::THOMSON,
                                                4, 3);


    eos              = new DummyEoS();

    // make a CDI, it should be empty
    std::string const matName("dummyMaterial");
    CDI cdi( matName );
    for (size_t i = 0; i < rtt_cdi::constants::num_Models; i++)
        for (size_t j = 0; j < rtt_cdi::constants::num_Reactions; j++)
        {
            // casting these is inherently dangerous, but because this is a
            // test I won't go nuts
            rtt_cdi::Model m    = static_cast<rtt_cdi::Model>(i);
            rtt_cdi::Reaction r = static_cast<rtt_cdi::Reaction>(j);

            if (cdi.isGrayOpacitySet(m,r))       ITFAILS;
            if (cdi.isMultigroupOpacitySet(m,r)) ITFAILS;
            if (cdi.isOdfmgOpacitySet(m,r)) ITFAILS;
        }
    if (cdi.isEoSSet()) ITFAILS;

    if( cdi.getMatID() == matName )
    {
        PASSMSG("Good, the material identifier was set and fetched correctly.");
    }
    else
    {
        FAILMSG("Oh-ho, the material identifier was not set and fetched correctly.");
    }

    // there should be no energy group boundaries set yet
    if (CDI::getFrequencyGroupBoundaries().empty())
    {
        PASSMSG("Good, no frequency group boundaries defined yet.");
    }
    else
    {
        FAILMSG("Oh-oh, frequency boundaries are defined.");
    }

    // there should be no opacity band boundaries set yet
    if (CDI::getOpacityCdfBandBoundaries().empty())
    {
        PASSMSG("Good, no opacity band boundaries defined yet.");
    }
    else
    {
        FAILMSG("Uh-oh, opacity band boundaries are defined.");
    }

    // now assign stuff to it
    cdi.setGrayOpacity(gray_planck_abs);
    cdi.setGrayOpacity(gray_iso_scatter);
    cdi.setMultigroupOpacity(mg_planck_abs);
    cdi.setMultigroupOpacity(mg_iso_scatter);
    cdi.setOdfmgOpacity(odfmg_planck_abs);
    cdi.setOdfmgOpacity(odfmg_iso_scatter);
    cdi.setEoS(eos);

    // check the energy group boundaries
    {
        vector<double> b1 = CDI::getFrequencyGroupBoundaries();
        vector<double> b2 = cdi.mg(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION)->
                            getGroupBoundaries();
        vector<double> b3 = cdi.mg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)->
                            getGroupBoundaries();

        // these should all be equal
        bool test = true;

        if (!soft_equiv(b1.begin(), b1.end(), b2.begin(), b2.end()))
            test = false;

        if (!soft_equiv(b1.begin(), b1.end(), b3.begin(), b3.end()))
            test = false;

        if (test)
        {
            PASSMSG("All multigroup data has consistent energy groups.");
        }
        else
        {
            FAILMSG("Multigroup data has inconsistent energy groups.");
        }
    }

    // check the odfmg energy group boundaries
    {
        vector<double> b1 = CDI::getFrequencyGroupBoundaries();
        vector<double> b2 = cdi.odfmg(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION)->
                            getGroupBoundaries();
        vector<double> b3 = cdi.odfmg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)->
                            getGroupBoundaries();

        // these should all be equal
        bool test = true;

        if (!soft_equiv(b1.begin(), b1.end(), b2.begin(), b2.end()))
            test = false;

        if (!soft_equiv(b1.begin(), b1.end(), b3.begin(), b3.end()))
            test = false;

        if (test)
        {
            PASSMSG("All odfmg data has consistent energy groups.");
        }
        else
        {
            FAILMSG("Odfmg data has inconsistent energy groups.");
        }
    }
	
    // check the odfmg opacity band boundaries
    {
        vector<double> b1 = CDI::getOpacityCdfBandBoundaries();
        vector<double> b2 = cdi.odfmg(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION)->
                            getBandBoundaries();
        vector<double> b3 = cdi.odfmg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)->
                            getBandBoundaries();

        // these should all be equal
        bool test = true;

        if (!soft_equiv(b1.begin(), b1.end(), b2.begin(), b2.end()))
            test = false;

        if (!soft_equiv(b1.begin(), b1.end(), b3.begin(), b3.end()))
            test = false;

        if (test)
        {
            PASSMSG("All odfmg data has consistent opacity bands.");
        }
        else
        {
            FAILMSG("Odfmg data has inconsistent opacity bands.");
        }
    }
	
    // catch an exception when we try to assign a multigroup opacity to CDI
    // that has a different frequency group structure
    bool caught = false;
    try 
    {
        cdi.setMultigroupOpacity(mg_diff_bound);
    }
    catch (const rtt_dsxx::assertion &ass)
    {
        ostringstream message;
        message << "Good, we caught the following exception: \n"
                << ass.what();
        PASSMSG(message.str());
        caught = true;
    }
    if (!caught) 
    {
        ostringstream message;
        message << "Failed to catch an exception for setting a different "
                << "frequency group structure.";
        FAILMSG(message.str());
    } 

    // catch an exception when we try to assign a odfmg opacity to CDI
    // that has a different frequency group structure with odfmg
    caught = false;
    try 
    {
        cdi.setOdfmgOpacity(odfmg_diff_bound);
    }
    catch (const rtt_dsxx::assertion &ass)
    {
        ostringstream message;
        message << "Good, we caught the following exception when trying to assign bad odfmg frequency structure: \n"
                << ass.what();
        PASSMSG(message.str());
        caught = true;
    }
    if (!caught) 
    {
        ostringstream message;
        message << "Failed to catch an exception for setting a different "
                << "frequency group structure with odfmg.";
        FAILMSG(message.str());
    } 
	
    // catch an exception when we try to assign a odfmg opacity to CDI
    // that has a different opacity band structure with odfmg
    caught = false;
    try 
    {
        cdi.setOdfmgOpacity(odfmg_diff_bound2);
    }
    catch (const rtt_dsxx::assertion &ass)
    {
        ostringstream message;
        message << "Good, we caught the following exception when trying to assign bad odfmg band structure: \n"
                << ass.what();
        PASSMSG(message.str());
        caught = true;
    }
    if (!caught) 
    {
        ostringstream message;
        message << "Failed to catch an exception for setting a different "
                << "opacity band structure.";
        FAILMSG(message.str());
    } 


    // make sure these are assigned
    if (cdi.isGrayOpacitySet(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION))
    {
        PASSMSG("Gray planck absorption set!");
    }
    else
    {
        FAILMSG("Gray planck absorption not set!");
    }

    if (cdi.isGrayOpacitySet(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING))
    {
        PASSMSG("Gray isotropic scattering set!");
    }
    else
    {
        FAILMSG("Gray isotropic scattering not set!");
    }

    if (cdi.isMultigroupOpacitySet(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION))
    {
        PASSMSG("Multigroup planck (in-group) absorption set!");
    }
    else
    {
        FAILMSG("Multigroup planck (in-group) absorption not set!");
    }

    if (cdi.isMultigroupOpacitySet(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING))
    {
        PASSMSG("Multigroup isotropic scattering set!");
    }
    else
    {
        FAILMSG("Multigroup isotropic scattering not set!");
    }

    if (cdi.isOdfmgOpacitySet(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION))
    {
        PASSMSG("Odfmg planck (in-group) absorption set!");
    }
    else
    {
        FAILMSG("Odfmg planck (in-group) absorption not set!");
    }

    if (cdi.isOdfmgOpacitySet(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING))
    {
        PASSMSG("Odfmg isotropic scattering set!");
    }
    else
    {
        FAILMSG("Odfmg isotropic scattering not set!");
    }


    if (cdi.isEoSSet())
    {
        PASSMSG("EoS set!");
    }
    else
    {
        FAILMSG("EoS not set!");
    }

    // catch some exceptions
    caught = false;
    try 
    {
        cdi.setGrayOpacity(gray_planck_abs);
    }
    catch (const rtt_dsxx::assertion &ass)
    {
        ostringstream message;
        message << "Good, we caught the following exception: \n"
                << ass.what();
        PASSMSG(message.str());
        caught = true;
    }
    if (!caught) 
    {
        FAILMSG("Failed to catch overset exception!");
    }

    caught = false;
    try
    {
        cdi.mg(rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION);
    }
    catch (const rtt_dsxx::assertion &ass)
    {
        ostringstream message;
        message << "Good, we caught the following exception: \n"
                << ass.what();
        PASSMSG(message.str());
        caught = true;
    }
    if (!caught)
    {
        FAILMSG("Failed to catch an illegal access exception!");
    }

    // check the cdi through a function call
    check_CDI(cdi);

    // reset and make sure we are empty
    cdi.reset();
    if (!cdi.isGrayOpacitySet(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION))
    {
        PASSMSG("Gray planck absorption unset!");
    }
    else
    {
        FAILMSG("Gray planck absorption is still set!");
    }

    if (!cdi.isGrayOpacitySet(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING))
    {
        PASSMSG("Gray isotropic scattering unset!");
    }
    else
    {
        FAILMSG("Gray isotropic scattering is still set!");
    }

    if (!cdi.isMultigroupOpacitySet(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION))
    {
        PASSMSG("Multigroup planck (in-group) absorption unset!");
    }
    else
    {
        FAILMSG("Multigroup planck (in-group) absorption is still set!");
    }

    if (!cdi.isMultigroupOpacitySet(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING))
    {
        PASSMSG("Multigroup isotropic scattering unset!");
    }
    else
    {
        FAILMSG("Multigroup isotropic scattering is still set!");
    }

    if (!cdi.isOdfmgOpacitySet(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION))
    {
        PASSMSG("Odfmg planck (in-group) absorption unset!");
    }
    else
    {
        FAILMSG("Odfmg planck (in-group) absorption is still set!");
    }

    if (!cdi.isOdfmgOpacitySet(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING))
    {
        PASSMSG("Odfmg isotropic scattering unset!");
    }
    else
    {
        FAILMSG("Odfmg isotropic scattering is still set!");
    }


    if (!cdi.isEoSSet())
    {
        PASSMSG("EoS unset!");
    }
    else
    {
        FAILMSG("EoS is still set!");
    }

    caught = false;
    try
    {
        cdi.mg(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION);
    }
    catch (const rtt_dsxx::assertion &ass)
    {
        ostringstream message;
        message << "Good, we caught the following exception: \n"
                << ass.what();
        PASSMSG(message.str());
        caught = true;
    }
    if (!caught)
    {
        FAILMSG("Failed to catch an illegal access exception!");
    }

    // there should be no energy group boundaries set yet
    if (CDI::getFrequencyGroupBoundaries().empty())
    {
        PASSMSG("Good, no frequency group boundaries defined after reset.");
    }
    else
    {
        FAILMSG("Oh-oh, frequency boundaries are defined after reset.");
    }

    if (rtt_cdi_test::passed)
    {
        PASSMSG("Fundamental CDI Operations are ok.");
        cout << endl;
    }
}

//---------------------------------------------------------------------------//

void test_planck_integration()
{
    // We have not defined any group structure yet; thus, the Insist will
    // always fire if an integration is requested over a non-existent group 
    bool caught = false;
    try
    {
        CDI::integratePlanckSpectrum( 1, 1.0 );
    }
    catch(const rtt_dsxx::assertion &ass)
    {
        ostringstream message;
        message << "Caught illegal Planck calculation exception: \n"
                << "\t" << ass.what();
        PASSMSG(message.str());
        caught = true;
    }
    if (!caught)
    {
        FAILMSG("Did not catch an exception for calculating Planck integral.");
    }


    // check some planck integrals
    double int_total = CDI::integratePlanckSpectrum(0.0, 100.0, 1.0);
    if (soft_equiv(int_total, 1.0, 3.5e-10))
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a total normalized Planck integral of "
                << setw(12) << setiosflags(ios::fixed) << int_total;
        PASSMSG(message.str());
    }
    else
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a total normalized Planck integral of "
                << setw(12) << setiosflags(ios::fixed) << int_total 
                << " instead of 1.0.";
        FAILMSG(message.str());
    }

    double int_1 = CDI::integratePlanckSpectrum(0.0, 5.0, 10.0);
    if (soft_equiv(int_1, .00529316, 1.0e-6))
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a normalized Planck integral of "
                << setw(12) << setiosflags(ios::fixed) << int_1;
        PASSMSG(message.str());
    }
    else
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a normalized Planck integral of "
                << setw(12) << setiosflags(ios::fixed) << int_1 
                << " instead of .00529316.";
        FAILMSG(message.str());
    }

    double int_2 = CDI::integratePlanckSpectrum(0.0, .50, 10.0);
    if (soft_equiv(int_2, 6.29674e-6, 1.0e-6))
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a normalized Planck integral of "
                << setw(12) << setiosflags(ios::fixed) << int_2;
        PASSMSG(message.str());
    }
    else
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a normalized Planck integral of "
                << setw(12) << setiosflags(ios::fixed) << int_2
                << " instead of 6.29674e-6.";
        FAILMSG(message.str());
    }

    double int_0 = CDI::integratePlanckSpectrum(0.0, 0.0, 10.0);
    if (!soft_equiv(int_0, 0.0, 1.e-6)) ITFAILS;

    double int_range = CDI::integratePlanckSpectrum(0.1, 1.0, 1.0);
    if (!soft_equiv(int_range, 0.0345683, 3.e-5)) ITFAILS;

    // Catch an illegal group assertion.
    CDI cdi;
    SP<const MultigroupOpacity> mg(
        new DummyMultigroupOpacity(rtt_cdi::SCATTERING, rtt_cdi::THOMSON)); 
    cdi.setMultigroupOpacity(mg);

    // Check the normalized planck integrals.
    if (CDI::getNumberFrequencyGroups() != 3) ITFAILS;
    if (CDI::getNumberOpacityBands() != 0) ITFAILS;

    double g1_integral = CDI::integratePlanckSpectrum(1, 1.0);
    double g2_integral = CDI::integratePlanckSpectrum(2, 1.0);
    double g3_integral = CDI::integratePlanckSpectrum(3, 1.0);
    double g_total     = CDI::integratePlanckSpectrum(   1.0);

    if (soft_equiv(g1_integral, 0.00528686, 1.e-6))
    {
        PASSMSG("Group 1 integral within tolerance.");
    }
    else
    {
        FAILMSG("Group 1 integral fails tolerance.");
    }

    if (soft_equiv(g2_integral, 0.74924, 1.e-6))
    {
        PASSMSG("Group 2 integral within tolerance.");
    }
    else
    {
        FAILMSG("Group 2 integral fails tolerance.");
    }

    if (soft_equiv(g3_integral, 0.245467, 1.e-6))
    {
        PASSMSG("Group 3 integral within tolerance.");
    }
    else
    {
        FAILMSG("Group 3 integral fails tolerance.");
    }

    if (soft_equiv(g_total, 0.999994, 1.e-6))
    {
        PASSMSG("Total integral over groups within tolerance.");
    }
    else
    {
        FAILMSG("Total integral over groups fails tolerance.");
    }

    // Test that a zero temperature returns a zero.
    if (soft_equiv(CDI::integratePlanckSpectrum(0.0, 100.0, 0.0), 0.0))
    {
        PASSMSG("Planck integral from hnu=0 to 100 at T=0 is zero: good!");
    }
    else
    {
        FAILMSG("Planck integral from hnu=0 to 100 at T=0 is not zero: BAD!");
    }

    if (soft_equiv(CDI::integratePlanckSpectrum(1, 0.0), 0.0))
    {
        PASSMSG("Planck integral of group 1 at T=0 is zero: good!");
    }
    else
    {
        FAILMSG("Planck integral of group 1 at T=0 is not zero: BAD!");
    }

    if (soft_equiv(CDI::integratePlanckSpectrum(2, 0.0), 0.0))
    {
        PASSMSG("Planck integral of group 2 at T=0 is zero: good!");
    }
    else
    {
        FAILMSG("Planck integral of group 2 at T=0 is not zero: BAD!");
    }

    if (soft_equiv(CDI::integratePlanckSpectrum(3, 0.0), 0.0))
    {
        PASSMSG("Planck integral of group 3 at T=0 is zero: good!");
    }
    else
    {
        FAILMSG("Planck integral of group 3 at T=0 is not zero: BAD!");
    }

    if (soft_equiv(CDI::integratePlanckSpectrum(0.0), 0.0))
    {
        PASSMSG("Planck integral over all groups at T=0 is zero: good!");
    }
    else
    {
        FAILMSG("Planck integral over all groups at T=0 is not zero: BAD!");
    }

    if (rtt_cdi_test::passed)
    {
        PASSMSG("All Planckian integral tests ok.");
        cout << endl;
    }


    // Compare the integration over all groups to the individual groups:

    std::vector<double> group_bounds ( CDI::getFrequencyGroupBoundaries() );
    if (group_bounds.size() != 4) ITFAILS;

    std::vector<double> planck;

    CDI::integrate_Planckian_Spectrum(group_bounds, 1.0, planck);

    for (int group_index = 1; group_index <= 3; ++group_index)
    {

        double planck_g = CDI::integratePlanckSpectrum(group_index, 1.0);

        if (!soft_equiv( planck[group_index-1], planck_g) ) ITFAILS;

    }

    // Test zero temperature special case
    CDI::integrate_Planckian_Spectrum(group_bounds, 0.0, planck);

    for (int group_index = 1; group_index <= 3; ++group_index)
    {

        double planck_g = CDI::integratePlanckSpectrum(group_index, 1.0);

        if (!soft_equiv( planck[group_index-1], planck_g) ) ITFAILS;

    }

    if (rtt_cdi_test::passed)
    {
        PASSMSG("Group-wise and Full spectrum Planckian and Rosseland integrals match.");
    }
    else
    {
        FAILMSG("Group-wise and Full spectrum Planckian and Rosseland integrals do not match.");
    }



}

//---------------------------------------------------------------------------//

void test_rosseland_integration()
{

    // Only report this as a failure if 1) the error was not caught AND 2)
    // the Require macro is available.
    bool dbc_require( DBC & 1 );
    if( dbc_require )
    {
        bool caught = false;
        try
        {
            CDI::integrateRosselandSpectrum(0, 1.0);
        }
        catch(const rtt_dsxx::assertion &ass)
        {
            ostringstream message;
            message << "Caught illegal Rosseland calculation exception: \n";
            PASSMSG(message.str());
            caught = true;
        }
        if (!caught)
        {
            FAILMSG("Did not catch an exception for calculating Rosseland integral.");
        }
        // catch our assertion
        caught = false;
        double P,R;
        try
        {
            CDI::integrate_Rosseland_Planckian_Spectrum(0, 1.0, P,R);
        }
        catch(const rtt_dsxx::assertion &ass)
        {
            ostringstream message;
            message << "Caught illegal Rosseland and Planckian "
                    << "calculation exception: \n";
            PASSMSG(message.str());
            caught = true;
        }
        if (!caught)
        {
            FAILMSG("Did not catch an exception for calculating Rosseland and Planckian integral.");
        }
    }

    // check some planck integrals
    double int_total = CDI::integrateRosselandSpectrum(0, 100.0, 1.0);
    if (soft_equiv(int_total, 1.0, 1.0e-7))
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a total normalized Rosseland integral of "
                << setw(12) << setiosflags(ios::fixed) << int_total;
        PASSMSG(message.str());
    }
    else
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a total normalized Rosseland integral of "
                << setw(12) << setiosflags(ios::fixed) << int_total 
                << " instead of 1.0.";
        FAILMSG(message.str());
    }

    double int_1 = CDI::integratePlanckSpectrum(0.1, 1.0, 1.0);
    if (soft_equiv(int_1, .0345683, 1.0e-5))
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a normalized Planck integral of "
                << setw(12) << setiosflags(ios::fixed) << int_1;
        PASSMSG(message.str());
    }
    else
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a normalized Planck integral of "
                << setw(12) << setiosflags(ios::fixed) << int_1 
                << " instead of .0345683";
        FAILMSG(message.str());
    }

    double int_2 = CDI::integrateRosselandSpectrum(0.1, 1.0, 1.0);
    if (soft_equiv(int_2, 0.01220025, 1.0e-5))
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a normalized Rosseland integral of "
                << setw(12) << setiosflags(ios::fixed) << int_2;
        PASSMSG(message.str());
    }
    else
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a normalized Rosseland integral of "
                << setw(12) << setiosflags(ios::fixed) << int_2
                << " instead of 0.01220025";
        FAILMSG(message.str());
    }

    double PL,ROS;
    CDI::integrate_Rosseland_Planckian_Spectrum(0.1, 1.0, 1.0, PL, ROS);
    if (soft_equiv(PL, .0345683, 1.e-5))
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a  normalized Planck integral for RosselandSpectrum"
                << setw(12) << setiosflags(ios::fixed) << int_2;
        PASSMSG(message.str());
    }
    else
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a  normalized Planck integral for RosselandSpectrum"
                << setw(12) << setiosflags(ios::fixed) << int_2
                << " instead of .0345683.";
        FAILMSG(message.str());
    }

    if (soft_equiv(ROS, 0.01220025, 1.e-5))
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a  normalized Rosseland integral for RosselandSpectrum"
                << setw(12) << setiosflags(ios::fixed) << int_2;
        PASSMSG(message.str());
    }
    else
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a  normalized Rosseland integral for RosselandSpectrum"
                << setw(12) << setiosflags(ios::fixed) << int_2
                << " instead of 0.01220025.";
        FAILMSG(message.str());
    }

    CDI::integrate_Rosseland_Planckian_Spectrum(0.1, 1.0, 1.0,PL,ROS);
    if (soft_equiv(PL, .0345683, 1.e-5))
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a  normalized Planck integral for RosselandPlanckianSpectrum"
                << setw(12) << setiosflags(ios::fixed) << int_2;
        PASSMSG(message.str());
    }
    else
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a normalized Planck integral for RosselandPlanckianSpectrum"
                << setw(12) << setiosflags(ios::fixed) << int_2
                << " instead of .0345683.";
        FAILMSG(message.str());
    }

    if (soft_equiv(ROS, 0.01220025, 1.e-5))
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a  normalized Rosseland integral for RosselandPlanckianSpectrum"
                << setw(12) << setiosflags(ios::fixed) << int_2;
        PASSMSG(message.str());
    }
    else
    {
        ostringstream message;
        message.precision(10);
        message << "Calculated a normalized Rosseland integral for RosselandPlanckianSpectrum"
                << setw(12) << setiosflags(ios::fixed) << int_2
                << " instead of 0.01220025.";
        FAILMSG(message.str());
    }

    // Test that a zero temperature returns a zero
    if (soft_equiv(CDI::integrateRosselandSpectrum(0.0, 100.0, 0.0), 0.0))
    {
        PASSMSG("Rosseland integral from hnu=0 to 100 at T=0 is zero: good!");
    }
    else
    {
        FAILMSG("Rosseland integral from hnu=0 to 100 at T=0 is not zero: BAD!");
    }
    CDI::integrate_Rosseland_Planckian_Spectrum(0.0, 100.0, 0.0, PL, ROS);
    if (soft_equiv(PL, 0.0))
    {
        PASSMSG("Rosseland call for Planck integral  at T=0 is zero: good!");
    }
    else
    {
        FAILMSG("Rosseland call for Planck integral at T=0 is not zero: BAD!");
    }
    if (soft_equiv(ROS, 0.0))
    {
        PASSMSG("Rosseland integral  at T=0 is zero: good!");
    }
    else
    {
        FAILMSG("Rosseland integral  at T=0 is not zero: BAD!");
    }

    // check the normalized planck integrals
    if (CDI::getNumberFrequencyGroups() != 3) ITFAILS;

    // First group
    CDI::integrate_Rosseland_Planckian_Spectrum(1, 1.0, PL, ROS);
    if (!soft_equiv(PL,  0.005286862763740451, 1.e-6)) ITFAILS;
    if (!soft_equiv(ROS, 0.00158258277444842, 1.e-5)) ITFAILS;

    if (rtt_cdi_test::passed)
    {
        PASSMSG("Group 1 Rosseland and Planck integrals ok.");
    }
    else
    {
        FAILMSG("Group 1 Rosseland and Planck integrals failed.");
    }

    // Second group
    CDI::integrate_Rosseland_Planckian_Spectrum(2, 1.0, PL, ROS);
    if (!soft_equiv(PL,  0.7492399297, 1.e-6)) ITFAILS;
    if (!soft_equiv(ROS, 0.5897280880, 1.e-6)) ITFAILS;

    if (rtt_cdi_test::passed)
    {
        PASSMSG("Group 2 Rosseland and Planck integrals ok.");
    }
    else
    {
        FAILMSG("Group 2 Rosseland and Planck integrals failed.");
    }



    // Third group
    CDI::integrate_Rosseland_Planckian_Spectrum(3, 1.0, PL, ROS);
    if (!soft_equiv(PL,  0.2454669108, 1.e-6)) ITFAILS;
    if (!soft_equiv(ROS, 0.4086877254, 1.e-6)) ITFAILS;

    if (rtt_cdi_test::passed)
    {
        PASSMSG("Group 3 Rosseland and Planck integrals ok.");
    }
    else
    {
        FAILMSG("Group 3 Rosseland and Planck integrals failed.");
    }



    // All groups:
    std::vector<double> group_bounds ( CDI::getFrequencyGroupBoundaries() );
    if (group_bounds.size() != 4) ITFAILS;

    std::vector<double> planck;
    std::vector<double> rosseland;

    CDI::integrate_Rosseland_Planckian_Spectrum(group_bounds, 1.0, planck, rosseland);

    for (int group_index = 1; group_index <= 3; ++group_index)
    {

        CDI::integrate_Rosseland_Planckian_Spectrum(group_index, 1.0, PL, ROS);

        if (!soft_equiv(planck   [group_index-1], PL )) ITFAILS;
        if (!soft_equiv(rosseland[group_index-1], ROS)) ITFAILS;

    }

    // Special case of zero temperature
    CDI::integrate_Rosseland_Planckian_Spectrum(3, 0.0, PL, ROS);
    if (!soft_equiv(PL,  0.0, 1.e-6)) ITFAILS;
    if (!soft_equiv(ROS, 0.0, 1.e-6)) ITFAILS;

    if (rtt_cdi_test::passed)
    {
        PASSMSG("Zero T Rosseland and Planck integrals ok.");
    }
    else
    {
        FAILMSG("Zero T Rosseland and Planck integrals failed.");
    }
        
    if (rtt_cdi_test::passed)
    {
        ostringstream msg;
        msg << "Group-wize and Full spectrum Planckian and Rosseland "
            << "integrals match.";
        PASSMSG( msg.str().c_str() );
    }
    else
    {
        ostringstream msg;
        msg << "Group-wize and Full spectrum Planckian and Rosseland "
            << "integrals do not match.";
        FAILMSG( msg.str().c_str() );
    }

    if (rtt_cdi_test::passed)
    {
        PASSMSG("All Rosseland and Rosseland/Planckian integral tests ok.");
        cout << endl;
    }
    return;
}

//---------------------------------------------------------------------------//

void printPkgVer()
{
    std::cout << "This is Draco package CDI.\n"
              << "Version: " <<  rtt_dsxx::release()
              << std::endl << std::endl;
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
        if (string(argv[arg]) == "--version")
        {
            cout << argv[0] << ": version " << rtt_dsxx::release() 
                 << endl;
            return 0;
        }

    try
    {
        // >>> UNIT TESTS
        printPkgVer();

        test_CDI();

        test_planck_integration();

        test_rosseland_integration();

    }
    catch (rtt_dsxx::assertion &ass)
    {
        cout << "While testing tCDI, " << ass.what()
             << endl;
        return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_cdi_test::passed) 
    {
        cout << "**** tCDI Test: PASSED" 
             << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tCDI." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tCDI.cc
//---------------------------------------------------------------------------//
