//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   units/test/tstUnitSystemEnums.cc
 * \author Kelly Thompson
 * \date   Wed Oct  8 13:50:19 2003
 * \brief  
 * \note   Copyright � 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <string>
#include <sstream>

#include "ds++/Soft_Equivalence.hh"
#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "../FundUnit.hh"
#include "../UnitSystemEnums.hh"
#include "units_test.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_enumValues()
{
    using std::string;
    using rtt_dsxx::soft_equiv;

    {
	int    const iSenVal(3);
	double const adSenVal[iSenVal] = { 0.0, 1.0, 100.0 };
	string const sSenVal( "NA,m,cm" );
	string const sSenVal2( "no length unit specified,meter,centimeter");

	if( rtt_units::num_Ltype == iSenVal )
	{
	    PASSMSG("num_Ltype has the expected value.");
	}
	else
	{
	    FAILMSG("num_Ltype does not have the expected value.");
	}

	if( soft_equiv( rtt_units::L_cf, rtt_units::L_cf+iSenVal, 
			adSenVal,        adSenVal+iSenVal ) )
	{
	    PASSMSG("L_cf has the expected values.");
	}
	else
	{
	    FAILMSG("L_cf does not have the expected values.");
	}

	if( rtt_units::L_labels == sSenVal )
	{
	    PASSMSG("L_labels has the expected values.");
	}
	else
	{
	    PASSMSG("L_labels does not have the expected values.");
	}

	if( rtt_units::L_long_labels == sSenVal2 )
	{
	    PASSMSG("L_long_labels has the expected values.");
	}
	else
	{
	    PASSMSG("L_long_labels does not have the expected values.");
	}
    }

//---------------------------------------------------------------------------//

    {
	int    const iSenVal(3);
	double const adSenVal[iSenVal] = { 0.0, 1.0, 1000.0 };
	string const sSenVal( "NA,kg,g" );
	string const sSenVal2( "no mass unit specified,kilogram,gram");

	if( rtt_units::num_Mtype == iSenVal )
	{
	    PASSMSG("num_Mtype has the expected value.");
	}
	else
	{
	    FAILMSG("num_Mtype does not have the expected value.");
	}

	if( soft_equiv( rtt_units::M_cf, rtt_units::M_cf+iSenVal, 
			adSenVal,        adSenVal+iSenVal ) )
	{
	    PASSMSG("M_cf has the expected values.");
	}
	else
	{
	    FAILMSG("M_cf does not have the expected values.");
	}

	if( rtt_units::M_labels == sSenVal )
	{
	    PASSMSG("M_labels has the expected values.");
	}
	else
	{
	    PASSMSG("M_labels does not have the expected values.");
	}

	if( rtt_units::M_long_labels == sSenVal2 )
	{
	    PASSMSG("M_long_labels has the expected values.");
	}
	else
	{
	    PASSMSG("M_long_labels does not have the expected values.");
	}
    }
//---------------------------------------------------------------------------//
    {
	int    const iSenVal(6);
	double const adSenVal[iSenVal] = { 0.0, 1.0, 1000.0, 1.0e6, 1.0e8, 1.0e9 };
	string const sSenVal( "NA,s,ms,us,sh,ns" );
	string const sSenVal2( "no time unit specified,second,milisecond,microsecond,shake,nanosecond");

	if( rtt_units::num_ttype == iSenVal )
	{
	    PASSMSG("num_ttype has the expected value.");
	}
	else
	{
	    FAILMSG("num_ttype does not have the expected value.");
	}

	if( soft_equiv( rtt_units::t_cf, rtt_units::t_cf+iSenVal, 
			adSenVal,        adSenVal+iSenVal ) )
	{
	    PASSMSG("t_cf has the expected values.");
	}
	else
	{
	    FAILMSG("t_cf does not have the expected values.");
	}

	if( rtt_units::t_labels == sSenVal )
	{
	    PASSMSG("t_labels has the expected values.");
	}
	else
	{
	    PASSMSG("t_labels does not have the expected values.");
	}

	if( rtt_units::t_long_labels == sSenVal2 )
	{
	    PASSMSG("t_long_labels has the expected values.");
	}
	else
	{
	    PASSMSG("t_long_labels does not have the expected values.");
	}
    }
//---------------------------------------------------------------------------//
    {
	int    const iSenVal(3);
	double const adSenVal[iSenVal] = { 0.0, 1.0, 1.0/1.1604449278e7 };
	string const sSenVal( "NA,K,keV" );

	if( rtt_units::num_Ttype == iSenVal )
	{
	    PASSMSG("num_Ttype has the expected value.");
	}
	else
	{
	    FAILMSG("num_Ttype does not have the expected value.");
	}

	if( soft_equiv( rtt_units::T_cf, rtt_units::T_cf+iSenVal, 
			adSenVal,        adSenVal+iSenVal ) )
	{
	    PASSMSG("T_cf has the expected values.");
	}
	else
	{
	    FAILMSG("T_cf does not have the expected values.");
	}

	if( rtt_units::T_labels == sSenVal )
	{
	    PASSMSG("T_labels has the expected values.");
	}
	else
	{
	    PASSMSG("T_labels does not have the expected values.");
	}
    }
//---------------------------------------------------------------------------//
    {
	int    const iSenVal(2);
	double const adSenVal[iSenVal] = { 0.0, 1.0 };
	string const sSenVal( "NA,Amp" );

	if( rtt_units::num_Itype == iSenVal )
	{
	    PASSMSG("num_Itype has the expected value.");
	}
	else
	{
	    FAILMSG("num_Itype does not have the expected value.");
	}

	if( soft_equiv( rtt_units::I_cf, rtt_units::I_cf+iSenVal, 
			adSenVal,        adSenVal+iSenVal ) )
	{
	    PASSMSG("I_cf has the expected values.");
	}
	else
	{
	    FAILMSG("I_cf does not have the expected values.");
	}

	if( rtt_units::I_labels == sSenVal )
	{
	    PASSMSG("I_labels has the expected values.");
	}
	else
	{
	    PASSMSG("I_labels does not have the expected values.");
	}
    }
//---------------------------------------------------------------------------//
    {
	int    const iSenVal(3);
	double const adSenVal[iSenVal] = { 0.0, 1.0, 57.295779512896171 };
	string const sSenVal( "NA,rad,deg" );

	if( rtt_units::num_Atype == iSenVal )
	{
	    PASSMSG("num_Atype has the expected value.");
	}
	else
	{
	    FAILMSG("num_Atype does not have the expected value.");
	}

	if( soft_equiv( rtt_units::A_cf, rtt_units::A_cf+iSenVal, 
			adSenVal,        adSenVal+iSenVal ) )
	{
	    PASSMSG("A_cf has the expected values.");
	}
	else
	{
	    FAILMSG("A_cf does not have the expected values.");
	}

	if( rtt_units::A_labels == sSenVal )
	{
	    PASSMSG("A_labels has the expected values.");
	}
	else
	{
	    PASSMSG("A_labels does not have the expected values.");
	}
    }
//---------------------------------------------------------------------------//
    {
	int    const iSenVal(2);
	double const adSenVal[iSenVal] = { 0.0, 1.0 };
	string const sSenVal( "NA,mol" );

	if( rtt_units::num_Qtype == iSenVal )
	{
	    PASSMSG("num_Qtype has the expected value.");
	}
	else
	{
	    FAILMSG("num_Qtype does not have the expected value.");
	}

	if( soft_equiv( rtt_units::Q_cf, rtt_units::Q_cf+iSenVal, 
			adSenVal,        adSenVal+iSenVal ) )
	{
	    PASSMSG("Q_cf has the expected values.");
	}
	else
	{
	    FAILMSG("Q_cf does not have the expected values.");
	}

	if( rtt_units::Q_labels == sSenVal )
	{
	    PASSMSG("Q_labels has the expected values.");
	}
	else
	{
	    PASSMSG("Q_labels does not have the expected values.");
	}
    }

//---------------------------------------------------------------------------//



    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    using std::cout;
    using std::endl;
    using std::string;
    using rtt_dsxx::release;

    // version tag
    for( int arg = 1; arg < argc; arg++ )
    {
	if( string( argv[arg] ) == "--version" )
	{
	    cout << argv[0] << ": version " << release() << endl;
	    return 0;
	}
    }

    try  //  >>>>> UNIT TESTS <<<<<
    {
	test_enumValues();
    }
    catch( rtt_dsxx::assertion &assert )
    {
	cout << "While testing tstUnitSystemEnums, " << assert.what() << endl;
	return 1;
    }

    // status of test
    cout << endl << "*********************************************" << endl;
    if( rtt_units_test::passed )
        cout << "**** tstUnitSystemEnums Test: PASSED" << endl;
    cout <<         "*********************************************" << endl << endl;
    
    cout << "Done testing tstUnitSystemEnums." << endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstUnitSystemEnums.cc
//---------------------------------------------------------------------------//
