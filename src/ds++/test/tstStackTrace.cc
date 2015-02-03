//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstStackTrace.cc
 * \author Kelly Thompson
 * \date   Thursday, Dec 19, 2013, 15:09 pm
 * \brief  Demonstrate/Test fpe_trap's print_stacktrace function.
 * \note   Copyright (C) 2013-2015 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../StackTrace.hh"
#include "../ScalarUnitTest.hh"
#include "../Release.hh"
#include <sstream>

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void sr2(std::string & msg)
{
    msg = rtt_dsxx::print_stacktrace( "myMessage" );
    return;
}

void sr1(std::string & msg)
{
    sr2(msg);
    return;
}

void runtest( rtt_dsxx::UnitTest & ut )
{
    std::cout << "Running tstStackTrace...\n\n"
              << "Requesting a trace..." << std::endl;


    // Create a stack trace.  It should look something like this:

    // Stack trace:
    // Signaling error: myMessage
    // Process /var/tmp/kgt/gcc-mpid-diag3/d/src/ds++/test/tstStackTrace
    
    // (PID:2849) ./tstStackTrace : sr1(std::string&)()+0x18
    // (PID:2849) ./tstStackTrace : runtest(rtt_dsxx::UnitTest&)()+0x54
    // (PID:2849) ./tstStackTrace : main()+0x50
    // (PID:2849) /lib64/libc.so.6 : __libc_start_main()+0xfd
    // (PID:2849) ./tstStackTrace() [0x4054f9] : ??
    // Stack trace: END (PID:2849)
    
    std::string trace;
    sr1(trace);
    std::cout << trace << std::endl;

    // store the trace in an ostring stream.
    std::ostringstream msg;
    msg << trace << std::endl;
    
    // Check for word counts
    bool const verbose(false);
    std::map< std::string, unsigned > words = ut.get_word_count( msg,
                                                                 verbose );

	// Expected values.
    if( words[ std::string("PID")]   != 1 ) ITFAILS;
    if( words[ std::string("Stack")] <  2 ) ITFAILS;
   
    if( ut.numFails==0 )
        PASSMSG("done with testing print_stacktrace.");
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_dsxx::ScalarUnitTest ut( argc, argv, rtt_dsxx::release );
    try { runtest(ut); }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
//  end of tstSafe_Divide.cc
//---------------------------------------------------------------------------//
