//----------------------------------*-C++-*----------------------------------//
// AMR_Interface.hh
// Thomas M. Evans
// Thu Jul 16 09:25:43 1998
//---------------------------------------------------------------------------//
// @> Rage Interface functions and classes
//---------------------------------------------------------------------------//

#ifndef __imc_AMR_Interface_hh__
#define __imc_AMR_Interface_hh__

//===========================================================================//
// class AMR_Interface - 
//
// Purpose : Interface functions to Rage's AMR-Eulerian Code.
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "imc/Names.hh"

//===========================================================================//
// F90 Functional Interface to Rage 
//===========================================================================//
// direct functional interface to F90 Rage code

extern "C"
{
    void rage_IMC_(int i);
}

//===========================================================================//
// class AMR_Interface
//===========================================================================//

IMCSPACE

class AMR_Interface
{
};

CSPACE

#endif                          // __imc_AMR_Interface_hh__

//---------------------------------------------------------------------------//
//                              end of imc/AMR_Interface.hh
//---------------------------------------------------------------------------//
