//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   lapack_wrap/Blas_Prototypes.hh
 * \brief  Header declaring BLAS prototypes
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_lapack_wrap_Blas_Prototypes_hh
#define rtt_lapack_wrap_Blas_Prototypes_hh

#include <lapack_wrap/config.h>

extern "C" {
// >>> LEVEL 1 BLAS

// y <- x
void FC_GLOBAL(scopy, SCOPY)(int *, float *, int *, float *, int *);
void FC_GLOBAL(dcopy, DCOPY)(int *, double *, int *, double *, int *);

// x <- ax
void FC_GLOBAL(sscal, SSCAL)(int *, float *, float *, int *);
void FC_GLOBAL(dscal, DSCAL)(int *, double *, double *, int *);

// dot <- x^T y
float FC_GLOBAL(sdot, SDOT)(int *, float *, int *, float *, int *);
double FC_GLOBAL(ddot, DDOT)(int *, double *, int *, double *, int *);

// y <- ax + y
void FC_GLOBAL(saxpy, SAXPY)(int *, float *, float *, int *, float *, int *);
void FC_GLOBAL(daxpy, DAXPY)(int *, double *, double *, int *, double *, int *);

// nrm2 <- ||x||_2
float FC_GLOBAL(snrm2, SNRM2)(int *, float *, int *);
double FC_GLOBAL(dnrm2, DNRM2)(int *, double *, int *);

} // end of extern "C"

#endif // rtt_lapack_wrap_Blas_Prototypes_hh

//------------------------------------------------------------------------------------------------//
// end of lapack_wrap/Blas_Prototypes.hh
//------------------------------------------------------------------------------------------------//
