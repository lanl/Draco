/*---------------------------------------------*-C-*----------------------------------------------*/
/*!
 * \file   c4_mpi.h
 * \author Thomas M. Evans
 * \brief  put the right includes for MPI header files
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
/*------------------------------------------------------------------------------------------------*/

#ifndef rtt_c4_c4_mpi_h
#define rtt_c4_c4_mpi_h

#include "c4/config.h"

/* defined in config/setupMPI.cmake */
#ifdef C4_MPI
#include <mpi.h>
#endif

#ifndef DRACO_MAX_PROCESSOR_NAME
#ifdef MPI_MAX_PROCESSOR_NAME
#define DRACO_MAX_PROCESSOR_NAME MPI_MAX_PROCESSOR_NAME
#endif
#endif

#endif /* __c4_c4_mpi_h__ */

/*------------------------------------------------------------------------------------------------*/
/* end of c4_mpi.h */
/*------------------------------------------------------------------------------------------------*/
