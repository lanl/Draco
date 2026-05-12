/*---------------------------------------------*-C-*----------------------------------------------*/
/*!
 * \file   rng/test/time_serial.h
 * \author Gabriel M. Rockefeller
 * \brief  time_serial header file.
 * \note   Copyright (C) 2012-2026 Triad National Security, LLC., All rights reserved.
 */
/*------------------------------------------------------------------------------------------------*/

#ifndef rng_time_serial_h
#define rng_time_serial_h

/* This is a Random123 file. Ignore advanced gcc warnings emanating from this
 * file. */
#ifdef __GNUC__
#pragma GCC system_header
#endif

#ifdef _MSC_FULL_VER
// - 4204 :: nonstandard extension used: non-constant aggregate initializer.
// - 4127 :: conditional expression is constant
#pragma warning(push)
#pragma warning(disable : 4204 4100 4127)
#endif

#include <Random123/aes.h>
#include <Random123/ars.h>
#include <Random123/philox.h>
#include <Random123/threefry.h>

#ifdef _MSC_FULL_VER
#pragma warning(pop)
#endif

#endif /* rng_time_serial_h */

/*------------------------------------------------------------------------------------------------*/
/* end of rng/time_serial.h */
/*------------------------------------------------------------------------------------------------*/
