//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   norms/Index_Proc.hh
 * \author Rob Lowrie
 * \date   Fri Jan 14 13:57:58 2005
 * \brief  Header for Index_Proc.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_norms_Index_Proc_hh
#define rtt_norms_Index_Proc_hh

#include <cstdlib>

namespace rtt_norms
{

//===========================================================================//
/*!
 * \struct Index_Proc
 * \brief  A Norms index class that contains the processor number.
 */
//===========================================================================//

struct Index_Proc 
{
    //! The index.
    size_t index;

    //! The processor number.
    size_t processor;

    // Constructor.  Allow auto-casting from index.
    Index_Proc(const size_t index_ = 0);

    // Equality operator.
    bool operator==(const Index_Proc &rhs) const;
};

} // end namespace rtt_norms

#endif // rtt_norms_Index_Proc_hh

//---------------------------------------------------------------------------//
//              end of norms/Index_Proc.hh
//---------------------------------------------------------------------------//
