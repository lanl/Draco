//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/SP.hh
 * \author Geoffrey Furnish, Thomas Evans
 * \date   Tue Feb  4 11:27:53 2003
 * \brief  Smart Point class file.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef RTT_ds_SP_HH
#define RTT_ds_SP_HH

#include <memory>

namespace rtt_dsxx {

template <class T> using SP = std::shared_ptr<T>;

} // end namespace rtt_dsxx

#endif // RTT_ds_SP_HH

//------------------------------------------------------------------------------------------------//
// end of ds++/SP.hh
//------------------------------------------------------------------------------------------------//
