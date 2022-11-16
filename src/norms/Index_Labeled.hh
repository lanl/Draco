//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   norms/Index_Labeled.hh
 * \author Rob Lowrie
 * \date   Fri Jan 14 13:57:58 2005
 * \brief  Header for Index_Labeled.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_norms_Index_Labeled_hh
#define rtt_norms_Index_Labeled_hh

#include "ds++/config.h"
#include <string>

namespace rtt_norms {

//================================================================================================//
/*!
 * \struct Index_Labeled
 * \brief  A Norms index class that contains processor number and a label.
 */
//================================================================================================//

struct Index_Labeled {
  //! The index.
  size_t index;

  //! The processor number.
  size_t processor;

  //! A string label.
  std::string label;

  //! Constructor.  Allow auto-casting from index.
  explicit Index_Labeled(const size_t index_ = 0, std::string label_ = "");

  //! Equality operator.
  bool operator==(const Index_Labeled &rhs) const;
};

} // end namespace rtt_norms

#endif // rtt_norms_Index_Labeled_hh

//------------------------------------------------------------------------------------------------//
// end of norms/Index_Labeled.hh
//------------------------------------------------------------------------------------------------//
