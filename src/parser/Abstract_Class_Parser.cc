//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/Abstract_Class_Parser.cc
 * \author Kent Budge
 * \brief  Define destructor for Abstract_Class_Parser_Base
 * \note   Copyright (C) 2012 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Abstract_Class_Parser.hh"

namespace rtt_parser
{
Abstract_Class_Parser_Base::c_string_vector Abstract_Class_Parser_Base::keys_;

//---------------------------------------------------------------------------//
Abstract_Class_Parser_Base::c_string_vector::~c_string_vector()
{
    unsigned const n = size();
    for (unsigned i=0; i<n; ++i)
    {
        delete[] operator[](i);
    }
}

} // end namespace rtt_parser


//---------------------------------------------------------------------------//
//              end of parser/Abstract_Class_Parser.cc
//---------------------------------------------------------------------------//
