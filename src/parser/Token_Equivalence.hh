//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   Token_Equivalence.hh
 * \author Kelly G. Thompson
 * \brief  Provide services for ApplicationUnitTest framework.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//

#ifndef rtt_Token_Equivalence_hh
#define rtt_Token_Equivalence_hh

#include "String_Token_Stream.hh"
#include "ds++/UnitTest.hh"

namespace rtt_parser {
//! Search tokens for existance of keyword.
void check_token_keyword(String_Token_Stream &tokens, std::string const &keyword,
                         rtt_dsxx::UnitTest &ut, unsigned const &occurance = 1);

//! Search tokens for keyword and compare associated value.  Report result to
//! UnitTest.
void check_token_keyword_value(String_Token_Stream &tokens, std::string const &keyword,
                               int const value, rtt_dsxx::UnitTest &ut,
                               unsigned const &occurance = 1);

//! Search tokens for keyword and compare associated value.  Report result to
//! UnitTest.
void check_token_keyword_value(String_Token_Stream &tokens, std::string const &keyword,
                               double const value, rtt_dsxx::UnitTest &ut,
                               unsigned const &occurance = 1);
} // namespace rtt_parser

#endif //  rtt_Token_Equivalence_hh

//--------------------------------------------------------------------//
// end of Token_Equivalence.hh
//--------------------------------------------------------------------//
