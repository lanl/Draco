//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/XGetopt.cc
 * \author Kelly Thompson, Katherine Wang
 * \date   Tuesday, Oct 27, 2015, 15:17 pm
 * \brief  Command line argument handling similar to getopt.
 * \note   Copyright (C) 2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "XGetopt.hh"
#include <algorithm>
#include <sstream>

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
/*!
 * \brief Convert a string into easy-to-parse vectors and strip colons.
 * \param[in] shortopts A string of character and optional colons. Each
 *            character represents one option that is registered. Characters
 *            followed by a colon represent options that require additional
 *            information to be read from the command line.
 * \return A vector of chars that represent all support options.
 *
 * - This function also populates the class member vshortopts_hasarg.
 * - Valid shortopts are single characters or a character followed by a colon.
 */
std::vector< char > XGetopt::decompose_shortopts( std::string const & shortopts )
{
    Require( shortopts.size() > 0 );
    std::vector< char > vso;
    for( size_t i = 0; i < shortopts.size(); ++i )
    {
        vso.push_back( shortopts[ i ] );
        if( i + 1 < shortopts.size() && shortopts[ i + 1 ] == std::string( ":" )[ 0 ] )
            this->vshortopts_hasarg[shortopts[ i ]] = true;
        else
            this->vshortopts_hasarg[shortopts[ i ]] = false;
    }
    Ensure( vso.size() > 0 && vso.size() <= shortopts.size() );
    return vso;
}


//---------------------------------------------------------------------------//
/*!
 * \brief Convert a map<char,string> into easy-to-parse vectors and strip
 *        colons.
 * \param[in] longopts A map<char,string> data structure that contains all
 *            supported options that are to be registered. Each entry pair
 *            represents one option that is registered, a short option that is
 *            used with a single dash (-h) and a long option that is used with a
 *            double dash (--help). Strings followed by a colon represent
 *            options that require additional information to be read from the
 *            command line.

 * \return A vector of chars that represent all support options.
 *
 * If the string value ends in a colon, ':', strip it from the map.  This data
 * has already been processed to indicate that an argument following the option
 * is expected.
 */
XGetopt::csmap XGetopt::store_longopts( csmap const & longopts_ )
{
    Require( longopts_.size() > 0 );
    csmap retValue;
    for( auto it=longopts_.begin(); it != longopts_.end(); ++it )
    {
        char        const key   = it->first;
        std::string const value = it->second;
        if( value[value.size()-1] == ':' )
            retValue[key] = value.substr(0,value.size()-1);
        else
            retValue[key] = value;
    }
    Ensure( retValue.size() == longopts_.size() );
    return retValue;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Convert long/short option map argument into state data for class.
 *
 * \param[in] longopts_ a map<char,string> that contains letters that represent
 *            valid command line arguments and associated long option names.
 * \return A vector<char> of all registered options
 *
 * - This function also populates the class member vshortopts_hasarg.
 * - Valid shortopts are single characters.
 * - Valid longopts are strings that are optionally followed by a colon to
 *   signify that a value is expected to follow the option.
 */
std::vector< char > XGetopt::decompose_longopts( csmap const & longopts_ )
{
    Require( longopts_.size() > 0 );
    std::vector< char > vso;
    for( auto it=longopts_.begin(); it != longopts_.end(); ++it )
    {
        char        const key   = it->first;
        std::string const value = it->second;
        Insist( std::find(vso.begin(),vso.end(),key) == vso.end(),
               "You cannot use the same single character command line argument more than once.");
        vso.push_back(key);
        if( value[value.size()-1] == ':' )
            this->vshortopts_hasarg[key] = true;
        else
            this->vshortopts_hasarg[key] = false;
    }
    Ensure( vso.size() == longopts_.size() );
    return vso;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Match provided command line arguments to registered options.
 *
 * Populates the member data:
 * - known_arguments
 * - known_arguments_values
 * - unknown_arguments
 */
void XGetopt::match_args_to_options()
{
    // Loop over all command line arguments
    for( size_t iarg=0; iarg < cmd_line_args.size(); ++iarg )
    {
        // stop processing command line arguments
        if( cmd_line_args[ iarg ] == std::string( "--" ) ) return;

        // try to match a short option ('-a')
        if( cmd_line_args[ iarg ].size() == 2 )
        {
            char const shortarg = cmd_line_args[ iarg ][1];
            if( std::find( vshortopts.begin(), vshortopts.end(), shortarg )
                != vshortopts.end() )
            {
                matched_arguments.push_back(cmd_line_args[ iarg ][1]);
                if( vshortopts_hasarg[shortarg] )
                {
                    bool cmd_line_option_requires_value(
                        iarg+1 < cmd_line_args.size() );
                    Insist( cmd_line_option_requires_value,
                            std::string("The command argument '-")
                            + shortarg +
                            "' expected a value to be provided.");
                    matched_arguments_values[shortarg]
                        = cmd_line_args[ ++iarg ];
                }
            }
            else
            {
                unmatched_arguments.push_back(cmd_line_args[ iarg ]);
            }
        }

        // consider string-based optons here.
        else if( cmd_line_args[ iarg ].substr( 0, 2 ) == std::string( "--" ) )
        {
            // The command line argument w/o the leading '--'
            std::string const longarg = cmd_line_args[ iarg ].
                                        substr(2,cmd_line_args[ iarg ].size());
            char shortarg('\0');

            // Find long argument match and its associated shortarg key.
            for( auto it=longopts.begin(); it!=longopts.end(); ++it )
            {
                if( it->second == longarg )
                {
                    shortarg = it->first;

                    // Save the match and associated data.
                    matched_arguments.push_back( shortarg );
                    if( vshortopts_hasarg[shortarg] )
                    {
                        bool cmd_line_option_requires_value(
                            iarg+1 < cmd_line_args.size() );
                        Insist( cmd_line_option_requires_value,
                                std::string("The command argument '--")
                                + longarg +
                                "' expected a value to be provided.");
                        matched_arguments_values[shortarg] =
                            cmd_line_args[ ++iarg ];
                    }
                    break;
                }
            }
            // if we get here, there was an argument that starts with '--' but
            // does not match a registered option:
            if( shortarg == '\0' )
                unmatched_arguments.push_back(cmd_line_args[ iarg ]);
        }

        // Provided argument is not '[-][A-z]' or does not start with '--' or is
        // not a value string associated with a preceeding argument.
        else
        {
            unmatched_arguments.push_back(cmd_line_args[ iarg ]);
        }
    }
    return;
}

// Comparator object that helps me find the longest string in a vector<string>
// structure.
struct strlenComparator
{
    strlenComparator( std::map< char, bool > const & vshortopts_hasarg_ )
        : vshortopts_hasarg( vshortopts_hasarg_ )
    { /* empty */ }
    std::map< char, bool > const vshortopts_hasarg;
    bool operator() (std::pair<char,std::string> const & a,
                     std::pair<char,std::string> const & b )
    {
        if( vshortopts_hasarg.at(b.first) )
            return false;
        else
            return (a.second).size() < (b.second).size();
    }
};

//---------------------------------------------------------------------------//
// File scope helper function.
// size_t get_max_len( XGetopt::csmap const & mymap, std::map<char,bool> include_mask; )
// {
//     size_t max_len=0;
//     for( auto it=mymap.begin(); it!=mymap.end();++it )
//     {
//         bool hasarg = vshortopts_hasarg.at(it->first);
//         if( hasarg && it->second.length() > max_len)
//             max_len = it->second.length();
//         }
// }

//---------------------------------------------------------------------------//
/*!
 * \brief Construct a help/use message that can be printed by the main program.
 *
 * \param[in] appName A string that will be printed as the name of the program.
 * \return A string that prings the program name and all registered options.
 */
std::string XGetopt::display_help( std::string const & appName ) const
{
    Require( appName.size() > 0 );

    // Build the output string starting with an empty ostringstream.
    std::ostringstream msg;
    msg << "\nUsage:"
        << "\n  " << appName << " [options]"
        << "\n\nOptions without arguments:\n";

    if( longopts.size() > 0 )
    {
        // find the longest longopt string:
        // size_t max_len = std::max_element(
        //     longopts.begin(), longopts.end(),
        //     strlenComparator(vshortopts_hasarg) )->second.size();

        // find the longest longopt string, excluding args with options:
        size_t max_len(0);
        for( auto it=longopts.begin(); it!=longopts.end();++it )
        {
            bool hasarg = vshortopts_hasarg.at(it->first);
            if( ! hasarg && it->second.length() > max_len)
                max_len = it->second.length();
        }

        // show options w/o arguments first
        for( auto it=longopts.begin(); it!=longopts.end();++it )
        {
            char        shortopt = it->first;
            std::string longopt  = it->second;
            bool        hasarg   = vshortopts_hasarg.at(shortopt);
            std::string helpmsg("");

            if( helpstrings.count(shortopt) > 0 )
                helpmsg = helpstrings.at(shortopt);

            // pad longopt so that they are all of the same length.
            if( longopt.length() < max_len )
                longopt = longopt + std::string(max_len-longopt.length(),' ');

            if( ! hasarg )
                msg << "   -" << shortopt << " | --" << longopt << " : "
                    << helpmsg << "\n";
        }

        // find the longest longopt string, including args with options:
        max_len=0;
        for( auto it=longopts.begin(); it!=longopts.end();++it )
        {
            bool hasarg = vshortopts_hasarg.at(it->first);
            if( hasarg && it->second.length() > max_len)
                max_len = it->second.length();
        }

        // show options that require arguments
        msg << "\nOptions requiring arguments:\n";
        for( auto it=longopts.begin(); it!=longopts.end();++it )
        {
            char        shortopt = it->first;
            std::string longopt  = it->second;
            bool        hasarg   = vshortopts_hasarg.at(shortopt);
            std::string helpmsg("");

            if( helpstrings.count(shortopt) > 0 )
                helpmsg = helpstrings.at(shortopt);

            if( hasarg )
            {
                // format the help string by replacing '\n' with
                // '\n'+hanging_indent.
                std::string hanging_indent(21+max_len,' ');
                std::string formatted_help_msg(helpmsg);
                size_t index(0);
                while (true)
                {
                    // find first line break.
                    index = formatted_help_msg.find('\n', index);
                    if (index == std::string::npos) break;
                    // add hanging indent.
                    formatted_help_msg.insert(++index, hanging_indent);
                    // Advance index forward so the next iteration doesn't pick
                    // it up as well.
                    index += 21;
                }

                msg << "   -" << shortopt << " | --" << longopt << " <value>"
                    << std::string(max_len-longopt.length(),' ') << " : "
                    << formatted_help_msg << "\n";
            }
        }
    }
    else
    {
        // show options w/o arguments first

        for( size_t i=0; i<vshortopts.size(); ++i )
        {
            if( vshortopts[i] == ':' ) continue;
            bool        hasarg   = vshortopts_hasarg.at(vshortopts[i]);
            std::string helpmsg("");
            if( helpstrings.count(vshortopts[i]) > 0 )
                helpmsg = helpstrings.at(vshortopts[i]);

            if( ! hasarg )
                msg << "   -" << vshortopts[i] << " : " << helpmsg << "\n";
        }

        // show options that require arguments

        msg << "\nOptions requiring arguments:\n";
        for( size_t i=0; i<vshortopts.size(); ++i )
        {
            char        shortopt = vshortopts[i];
            bool        hasarg   = vshortopts_hasarg.at(shortopt);
            std::string helpmsg("");
            if( helpstrings.count(shortopt) > 0 )
                helpmsg = helpstrings.at(shortopt);

            if( hasarg )
            {
                // format the help string by replacing '\n' with
                // '\n'+hanging_indent.
                std::string hanging_indent(16,' ');
                std::string formatted_help_msg(helpmsg);
                size_t index(0);
                while (true)
                {
                    // find first line break.
                    index = formatted_help_msg.find('\n', index);
                    if (index == std::string::npos) break;
                    // add hanging indent.
                    formatted_help_msg.insert(++index, hanging_indent);
                    // Advance index forward so the next iteration doesn't pick
                    // it up as well.
                    index += hanging_indent.size() ;
                }

                msg << "   -" << shortopt << " <value>" << " : "
                    << formatted_help_msg << std::endl;
            }
        }
    }

    msg << "   -- : Stop processing any further arguments." << std::endl;
    return msg.str();
}

} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
// end of XGetopt.cc
//---------------------------------------------------------------------------//
