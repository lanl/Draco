//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstAbstract_Class_Parser.cc
 * \author Kent G. Budge
 * \date   Tue Nov  9 14:34:11 2010
 * \brief  Test the Abstract_Class_Parser template
 * \note   Copyright (C) 2006-2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "parser/Class_Parse_Table.hh"
#include "parser/Abstract_Class_Parser.hh"
#include "parser/utilities.hh"
#include "parser/File_Token_Stream.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_parser;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
/*
 * The following would typically be declared in a file associated with the
 * Parent class.
 */

class Parent
{
  public:

    virtual ~Parent(){}

    virtual string name() = 0;
};

//---------------------------------------------------------------------------//
namespace rtt_parser
{
template<>
class Class_Parse_Table<Parent>
{
  public:

    // TYPEDEFS

    typedef Parent Return_Class;

    // MANAGEMENT

    Class_Parse_Table()
    {
        child_.reset();
        current_ = this;
    }

    // SERVICES

    Parse_Table const &parse_table() const { return parse_table_; }

    bool allow_exit() const { return false; }

    void check_completeness(Token_Stream &tokens)
    {
        if (child_ == SP<Parent>())
        {
            tokens.report_semantic_error("no parent specified");
        }
    }

    SP<Parent> create_object()
    {
        return child_;
    }

    // STATICS

    static void register_model(string const &keyword,
                               SP<Parent> parse_function(Token_Stream&) )
    {
        Abstract_Class_Parser<Parent,
                              get_parse_table_,
                              get_parsed_object_>::register_child(keyword, parse_function);
    }

  private:

    // IMPLEMENTATION

    static Parse_Table &get_parse_table_()
    {
        return parse_table_;
    }

    static SP<Parent> &get_parsed_object_()
    {
        return child_;
    }

    // DATA

    static SP<Parent>         child_;
    static Class_Parse_Table *current_;
    static Parse_Table        parse_table_;
};

SP<Parent>                 Class_Parse_Table<Parent>::child_;
Class_Parse_Table<Parent> *Class_Parse_Table<Parent>::current_;
Parse_Table                Class_Parse_Table<Parent>::parse_table_;

//---------------------------------------------------------------------------//
/*
 * Specialization of the parse_class function for T=Parent
 */
template<>
SP<Parent> parse_class<Parent>(Token_Stream &tokens)
{
    return parse_class_from_table<Class_Parse_Table<Parent> >(tokens);
}

} // namespace rtt_parser

//---------------------------------------------------------------------------//
/*
 * The following is what you would expect to find in a file associated with
 * the Son class.
 */
class Son : public Parent
{
  public:

    virtual string name(){ return "son"; }

    Son(double /*snip_and_snails*/){}
};

double parsed_snips_and_snails;

void parse_snips_and_snails(Token_Stream &tokens, int)
{
    if (parsed_snips_and_snails>=0.0)
    {
        tokens.report_semantic_error("snips and snails already specified");
    }

    parsed_snips_and_snails = parse_real(tokens);
    if (parsed_snips_and_snails<0.0)
    {
        tokens.report_semantic_error("snips and snails must not be "
                                     "negative");
        parsed_snips_and_snails = 2;
    }
}

namespace rtt_parser
{
template<>
class Class_Parse_Table<Son>
{
  public:

    // TYPEDEFS

    typedef Son Return_Class;

    // MANAGEMENT

    Class_Parse_Table()
    {
        if (!parse_table_is_initialized_)
        {
            const Keyword keywords[] =
                {
                    {"snips and snails", parse_snips_and_snails, 0, ""},
                };

            const unsigned number_of_keywords = sizeof(keywords)/sizeof(Keyword);
            parse_table_.add(keywords, number_of_keywords);

            parse_table_is_initialized_ = true;
        }

        parsed_snips_and_snails = -1;
        current_ = this;
    }

    // SERVICES

    bool allow_exit() const { return false; }

    Parse_Table const &parse_table() const { return parse_table_; }

    void check_completeness(Token_Stream &tokens)
    {
        if (parsed_snips_and_snails == -1)
        {
            tokens.report_semantic_error("no snips and snails specified");
        }
    }

    SP<Son> create_object()
    {
        SP<Son> Result(new Son(parsed_snips_and_snails));
        return Result;
    }

  private:

    // STATIC

    static Class_Parse_Table *current_;
    static Parse_Table parse_table_;
    static bool parse_table_is_initialized_;
};

//---------------------------------------------------------------------------//
Class_Parse_Table<Son> *Class_Parse_Table<Son>::current_;
Parse_Table             Class_Parse_Table<Son>::parse_table_;
bool Class_Parse_Table<Son>::parse_table_is_initialized_ = false;

//---------------------------------------------------------------------------//
template<>
SP<Son> parse_class<Son>(Token_Stream &tokens)
{
    return parse_class_from_table<Class_Parse_Table<Son> >(tokens);
}

} // end namespace rtt_parser


//---------------------------------------------------------------------------//
/*
 * The following is what you would expect to find in a file associated with
 * the Daughter class.
 */

class Daughter : public Parent
{
  public:

    virtual string name(){ return "daughter"; }

    Daughter(double /*sugar_and_spice*/){}
};

double parsed_sugar_and_spice;

void parse_sugar_and_spice(Token_Stream &tokens, int)
{
    if (parsed_sugar_and_spice>=0.0)
    {
        tokens.report_semantic_error("sugar and spice already specified");
    }

    parsed_sugar_and_spice = parse_real(tokens);
    if (parsed_sugar_and_spice<0.0)
    {
        tokens.report_semantic_error("sugar and spice must not be "
                                     "negative");
        parsed_sugar_and_spice = 2;
    }
}

namespace rtt_parser
{
template<>
class Class_Parse_Table<Daughter>
{
  public:

    // TYPEDEFS

    typedef Daughter Return_Class;

    // MANAGEMENT

    Class_Parse_Table()
    {
        if (!parse_table_is_initialized_)
        {
            const Keyword keywords[] =
                {
                    {"sugar and spice", parse_sugar_and_spice, 0, ""},
                };

            const unsigned number_of_keywords = sizeof(keywords)/sizeof(Keyword);
            parse_table_.add(keywords, number_of_keywords);
            parse_table_is_initialized_ = true;
        }

        parsed_sugar_and_spice = -1;
        current_ = this;
    }

    // SERVICES

    bool allow_exit() const { return false; }

    Parse_Table const &parse_table() const { return parse_table_; }

    void check_completeness(Token_Stream &tokens)
    {
        if (parsed_sugar_and_spice == -1)
        {
            tokens.report_semantic_error("no sugar and spice specified");
        }
    }

    SP<Daughter> create_object()
    {
        SP<Daughter> Result(new Daughter(parsed_sugar_and_spice));
        return Result;
    }

  private:

    // STATIC

    static Class_Parse_Table *current_;
    static Parse_Table parse_table_;
    static bool parse_table_is_initialized_;
};

Class_Parse_Table<Daughter> *Class_Parse_Table<Daughter>::current_;
Parse_Table Class_Parse_Table<Daughter>::parse_table_;
bool Class_Parse_Table<Daughter>::parse_table_is_initialized_;

//---------------------------------------------------------------------------//
template<>
SP<Daughter> parse_class<Daughter>(Token_Stream &tokens)
{
    return parse_class_from_table<Class_Parse_Table<Daughter> >(tokens);
}

} // end namespace rtt_parser

/* the followingn would typically live in some client file for Parent. */

SP<Parent> parent;

//static
void parse_parent(Token_Stream &tokens, int)
{
    parent = parse_class<Parent>(tokens);
}

const Keyword top_keywords[] =
{
    {"parent", parse_parent, 0, ""},
};

const unsigned number_of_top_keywords = sizeof(top_keywords)/sizeof(Keyword);

Parse_Table top_parse_table(top_keywords, number_of_top_keywords);

SP<Parent> parse_son(Token_Stream &tokens)
{
    return parse_class<Son>(tokens);
}

SP<Parent> parse_daughter(Token_Stream &tokens)
{
    return parse_class<Daughter>(tokens);
}

//---------------------------------------------------------------------------//
/* Test all the above */

void test(UnitTest &ut)
{
    Class_Parse_Table<Parent>::register_model("son",
                                              parse_son );

    Class_Parse_Table<Parent>::register_model("daughter",
                                              parse_daughter );

    // Build path for the input file
    string const sadInputFile( ut.getTestSourcePath()
			       + std::string("sons_and_daughters.inp") );

    File_Token_Stream tokens( sadInputFile );

    top_parse_table.parse(tokens);

    cout << parent->name() << endl;

    if (tokens.error_count()==0 &&
        parent!=SP<Parent>() &&
        parent->name()=="son")
    {
        PASSMSG("Parsed son correctly");
    }
    else
    {
        FAILMSG("Did NOT parse son correctly");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try { test(ut); }
    UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstAbstract_Class_Parser.cc
//---------------------------------------------------------------------------//
