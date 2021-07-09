# "Why we did it that way": Parsers

* _Kent G. Budge_

## Introduction

In this, the first installment of "Why we did it that way", I will discuss why the parser framework for our experimental drivers was written the way it was, what changes were made based on experience, and how we might change it in the future.

## Motivation

I desired to allow our test drivers to read input decks that were as human-friendly as possible. To me, this meant that they should be insensitive to white space, should use a reasonably intuitive grammar, and should not throw surprises at the user. I also wanted the parser system to be reasonably easy to maintain and to match to the data structures used to contain a problem specification.

I was influenced in this by my experience with various programming languages. I did not want the kind of sensitivity to white space evident in older versions of Fortran or in input files read by older programs. I wanted something more C-like in its treatment of white space. I was also influenced by my readings on compiler theory, which is well developed for scanners and parsers.

I also wished to avoid dependency on outside tools or libraries. This led me to rule out using flex and bison to automatically generate parsers. At the time, C++ support in these tools was quite poor (it still isn't excellent) and that was another consideration.

I drew on my previous experience with developing parsers for input files for the Alegra code at Sandia. As a result, in its original set of features, the Capsaicin parser library was quite similar to the Alegra parser library.

## The original design

Compiler developers make a distinction between a scanner, which reads the raw characters of a file and breaks them into tokens, and parsers, which take a stream of tokens and parses them into grammatical structures. I adopted this distinction, with `Token_Stream` and its children representing different flavors of scanners, and `Parse_Table` representing a particular kind of parser.

### The scanner

`Token_Stream` takes a stream of characters, provided by virtual functions implemented in its child classes, and turns these into tokens. Originally, the token types were integers, real numbers, character strings, keywords, EXIT, ERROR, and OTHER. Integer and real tokens are similar to their C/C++ counterparts. Manifest strings were limited to double quote (") delimited strings with no escape characters, so they were a simplified version of string constants in C/C++. The EXIT token was needed to indicate end of file. The ERROR token was reserved for mysterious errors in reading the character stream, such as file system I/O failures. Any character not conforming to one of these token types, and not defined as a whitespace character, was returned as OTHER. Each token, in addition to being labeled with a type, carried its associated text, such as "37" for an integer or "3.1416e0" for a real number.

I chose to define keywords as sequences of C/C++ identifiers, such as

    gray source boundary

A keyword can contain digits, but not as the first character of an identifier within the keyword. I chose to define a keyword as a _sequence_ of identifiers mostly as a matter of taste; I preferred

    gray source boundary
to

    gray_source_boundary
However, this introduced the difficulty of how to separate adjacent keywords.

I solved this problem with my definitions of whitespace. I chose to distinguish _breaking_ and _nonbreaking_ whitespace. Nonbreaking whitespace, which consists of spaces and tabs, is permitted between identifiers within a single keyword. Breaking whitespace was defined as things like line breaks, commas, colons, semicolons, and the equals sign. Thus

    gray source boundary
is a single identifier while

     gray, source boundary
or

    gray = source boundary
or

    gray
    boundary source

are all scanned as two keywords, `gray` followed by `source boundary`. The comma or equals sign is discarded as whitespace. Likewise, whitespace between identifiers is reduced by the scanner to a single space character, so

    gray        source       boundary

becomes

    gray source boundary

after scanning. This means that a parser needs only deal with one possible "spelling" of a keyword.

My choice of whitespace characters also meant that

    time step 1.0
    time step, 1.0
    time    step = 1.0

all had exactly the same meaning.  I added options to `Token_Stream` to allow a client to require the scanner to treat _all_ white space as breaking (so that the definition of a keyword becomes identical to a single C/C++ identifier), or to choose a different set of nonbreaking whitespace characters. To my knowledge, no one has ever made use these, but they're still part of the class.

Likewise, the scanner is case-sensitive as written, but there are provisions for making the scanner case-insensitive. Again, this capability has never actually been used.

Early in development, I added an additional token, the END token. This was added to give the `Parse_Table` class a distinctive token to flag the end of a keyword block. The token is spelled as `end` so it looks like a keyword, but it has special semantic meaning. It is, in effect, the only reserved keyword in the Capsaicin parser library.

I also eventually added common C/C++ operators, such as `>=`, as special cases of OTHER. This was to allow a `Token_Stream` to support expression parsing.

With `Token_Stream` abstracting a string of characters that are to be converted to grammatical tokens, the concrete children of `Token_Stream` provide the stream of characters in different ways, and also provide different ways of reporting errors to the client.

The most straightforward is `File_Token_Stream`, which is takes its stream of characters from a disk file read in serial. A single MPI rank can attach a `File_Token_Stream` to a specified file with no reference to any other MPI rank, and scan the file asynchronously.

It's possible for all MPI ranks to open the same input file and read it independently using a `File_Token_Stream`, but if all MPI ranks are reading the same input file (as is typically the case for the main input file provided to a test driver), you could end up with 10,000 ranks simultaneously trying to read the same file. This will almost certainly cause serious problems. For the case where all ranks are reading the same input synchronously, I wrote the `Parallel_File_Token_Stream` class. Only rank 0 actually opens the disk file and reads its contents; the characters read from the file are then communicated to all other processors using MPI calls. Each processor then independently (but synchronously) scans the characters to convert them to tokens. This avoids having to communicate any structure more complicated that a block of characters across MPI. In particular, I don't have to communicate either scanned tokens or problem specification structures parsed from them.

I added `Console_Token_Stream` a bit later to support reading a stream of characters directly from the console.  I also found it useful, particularly for test code, to define `String_Token_Stream` for which the string of characters is usually defined explicitly at the point of construction.

I later inserted `Text_Token_Stream` between `Token_Stream` and all its (present) children, and moved everything associated with converting a stream of characters to a stream of tokens into `Text_Token_Stream`. `Token_Stream` is now a very simple abstract class representing a stream of tokens, without reference to how they are obtained. I don't recall what I had in mind with this, unless it was some half-baked thoughts about supporting a GUI from which a parser could read its tokens. But now that the distinction is there, we might as well keep it, since it's a well-defined distinction that has not proven problematic (even if we do not currently take advantage of it.)

When a client (parsing code) is handed a `Token_Stream`, it has essentially three operations it can perform on the token stream. It can do a `shift()`, which has the effect of popping the next token off the stream and returning it to the client. The client can do a `lookahead()`, which returns the next token on the stream without popping it off the stream. These strange names (shift and lookahead) are straight out of compiler theory.

The client can also report an error to the stream. It may seem odd for the source of input to also be the destination for error reports, but it reflects the fact that different sources of tokens require different forms of error reporting. With `File_Token_Stream`, the error goes to console in the obvious way, along with information on which line in the input file triggered the error. (This is provided by the `File_Token_Stream`, so the client doesn't have to worry about it.) With `Parallel_File_Token_Stream`, only processor 0 reports the error to console, since all processors are parsing the same input synchronously and will detect the same errors. With `Console_Token_Stream`, the error is reported to console and then the entire token stream is flushed -- the guy typing stuff live into his terminal gets to try again.

I later added capability to push a token back onto the stream, though I don't think it's ever been used in our drivers. Likewise, there is capability to look more than one token ahead, using (for example) `lookahead(3)`  to see the fourth token presently in the stream. (`lookahead(0)` is identical to `lookahead()` and we're counting like mathematicians: 0, 1, 2, 3 ...) So far as I know, the lookahead more than one token has never been used in our drivers. (Lookahead at the next token gets used all the time.)

### Parsing

Programming language grammar is a weird and wonderful field in which I've occasionally dabbled. C has a slightly weird version of what is called an LALR(1) grammar, which is a large family of grammars for which one can conveniently generate a parser using the original `yacc`. It's slightly weird because of the so-called "scanner trick", which turns a grammar that is not strictly LALR(1) into something that yacc can parse as it it _was_ LALR(1). C++ infamously departed from LALR(1) with its very large number of ambiguous grammatical constructs, and so has what is called a GLR grammar. GLR is basically any grammar that is parsed from left to right.  It turns out `bison` has an option permitting it to generate a parser for any GLR language, but you have to add all kinds of stuff to decide which of several possible ways of parsing a construct is the right one.

#### Parse_Table

I decided early on that we would organize at least the top level of our parsers as keyword-based parsers. That is, every construct in our input grammar would begin with a keyword that would immediately identify the construct. This is not terribly restrictive; there are real general-purpose programming languages (Eiffel comes to mind) which are close to keyword-based. So a `Parse_Table` would operate by scanning the next token off the stream, which would always be a keyword if the token stream was well formed. The keyword would be sufficient to identify a particular parse function that would be called to complete the parsing of the construct. So if the parser saw

    sn
it would immediately call `parse_sn` which would parse the rest of the Sn specification. A parser of this kind would simply do a lookup on the keyword to find the right parse function to call.

I left unspecified, and it remains unspecified, just what each parse function called by a keyword would do. The parse function is passed the token stream, and it is allowed to pull off as many additional tokens as it wishes and to do whatever it wants with them.  Thus a programmer can insert a very general parser underneath the opening keyword, hand-crafted however he likes -- though this is not usual practice and is not really recommended.

#### Expressions

At some point, I decided it would be useful in some contexts to specify something as an expression, including expressions in space or time (for space- or time-dependent sources, for example.) I chose a C/C++ expression grammar, for which techniques for efficient parsing are very well established. Thus

    source boundary, boundary flag = 1, rate = 1.0 + 3.0*t
would specify a time-dependent source boundary for a test problem. The `t` in this context now has builtin semantics ("time"), but it is not reserved; you can still have `t` as a keyword in a parse table. The context determines now it is interpreted.

#### Utilities

I found it useful to define a number of classes like

    unsigned parse_unsigned_integer(Token_Stream &);
that would expect a particular kind of construct (in this case, an integer) and automatically report an error if the next token was not an INTEGER or if its value was negative. Likewise

    double parse_real(Token_Stream &);
would accept either an INTEGER or a REAL token, casting the former to double if necessary.

This led fairly naturally to functions like

    bool at_real(Token_Stream &);
to check whether the next token on the `Token_Stream` can be read as a real number.

Another wrinkle was

    double parse_quantity(Token_Stream &, Unit &);
which served two purposes. First, it allowed an input file to contain quantities with specified units, such as

    length = 2.3 cm
which would be automatically converted to Capsaicin's internal units (SI by default).  Second, it reported an error if a quantity specified this way had the wrong units. Detecting user input errors as early as possible is best practice.

We did not use the expression parser for parsing unit expressions, because we wanted to be able to write things like

    action = 3.9 J-s

corresponding to how such things are typically written by humans. The unit expression here is not joules minus seconds (which would be ill-formed) but joule-seconds, the unit of action. So we wrote a special expression parser for unit expressions, which also recognizes common unit names. So both `J` and `s` have special semantics in a unit expression, but not elsewhere.

Finally,

    real parse_temperature(Token_Stream &);
allowed a temperature to be specified in _either_ conventional temperature units or in energy units, so

    initial T = 11605 K
and

    initial T = 1 eV
are both supported in an input file and mean roughly the same thing. The second is preferred if you really mean 1 eV since the value will be taken from the Draco table of physical constants.

#### Class_Parse_Table

I came to Los Alamos in 2002 with a strong prejudice against templates. They were still relatively new and quite buggy in C++ at that time. Different compilers instantiated templates in different ways and keeping it all straight on multiple platforms was just a pain in the posterior.

Within a few years, templates became much more reliable and much better standardized. The standard template library, which was still very new in 2002, did a lot to help drive this. We began looking for opportunities to improve the code base by using templates intelligently.

This led me to define a `parse_class<Classname>` template that would automate a lot of the work of creating a data structure to hold the information from a parsed keyword block, parse the block, and return the data structure to the user. This was based off of `Parse_Table` but with a lot of standardization of how things were done.

The `parse_class<Classname>` template has a default implementation that makes use of a `Class_Parser<Classname>` which in turn uses `Parse_Table`. `Class_Parser<Classname>` can be defined how a developer wishes, as long as it provides the members expected by `parse_class<Classname>`.  The full documentation is in the class header.

Much of the complexity of parsing at present is that the more recent refinements to the parse library have never been retrofitted to older classes. This means that "best practice" may not be clear to new developers. However, it's also true that I've not been able to make the `Class_Parser<Classname>` approach as lightweight as I would like.

#### Abstract_Class_Parse_Table

This is probably the most difficult part of the library. It was written to address the question of how to parse variants of an abstract class. For example,

    transport
      sn
alerts the top level parser in our drivers that we are about to specify a transport model (in this case, Sn). We did not want to go with just

    sn
at the top level because this meant that, if we added a new transport model, we'd have to go in and modify the top level parse table for every driver that might use that class. By nesting the concrete model within an abstract wrapper, we have a function for parsing a transport model that is the only thing that has to change when a new model is added. All drivers will automatically get the new model.

I may have outsmarted myself on this one. I wanted it make it possible for a developer to add a local experimental model without having to change the shared code base to include his new model keyword. So I added a registration function to allow a client to add a new model. All the messiness of `Abstract_Class_Parse_Table` comes from this cool, but ultimately not really needed, requirement.  (Yes, I recognize that 'not really needed requirement' is an oxymoron.)

#### Unit systems

Capsaicin originally converted all input to SI units and used these consistently in its calculations. This proved difficult when debugging problems that were more naturally expressed in other unit systems (such as cgs). We added the ability for Capsaicin to pick an arbitrary unit system to use consistently, with SI as the default. This was specified in input files with a leading

    unit system, cgs-K
This had to appear before any other input. Otherwise we'd have to retroactively convert input already parsed, which was a horrible thought. We also added the option, strictly for machine-generated decks, to allow quantities normally requiring units to use implied units as specified by the leading `unit system` keyword. Users really should not do this!

## Move to Draco

The scanner and parser classes we originally wrote for Capsaicin are now in Draco. This is because we wrote other classes, such as quadrature classes and material models, that were also excellent candidates for shared projects but which came with their own parsers based on the Capsaicin library.

## Our parsers and host codes

Our parser library was written for use within our own project, in our test drivers. The test drivers are essential for performing the large number of regression tests we run nightly and weekly and before merging pull requests, and are also essential for our ongoing research work in transport methods.  This makes the parser library an essential part of our code base.

However, we use the parsers even in our flat interfaces, which may not seem like an obvious thing to do. Rather than simply construct the objects we need directly with a constructor, we tend to use the parameters passed in by the host code to build an input deck internally as a `String_Token_Stream`, then parse this to actually create the objects.

There are two reasons for this. The first is that, if there are semantic errors in the parameters passed in by the host, these are processed the same way we would process them in our drivers. We don't have to duplicate that diagnostic code. The second reason for using our parsers in flat interfaces is that it gives us the option of dumping out the internal input deck. We can then inspect it to see exactly what our code thinks it's been asked to solve. We could also feed the deck directly to our test drivers for debugging.

Is this really a good idea? Randy Baker thought so, and I think the case can certainly be made. But it adds  complexity, and it means we have to pull more object files into a linked executable.

## The future

> _Greetings, my friend. We are all interested in the future, for that is where you and I are going to spend the rest of our lives. And remember my friend, future events such as these will affect you in the future._
> -- Plan 9 from Outer Space
>
Our parsers are not as simple as we would like. It would be nice if adding a new integer parameter (for example) to a keyword block was as simple as adding a line like

    {"new parameter", integer, "this new parameter does x"}
to a table somewhere, and then

    int new_parameter;
was magically added to the class that holds the parse specification and the code to use ``parse_integer`` to parse this new parameter was magically added to the code base. I say "magically" because I am not confident this can be done entirely within C++. Perhaps (we are smart) but there might be no way to accomplish this other than to go outside the language and script the generation of parser code. Then there's the issue of what to do with parameters that are themselves structures based on a keyword block, and the issue of how to enforce semantics in a way that is informative to the user.  It's possible what we have now is the best compromise between maintainability, flexibility, and capability. Life is full of trade offs.
