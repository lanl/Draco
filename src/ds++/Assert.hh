//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file  ds++/Assert.hh
 * \brief Header file for Draco specific exception class definition (rtt_dsxx::assertion). Also
 *        define Design-by-Contract macros.
 * \note  Copyright (C) 2010-2021 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef RTT_dsxx_Assert_HH
#define RTT_dsxx_Assert_HH

#include "DracoTerminal.hh"
#include "ds++/config.h"
#include <stdexcept>
#include <string>

// This warning is issued because we chose to derive from std::logic_error (IS-A) instead of
// creating a class that 'HAS-A' std::logic_error.object.
#if defined(MSVC)
#pragma warning(push)
// non dll-interface class 'std::logic_error' used as base for dll-interface class
// 'rtt_dsxx::assertion'
#pragma warning(disable : 4275)
#endif

// DBC should be defined in ds++/config.h and can be set via cmake as DRACO_DBC_LEVEL=[0-31]
#if !defined(DBC)
#define DBC 7
#endif

namespace rtt_dsxx {

//================================================================================================//
/*!
 * \class assertion
 *
 * Exception notification class for Draco specific assertions.
 *
 * This class is derived from std::runtime_error.  In fact, this class provides no significant
 * change in functionality from std::runtime_error.  This class provides the following features in
 * addition to those found in std::runtime_error:
 *
 * 1. rtt_dsxx::assertion does provide an alternate constructor that allows us to automatically
 *    insert file name and line location into error messages.
 *
 * 2. It provides a specialized form of std::runtime_error.  This allows Draco code to handle Draco
 *    specific assertions differently from generic C++ or STL exceptions.  For example
 *
 * \code
 *    try
 *    {
 *       throw rtt_dsxx::assertion( "My error message." );
 *    }
 *    catch ( rtt_dsxx::assertion &a )
 *    {
 *       // Catch Draco exceptions first.
 *       cout << a.what() << endl;
 *       exit(1);
 *    }
 *    catch ( std::runtime_error &e )
 *    {
 *       // Catch general runtime_errors next
 *       cout << e.what() << endl;
 *    }
 *    catch ( ... )
 *    {
 *       // Catch everything else
 *        exit(1);
 *    }
 * \endcode
 *
 * \note Assertion should always be thrown as objects on the stack and caught as references.
 */
/*!
 * \example ds++/test/tstAssert.cc
 * Assertion and DBC examples.
 */
/*!
 * \section assertRevisionHistory Revision history:
 * \subsection date11March2003 11 March 2003 (Kelly Thompson)
 *
 * Major changes for class rtt_dsxx::assertion:
 *
 * Let rtt_dsxx::assertion derive from std::runtime_error.  G. Furnish intended for the design to
 * follow this model but the C++ compilers at time did not have full support for \c "stdexcept".
 * API is unchanged but the guts of the class are significantly different (i.e.: use copy and
 * assignment operators from base class, use string instead of \c const \c char \c *, etc.).
 *
 * \subsection date28July1997 28 July 1997
 *
 * Total rewrite of the ds++ Assertion facility.  In the old formulation, Assert and Insist expanded
 * into inline if statements, for which the body were blocks which both contained output to cerr,
 * and which threw an exception.  This has proven problematic because it meant that \c Assert.hh
 * included \c iostream.h, and also because throw was inlined (not a problem for Assert per se,
 * since Assert was generally shut off for optimized codes, but still a problem for Insist).  These
 * effects combined to present problems both at compile time and, possibly, at run time.
 *
 * The new formulation, motivated by suggestions from Arch Robison and Dave Nystrom, uses a call to
 * an out-of-line function to do the actual logging and throwing of exceptions.  This permits
 * removing \c iostream.h from the inclusion graph for translation units including \c Assert.hh, and
 * avoids inlining of exception throwing code.  Moreover, it also allows us to get \c __FILE__ into
 * play, which we couldn't do before because there was no way to do string concatenation with it,
 * and you couldn't do free store management effectively in an inlined macro.  By going out to a
 * global function, we pick up the ability to formulate a more complete picture of the error, and
 * provide some optimization capability (both compile time and run time).
 *
 * Note also that at this juncture, we go ahead and drop all support for compilers which are
 * incapable of compiling exception code.  From here forward, exceptions are assumed to be
 * available.  And thereby we do our part to promote "standard C++".
 *
 * \sa http://akrzemi1.wordpress.com/2013/01/04/preconditions-part-i/
 */
//================================================================================================//

class assertion : public std::logic_error {
public:
  /*!
   * \brief Default constructor for ds++/assertion class.
   *
   * This constructor creates a ds++ exception object.  This object is derived form
   * std::runtime_error and has identical functionality.  The principal purpose for this class is to
   * provide an exception class that is specialized for Draco.  See the notes on the overall class
   * for more details.
   *
   * \param msg The error message saved with the exception.
   */
  explicit assertion(std::string const &msg) : std::logic_error(msg) { /* empty */
  }

  /*!
   * \brief Specialized constructor for rtt_dsxx::assertion class.
   *
   * This constructor creates a ds++ exception object.  This object is derived form
   * std::runtime_error and has identical functionality.  This constructor is specialized for use by
   * Draco DbC commands (Require, Ensure, Check, and Insist).  It forms the error message from the
   * test condition and the file and line number of the DbC command.
   *
   * \param cond The expression that failed a DbC test.
   * \param file The source code file name that contains the DbC test.
   * \param line The source code line number that contains the DbC test.
   *
   * \sa \c --with-dbc[=level], Require, Ensure, Check, Insist (\ref Draco_DBC)
   */
  assertion(std::string const &cond, std::string const &file, int const line)
      : std::logic_error(build_message(cond, file, line)) { /* empty */
  }

  //! Helper function to build error message that includes source file name and line number.
  static std::string build_message(std::string const &cond, std::string const &file,
                                   int const line);
};

//------------------------------------------------------------------------------------------------//
// FREE NAMESPACE FUNCTIONS
//------------------------------------------------------------------------------------------------//

//! Throw a rtt_dsxx::assertion for Require, Check, Ensure.
[[noreturn]] void toss_cookies(std::string const &cond, std::string const &file, int const line);

[[noreturn]] void toss_cookies_ptr(char const *const cond, char const *const file, int const line);

//! Throw a rtt_dsxx::assertion if condition fails
void check_cookies(bool cond, char const *cond_text, char const *file, int line);

//! Print error w/o throw
void show_cookies(std::string const &cond, std::string const &file, int const line);

//! Throw a rtt_dsxx::assertion for Insist.
[[noreturn]] void insist(std::string const &cond, std::string const &msg, std::string const &file,
                         int const line);

//! Pointer version of insist
[[noreturn]] void insist_ptr(char const *const cond, char const *const msg, char const *const file,
                             int const line);

#if defined __CUDA_ARCH__ && defined USE_CUDA
/*! \brief A special version of insist that does not throw.  Useful for GPU code. \sa
 *         device/config.h.in */
__host__ __device__ inline void no_exception_insist(char const *const cond, char const *const msg,
                                                    char const *const file, int const line) {
  printf("Insist: %s, failed in %s, line %d.\n", cond, file, line);
  printf("The following message was provided: \"%s\"\n", msg);
  return;
}
#endif

#if DBC & 16
//! Check version of insist
void check_insist(bool cond, char const *const condstr, std::string const &msg,
                  char const *const file, int const line);

//! Check Pointer version of insist
void check_insist_ptr(bool cond, char const *const condstr, char const *const msg,
                      char const *const file, int const line);
#endif

//! Add hostname and pid to error messages.
std::string verbose_error(std::string const &message);

} // namespace rtt_dsxx

//------------------------------------------------------------------------------------------------//
/*!
 * \page Draco_DBC Using the Draco Design-by-Contract Macros
 *
 * \section ddbc Using the Draco Design-by-Contract Macros
 *
 * The assertion macros are intended to be used for validating preconditions which must be true in
 * order for following code to be correct, etc.  For example,
 *
 * \code
 * Assert( x > 0. );
 * y = sqrt(x);
 * \endcode
 *
 * If the assertion fails, the code should just bomb.  Philosophically, it should be used to ferret
 * out bugs in preceding code, making sure that prior results are within reasonable bounds before
 * proceeding to use those results in further computation, etc.
 *
 * These macros are provided to support the Design By Contract formalism.  The activation of each
 * macro is keyed off a bit in the DBC macro which can be specified on the command line:
 *
 \verbatim
 Bit     DBC macro affected
 ---     ------------------
 0      Require
 1      Check
 2      Ensure, Remember
 3      (nothrow versions of above macros) [8-15]
 4      (deferred exception throwing versions of above macros)
 \endverbatim
 *
 * So for instance, \c -DDBC=7 turns them all on, \c -DDBC=0 turns them all off, and \c -DDBC=1
 * turns on \c Require but turns off \c Check and \c Ensure.  The default is to have them all
 * enabled.  The fourth bit (3) can be toggled on to disable the C++ thrown exception while keeping
 * all of the DBC checks and messages active.
 *
 * The \c Insist macro is akin to the \c Assert macro, but it provides the opportunity to specify an
 * instructive message.  The idea here is that you should use Insist for checking things which are
 * more or less under user control.  If the user makes a poor choice, we "insist" that it be
 * corrected, providing a corrective hint.
 *
 * \note We provide a way to eliminate assertions, but not insists.  The idea is that \c Assert is
 * used to perform sanity checks during program development, which you might want to eliminate
 * during production runs for performance sake.  Insist is used for things which really must be
 * true, such as "the file must've been opened", etc.  So, use \c Assert for things which you want
 * taken out of production codes (like, the check might inhibit inlining or something like that),
 * but use Insist for those things you want checked even in a production code.
 *
 * Special code for CUDA.
 *
 * If __CUDA_ARCH__ (processing with nvcc) and USE_CUDA=ON, then alter the behavior of the DbC macros
 * because cuda code cannot throw.
 */
/*!
 * \def Require(condition)
 *
 * Pre-condition checking macro.  On when DBC & 1 is true.
 */
/*!
 * \def Check(condition)
 *
 * Intra-scope checking macro.  On when DBC & 2 is true.
 */
/*!
 * \def Ensure(condition)
 *
 * Post-condition checking macro.  On when DBC & 4 is true.
 */
/*!
 * \def Remember(code)
 *
 * Add code to compilable code.  Used in the following manner:
 * \code
 *     Remember (int old = x;)
 *     // ...
 *     Ensure (x == old);
 * \endcode
 * On when DBC & 4 is true.
 */
/*!
 * \def Insist(condition, message)
 *
 * Inviolate check macro.  Insist is always on.
 */
/*!
 * \def Insist_ptr(condition, message)
 *
 * Same as Insist, except that it uses char pointers, rather than strings.  This is more efficient
 * when inlined.
 */
//------------------------------------------------------------------------------------------------//

// clang-format off

//------------------------------------------------------------------------------------------------//
/* No-throw versions of DBC [8-15]
 *
 * Eventually, we want the DBC to work in GPU/Cuda code, but for now just disable DBC.
 */
//------------------------------------------------------------------------------------------------//
#if ( DBC & 8 ) || ( defined __CUDA_ARCH__ && defined USE_CUDA )

#if ( DBC & 1 ) && !( defined __CUDA_ARCH__ && defined USE_CUDA )
#define REQUIRE_ON
#define Require(c) if (!(c)) rtt_dsxx::show_cookies( #c, __FILE__, __LINE__ )
#else
#define Require(c)
#endif

#if ( DBC & 2 ) && !( defined __CUDA_ARCH__ && defined USE_CUDA )
#define CHECK_ON
#define Check(c) if (!(c)) rtt_dsxx::show_cookies( #c, __FILE__, __LINE__ )
#define Assert(c) if (!(c)) rtt_dsxx::show_cookies( #c, __FILE__, __LINE__ )
#else
#define Check(c)
#define Assert(c)
#endif

#if ( DBC & 4 ) && !( defined __CUDA_ARCH__ && defined USE_CUDA )
#define ENSURE_ON
#define Ensure(c) if (!(c)) rtt_dsxx::show_cookies( #c, __FILE__, __LINE__ )
#else
#define Ensure(c)
#endif

//------------------------------------------------------------------------------------------------//
// Always on
//------------------------------------------------------------------------------------------------//
#if ( defined __CUDA_ARCH__  && defined USE_CUDA )
#define Insist_device(c, m) if(!(c)) rtt_dsxx::no_exception_insist( #c, m, __FILE__, __LINE__)
#else
#define Insist_device(c,m) if (!(c)) rtt_dsxx::insist( #c, m, __FILE__, __LINE__ )
#endif

#define Insist(c,m) if (!(c)) rtt_dsxx::insist( #c, m, __FILE__, __LINE__ )
#define Insist_ptr(c,m) if (!(c)) rtt_dsxx::insist_ptr( #c, m, __FILE__, __LINE__ )

#elif DBC & 16

//------------------------------------------------------------------------------------------------//
// Regular (exception throwing) versions of DBC, but with the check deferred.  This eliminates
// numerous untestable branches for coverage analysis, but can be costly in run time.
//------------------------------------------------------------------------------------------------//

// Note that the !!(c) is needed to force the argument to a bool for certain objects, notably
// instances of the stream classes, which provide an operator!  but no direct conversion to bool.

// Also note that Bad_Case is always off in this selection.

#if DBC & 1
#define REQUIRE_ON
#define Require(c) rtt_dsxx::check_cookies(!!(c), #c, __FILE__, __LINE__)
#else
#define Require(c)
#endif

#if DBC & 2
#define CHECK_ON
#define Check(c) rtt_dsxx::check_cookies(!!(c), #c, __FILE__, __LINE__)
#define Assert(c) rtt_dsxx::check_cookies(!!(c), #c, __FILE__, __LINE__)
#define Bad_Case(m)
#else
#define Check(c)
#define Assert(c)
#define Bad_Case(m)
#endif

#if DBC & 4
#define ENSURE_ON
#define Ensure(c) rtt_dsxx::check_cookies(!!(c), #c, __FILE__, __LINE__)
#else
#define Ensure(c)
#endif

//------------------------------------------------------------------------------------------------//
// Always on
//------------------------------------------------------------------------------------------------//
#define Insist_device(c, m) rtt_dsxx::check_insist(!!(c), #c, m, __FILE__, __LINE__)
#define Insist(c, m) rtt_dsxx::check_insist(!!(c), #c, m, __FILE__, __LINE__)
#define Insist_ptr(c,m) rtt_dsxx::check_insist_ptr( !!(c), #c, m, __FILE__, __LINE__ )

#else // not DBC & 8 or DBC & 16

//------------------------------------------------------------------------------------------------//
// Regular (exception throwing) versions of DBC
//------------------------------------------------------------------------------------------------//
#if DBC & 1
#define REQUIRE_ON
#define Require(c) if (!(c)) rtt_dsxx::toss_cookies( #c, __FILE__, __LINE__ )
#else
#define Require(c)
#endif

#if DBC & 2
#define CHECK_ON
#define Check(c) if (!(c)) rtt_dsxx::toss_cookies( #c, __FILE__, __LINE__ )
#define Assert(c) if (!(c)) rtt_dsxx::toss_cookies( #c, __FILE__, __LINE__ )
#define Bad_Case(m) default: rtt_dsxx::toss_cookies( m, __FILE__, __LINE__ )
#else
#define Check(c)
#define Assert(c)
#define Bad_Case(m)
#endif

#if DBC & 4
#define ENSURE_ON
#define Ensure(c) if (!(c)) rtt_dsxx::toss_cookies( #c, __FILE__, __LINE__ )
#else
#define Ensure(c)
#endif

//------------------------------------------------------------------------------------------------//
// Always on
//------------------------------------------------------------------------------------------------//
#define Insist_device(c,m) if (!(c)) rtt_dsxx::insist( #c, m, __FILE__, __LINE__ )
#define Insist(c,m) if (!(c)) rtt_dsxx::insist( #c, m, __FILE__, __LINE__ )
#define Insist_ptr(c,m) if (!(c)) rtt_dsxx::insist_ptr( #c, m, __FILE__, __LINE__ )

#endif // DBC & 8

//------------------------------------------------------------------------------------------------//
// If any of DBC is on, then make the remember macro active and the NOEXCEPT inactive.
//------------------------------------------------------------------------------------------------//
#if DBC && !( defined USE_CUDA && defined __CUDA_ARCH__ )
#define REMEMBER_ON
#define Remember(c) c
#define NOEXCEPT
#define NOEXCEPT_C(c)
#else
#define Remember(c)
#define NOEXCEPT noexcept
#define NOEXCEPT_C(c) noexcept(c)
#endif

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Define a macro that disables potential exception throws for optimized (Release) code.
 *
 * Example:
 *
 * \code
 * double get_db(Vec3 const &, Vec3 const &) const ONLY_DBC_THROWS;
 * \endcode
 *
 * Issues:
 * - C++11 - Dynamic exception specifications are deprecated until C++17 except on lambda-declarator
 *           or on a function declarator that is the top-level (until C++17) declarator of a
 *           function, variable, or non-static data member, whose type is a function type, a pointer
 *           to function type, a reference to function type, a pointer to member function type. It
 *           may appear on the declarator of a parameter or on the declarator of a return type.  \sa
 *           http://en.cppreference.com/w/cpp/language/except_spec */
//------------------------------------------------------------------------------------------------//

#if DBC
#define ONLY_DBC_THROWS
#else
#define ONLY_DBC_THROWS noexcept
#endif

//------------------------------------------------------------------------------------------------//
/*! Ensure all possible Draco DbC macros are defined.  If not already defined, then define them as
 *  no-op. */
//------------------------------------------------------------------------------------------------//
#ifndef Require
#define Require(c)
#endif
#ifndef Check
#define Check(c)
#endif
#ifndef Assert
#define Assert(c)
#endif
#ifndef Ensure
#define Ensure(c)
#endif
#ifndef Insist
#define Insist(c,m)
#endif
#ifndef Insist_device
#define Insist_device(c,m)
#endif
#ifndef Insist_ptr
#define Insist_ptr(c,m)
#endif
#ifndef Bad_Case
#define Bad_Case(m)
#endif
#ifndef Remember
#define Remember(c)
#endif

//------------------------------------------------------------------------------------------------//
#if defined(MSVC)
#pragma warning(pop)
#endif

// clang-format on

#endif // RTT_dsxx_Assert_HH

//------------------------------------------------------------------------------------------------//
// end of ds++/Assert.hh
//------------------------------------------------------------------------------------------------//
