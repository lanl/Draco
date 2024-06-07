#--------------------------------------------*-cmake-*---------------------------------------------#
# file   config/generate_dll_declspec.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2016 Feb 13
# brief  Generate dll_declspec.h used to define DLL_PUBLIC_<pkg> definitions.
# note   Copyright (C) 2015-2024 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

include_guard(GLOBAL)

# Generate a file used to define DLL_PUBLIC_pkg> definitions
function(generate_dll_declspec dir components)

  string(REPLACE "+" "x" safedir ${dir})

  set(dll_declspec_content
      "/*------------------------------------------*-C-*------------------------------------------*/
/*!
 * file   dll_declspec.h
 * brief  Defined macros that are used as declarators to control dllexport
 *        or dllimport linkage for dll files.
 * note   Copyright (C) 2015-2024 Triad National Security, LLC., All rights reserved.
 *
 * Provide toolset for hiding symbols from shared libraries.  By default MSVC
 * hides everying unless they are marked as 'dllexport'.  By default
 * Linux-based linkers export everything unless they are marked as
 * 'hidden'. The macros DLL_PUBLIC_<pkg> and DLL_LOCAL_<pkg> are provided to
 * help make symbol export portable between Linux and Windows.
 *
 * - http://gcc.gnu.org/wiki/Visibility
 * - https://msdn.microsoft.com/en-us/library/3y1sfaz2.aspx
 * - http://www.cmake.org/Wiki/BuildingWinDLL
 *
 * For MSVC, you must also mark functions as dllexport or dllimport.  When
 * compiling a file to place it in a dll, the symbols must be marked
 * 'dllexport', but when including a header that is associated with an
 * external .dll (include Array.hh and linking against -ldsxx.so), the same
 * symbol must be marked as 'dllimport.' The macro '<target>_EXPORTS>' will
 * be used to control how the symbols are marked.  This CPP definition is set
 * by cmake.
 *
 * Example use of DLL_PUBLIC_foo CPP Symbol:
 *
 * extern \"C\" DLL_PUBLIC_foo void function(int a);
 * class DLL_PUBLIC_foo SomeClass
 * {
 *    int c;
 *    DLL_LOCAL_foo void privateMethod();  // Only for use within this DSO
 *  public:
 *    Person(int _c) : c(_c) { }
 *    static void foo(int a);
 * };
 *
 * // dsxx/isFinite_pt.cc file:
 * template DLL_PUBLIC_foo
 * bool isInfinity( float const & x );
 */
/*------------------------------------------------------------------------------------------------*/

#ifndef rtt_${safedir}_config_h
#error \"Do not call this file directly. Call ${dir}/config.h instead.\"
#endif

#ifndef rtt_${safedir}_declspec_h
#define rtt_${safedir}_declspec_h")

  if(CMAKE_COMPILER_IS_GNUCXX)
    set(dll_declspec_content
        "${dll_declspec_content}

#define CMAKE_COMPILER_IS_GNUCXX 1

")
  endif()

  set(dll_declspec_content
      "${dll_declspec_content}

#ifdef DRACO_SHARED_LIBS /* building shared or static libs? */

  /* Windows ---------------------------------------- */
  #if defined _WIN32 || defined __CYGWIN__
")

  foreach(pkg ${components})
    set(dll_declspec_content
        "${dll_declspec_content}
    #ifdef Lib_${pkg}_EXPORTS /* building dll or importing symbols? */
      #ifdef __GNUC__
        #define DLL_PUBLIC_${pkg} __attribute__((dllexport))
      #else
         /* Note: actually gcc seems to also supports this syntax. */
        #define DLL_PUBLIC_${pkg} __declspec(dllexport)
      #endif /* __GNUC__ */
    #else
      #ifdef __GNUC__
        #define DLL_PUBLIC_${pkg} __attribute__((dllimport))
      #else
        /* Note: actually gcc seems to also supports this syntax. */
        #define DLL_PUBLIC_${pkg} __declspec(dllimport)
      #endif /* __GNUC__ */
    #endif
")
  endforeach()

  set(dll_declspec_content "${dll_declspec_content}
    #define DLL_LOCAL")

  # set( dll_declspec_content "${dll_declspec_content} #else /* DRACO_SHARED_LIBS */ /* when
  # building static libraries, no need to export symbols. */")

  # foreach(pkg ${components} ) set( dll_declspec_content "${dll_declspec_content} #define
  # DLL_PUBLIC_${pkg}") endforeach()

  # set( dll_declspec_content "${dll_declspec_content} #define DLL_LOCAL

  # #endif /* DRACO_SHARED_LIBS */")

  set(dll_declspec_content
      "${dll_declspec_content}

  /* Linux ---------------------------------------- */
  #else
    #if defined CMAKE_COMPILER_IS_GNUCXX && __GNUC__ >= 4
/* GCC's attributes don't work exactly like MSVC's declspec for the
way we implement explicit instantiations.  In Jayenne, we often have
class templates defined in the package directory and the explicit
instantiations are done in the test directory. Neither GCC nor MSVC
expect the instantiations to be defined in a library other than the
main package where the class template is declared.  Under this model,
MSVC requires the explicit instantiation to be exported from the test
library but GCC ignores these attributes (with a warning) because
declaration does not have the same attributes.  Because of this
distiction, I am disabling setting attributes when GCC is used.
(KT 2016-03-04) */
")

  foreach(pkg ${components})
    set(dll_declspec_content
        "${dll_declspec_content}
      // #define DLL_PUBLIC_${pkg} __attribute__ ((visibility(\"default\")))
      #define DLL_PUBLIC_${pkg} ")
  endforeach()

  set(dll_declspec_content
      "${dll_declspec_content}
      #define DLL_LOCAL  __attribute__ ((visibility(\"hidden\")))
      /* These defines are taken from dsxx/Compiler.hh written by Paul Henning */
      /* #define HIDE_FUNC __attribute__ ((visibility (\"hidden\"))) */
      /* #define EXPORT_FUNC __attribute__ ((visibility (\"default\"))) */
      /* #define HIDE_CLASS __attribute ((visibility (\"hidden\"))) */
      /* #define EXPORT_CLASS __attribute ((visibility (\"default\"))) */
    #else
")

  foreach(pkg ${components})
    set(dll_declspec_content "${dll_declspec_content}
      #define DLL_PUBLIC_${pkg}")
  endforeach()

  set(dll_declspec_content
      "${dll_declspec_content}
      #define DLL_LOCAL
    #endif /* __GNUC__ >=4 */
  #endif /* Win or Linux */

/* When building static libraries, no need to export symbols. ---------------*/
#else /* DRACO_SHARED_LIBS */
")
  foreach(pkg ${components})
    set(dll_declspec_content "${dll_declspec_content}
  #define DLL_PUBLIC_${pkg}")
  endforeach()

  set(dll_declspec_content
      "${dll_declspec_content}
  #define DLL_LOCAL

#endif /* DRACO_SHARED_LIBS */
")

  set(dll_declspec_content
      "${dll_declspec_content}
#endif /* __${safedir}_declspec_h__ */

/*------------------------------------------------------------------------------------------------*/
/* end of dll_declspec.h */
/*------------------------------------------------------------------------------------------------*/
")

  file(WRITE ${PROJECT_BINARY_DIR}/${dir}/dll_declspec.h ${dll_declspec_content})

endfunction()

# -------------------------------------------------------------------------------------------------#
# end config/generate_dll
# -------------------------------------------------------------------------------------------------#
