#--------------------------------------------*-cmake-*---------------------------------------------#
# file   config/unix-clang.cmake
# brief  Establish flags for Unix clang
# note   Copyright (C) 2010-2020 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

# Note: In config/compilerEnv.cmake, the build system sets flags for
# 1) the language standard (C++14, C99, etc)
# 2) interprocedural optimization.

#
# Compiler flag checks
#

# Debug flags to consider adding:
# http://clang.llvm.org/docs/UsersManual.html#options-to-control-error-and-warning-messages
# -fdiagnostics-show-hotness
#
# valgrind like options - https://clang.llvm.org/docs/AddressSanitizer.html
#      '-g -fsanitize=address -fno-omit-frame-pointer'
#      must use clang++ for linking
#      suppressions: LSAN_OPTIONS=suppressions=MyLSan.supp
#      human readable: ASAN_SYMBOLIZER_PATH=/usr/local/bin/llvm-symbolizer ./a.out

if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0 AND NOT MSVC )
  message( FATAL_ERROR "Draco requires LLVM clang version >= 6.0." )
endif()

#
# Compiler Flags
#
if( NOT CXX_FLAGS_INITIALIZED )
  set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

  string( APPEND CMAKE_C_FLAGS " -Wcast-align -Wpointer-arith -Wall -Wno-long-long"
    " -Wreserved-id-macro -pedantic" )
  if (NOT ${CMAKE_GENERATOR} MATCHES Xcode AND HAS_MARCH_NATIVE)
    set( CMAKE_C_FLAGS             "${CMAKE_C_FLAGS} -march=native" )
  endif()
  set( CMAKE_C_FLAGS_DEBUG          "-g -fno-inline -O0 -Wextra -DDEBUG")
  set( CMAKE_C_FLAGS_RELEASE        "-O3 -funroll-loops -DNDEBUG" )
  set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )
  set( CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -g -Wextra -funroll-loops" )

  # Suppress warnings about typeid() called with function as an argument. In this case, the function
  # might not be called if the type can be deduced.
  string( APPEND CMAKE_CXX_FLAGS " ${CMAKE_C_FLAGS} -Wno-undefined-var-template"
    " -Wno-potentially-evaluated-expression" )
  if( DEFINED CMAKE_CXX_COMPILER_WRAPPER AND "${CMAKE_CXX_COMPILER_WRAPPER}" STREQUAL "CrayPrgEnv" )
    string(APPEND CMAKE_CXX_FLAGS " -stdlib=libstdc++")
  else()
    string(APPEND CMAKE_CXX_FLAGS " -stdlib=libc++")
  endif()

  set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG} -Woverloaded-virtual")
  # Tried to use -fsanitize=safe-stack but this caused build issues.
  set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}")
  set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_RELEASE}")
  set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" )

endif()

#--------------------------------------------------------------------------------------------------#
# Toggle for OpenMP support
if( OpenMP_C_FLAGS )
  toggle_compiler_flag( OPENMP_FOUND "${OpenMP_C_FLAGS}" "C" "" )
endif()
if( OpenMP_CXX_FLAGS )
  toggle_compiler_flag( OPENMP_FOUND "${OpenMP_CXX_FLAGS}" "CXX" "" )
endif()
# Note: adding openmp option to EXE_LINKER will break MPI detection for gfortran when running with
#       clang++/clang/gfortran.

#--------------------------------------------------------------------------------------------------#
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_CXX_FLAGS)
force_compiler_flags_to_cache()

#--------------------------------------------------------------------------------------------------#
# End config/unix-clang.cmake
#--------------------------------------------------------------------------------------------------#
