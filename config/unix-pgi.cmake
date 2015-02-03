#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-pgi.cmake
# brief  Establish flags for Linux64 - Portland Group C/C++
# note   Copyright (C) 2010-2015 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

set( draco_isPGI 1 CACHE BOOL "Are we using PGI CXX? -> ds++/config.h" )

# NOTE: You may need to set TMPDIR to a location that the compiler can
#       write temporary files to.  Sometimes the default space on HPC
#       is not sufficient.

#
# Sanity Checks
#

if( BUILD_SHARED_LIBS )
  message( FATAL_ERROR "Feature not available - yell at KT." )
endif( BUILD_SHARED_LIBS )

if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.0 )
  message( FATAL_ERROR "PGI must be version 13.0 or later to support C++ features used by this software.")
endif()

#
# Compiler Flags
#
include(platform_checks)
query_openmp_availability()


if( NOT CXX_FLAGS_INITIALIZED )
  set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

  set( CMAKE_C_FLAGS                "-Kieee -Mdaz -pgf90libs" )
  set( CMAKE_C_FLAGS_DEBUG          "-g -O0")
  set( CMAKE_C_FLAGS_RELEASE        "-O3 -DNDEBUG" )
  set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )
  set( CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -gopt" )
  set( CMAKE_CXX_FLAGS              "${CMAKE_C_FLAGS} --no_using_std --no_implicit_include --diag_suppress 940 --diag_suppress 11 --diag_suppress 450 -DNO_PGI_OFFSET" )

  # Extra flags for pgCC-11.2+
  # --nozc_eh    (default for 11.2+) Use low cost exception handling. This
  #              option appears to break our exception handling model resulting
  #              in SEGV.
  # This may be related to PGI bug 1858 (http://www.pgroup.com/support/release_tprs.htm).
  if( "${DBS_CXX_COMPILER_VER_MAJOR}.${DBS_CXX_COMPILER_VER_MINOR}" GREATER 10 )
    if( NOT "${CMAKE_CXX_FLAGS}" MATCHES "--nozc_eh" )
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --nozc_eh" )
    endif()
  endif()

  set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}")
  set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE} -Munroll=c:10 -Mautoinline=levels:10 -Mvect=sse -Mflushz -Mlre")

  # -Mipa=fast,inline
  # -tp x64      Create a PGI Unified Binary which functions correctly
  #              on and is optimized for both Intel and AMD processors.
  # -Mprefetch   Control generation of prefetch instructions to improve
  #              memory performance in compute-intensive loops.

  # -Mnoframe (we use this to debug crashed programs).
  # -Mcache_align (breaks some tests in wedgehog)
  # -Msafeptr (breaks some array operations in MatRA).
  set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_RELEASE}")
  set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -gopt" )

endif()

##---------------------------------------------------------------------------##
# Ensure cache values always match current selection
##---------------------------------------------------------------------------##
set( CMAKE_C_FLAGS                "${CMAKE_C_FLAGS}"                CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )

set( CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS}"                CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )

# Do we add '-ansi -pedantic'?
toggle_compiler_flag( DRACO_ENABLE_C99         "-c99"   "C"   "" )
toggle_compiler_flag( DRACO_ENABLE_STRICT_ANSI "-Xa -A" "CXX" "" )
toggle_compiler_flag( DRACO_ENABLE_CXX11 "--c++11" "CXX" "")
toggle_compiler_flag( OPENMP_FOUND ${OpenMP_C_FLAGS} "C;CXX;EXE_LINKER" "" )

#------------------------------------------------------------------------------#
# End config/unix-pgi.cmake
#------------------------------------------------------------------------------#
