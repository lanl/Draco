#--------------------------------------------*-cmake-*---------------------------------------------#
# file   config/apple-clang.cmake
# brief  Establish flags for Apple OSX
# note   Copyright (C) 2016-2020 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

include_guard(GLOBAL)

#
# Compiler Flags
#

# Flags from Draco autoconf build system:

if( NOT CXX_FLAGS_INITIALIZED )
   set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

   string(APPEND CMAKE_C_FLAGS " -Wcast-align -Wpointer-arith -Wall -Wno-long-long -pedantic" )
   if (NOT ${CMAKE_GENERATOR} MATCHES Xcode AND HAS_MARCH_NATIVE)
     string(APPEND CMAKE_C_FLAGS " -march=native" )
   endif()
   set( CMAKE_C_FLAGS_DEBUG          "-g -fno-inline -O0 -Wextra -DDEBUG")
   set( CMAKE_C_FLAGS_RELEASE        "-O3 -funroll-loops -DNDEBUG" )
   set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )
   set( CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -g -Wextra -funroll-loops" )

   string(APPEND CMAKE_CXX_FLAGS " ${CMAKE_C_FLAGS} -stdlib=libc++ -std=c++11" )
   set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG} -Woverloaded-virtual")
   set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}")
   set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_RELEASE}")
   set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" )

endif()

#--------------------------------------------------------------------------------------------------#
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_C_FLAGS)
deduplicate_flags(CMAKE_CXX_FLAGS)
force_compiler_flags_to_cache("C;CXX")

#--------------------------------------------------------------------------------------------------#
# End config/apple-clang.cmake
#--------------------------------------------------------------------------------------------------#
