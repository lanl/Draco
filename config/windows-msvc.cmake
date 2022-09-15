# --------------------------------------------*-cmake-*------------------------------------------- #
# file   windows-msvc.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 June 5
# brief  Establish flags for Visual Studio on Windows. These settings are shared by cl.exe and
#        clang-cl.exe.
# note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

# Useful reference information:
# https://docs.microsoft.com/en-us/cpp/intrinsics/cpuid-cpuidex?view=vs-2017

#
# Sanity Checks
#
if(NOT ${CMAKE_GENERATOR} MATCHES "Visual Studio" AND NOT ${CMAKE_GENERATOR} MATCHES
                                                      "NMake Makefiles")
  message(FATAL_ERROR "config/windows-msvc.cmake must be taught to build for this compiler "
                      "(CMAKE_GENERATOR = ${CMAKE_GENERATOR}). Yell at kt for help on this error.")
endif()

#
# Compiler flag checks
#

# This is required to provide compatibility between MSVC and MinGW generated libraries.
if(DRACO_SHARED_LIBS)
  set(CMAKE_GNUtoMS
      ON
      CACHE BOOL "Compatibility flag for MinGW/MSVC." FORCE)
endif()

# Extra setup (ds++/config.h) for MSVC 1. Allow M_PI to be found via <cmath>
set(_USE_MATH_DEFINES 1)

# Note that in component_macros.cmake, build targets are assigned the property
# WINDOWS_EXPORT_ALL_SYMBOLS=TRUE.  See
# blog.kitware.com/create-dlls-on-windows-without-declspec-using-new-cmake-export-all-feature/ for
# details. However, we still use DLL_PUBLIC macros to control public visibility of global variables.
# See dll_declspec.h generated by ds++'s CMakeLists.txt.

set(MD_or_MT "MD")
set(MD_or_MT_debug "${MD_or_MT}d")
if(DEFINED DEBUG_RUNTIME_EXT AND "${DEBUG_RUNTIME_EXT}" STREQUAL "d")
  set(MD_or_MT_debug "${MD_or_MT}${DEBUG_RUNTIME_EXT} /RTC1")
endif()

set(numproc $ENV{NUMBER_OF_PROCESSORS})
if("${numproc}notfound" STREQUAL "notfound")
  set(numproc 1)
endif()

if(NOT CXX_FLAGS_INITIALIZED)

  # Alternative for per-directory, or per-target specific flags:
  # add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/MP>")

  # Notes on options:
  #
  # * /wd 4251 disable warning #4251: 'identifier' : class 'type' needs to have  dll-interface to be
  #   used by clients of class 'type2'
  # * /wd 5105 After upgrading to VS2019 version 16.8.0 preview 3, warning C5105 is issued from many
  #   system headers like stdio.h. Suppress for now.
  # * /arch:[SSE|SSE2|AVX|AVX2|IA32]
  # * /fsanitize=address
  # * /RTC1 - check for uninitialized variables.
  # * /JMC  - just my code debugging.

  string(APPEND CMAKE_C_FLAGS " /nologo /Gy /fp:precise /DWIN32 /D_WINDOWS /MP /wd4251")
  if(HAVE_HARDWARE_AVX2)
    string(APPEND CMAKE_C_FLAGS " /arch:AVX2")
  endif()
  set(CMAKE_C_FLAGS_DEBUG "/${MD_or_MT_debug} /Od /JMC /Zi /DDEBUG /D_DEBUG")
  set(CMAKE_C_FLAGS_RELEASE "/${MD_or_MT} /O2 /DNDEBUG")
  set(CMAKE_C_FLAGS_MINSIZEREL "/${MD_or_MT} /O1 /DNDEBUG")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "/${MD_or_MT} /O2 /Zi /DDEBUG")

  # Suppress some MSVC warnings about "unsafe" pointer use.
  if(MSVC_VERSION GREATER 1399)
    string(APPEND CMAKE_C_FLAGS
           " /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0")
  endif()

  # If building static libraries, include debugging information in the library.
  if(${DRACO_LIBRARY_TYPE} MATCHES "STATIC")
    string(APPEND CMAKE_C_FLAGS_DEBUG " /Z7")
  endif()

  # /Zc:__cplusplus - enables the __cplusplus preprocessor macro to report an updated value for
  # recent C++ language standards
  string(APPEND CMAKE_CXX_FLAGS " ${CMAKE_C_FLAGS} /EHa /Zc:__cplusplus")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "/${MD_or_MT} /O1 /DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/${MD_or_MT} /O2 /Zi /DDEBUG")

  # Don't warn about missing pdb files.  We won't necessarily have these for TPLs. Link options are
  # applied in component_macros.cmake near calls to add_libarary or add_executable.
  set(DRACO_LINK_OPTIONS "$<$<CONFIG:DEBUG>:/ignore:4099>")

endif()

#
# Toggle compiler flags for optional features
#
if(NOT "${OpenMP_C_FLAGS}x" STREQUAL "x")
  toggle_compiler_flag(OPENMP_FOUND ${OpenMP_C_FLAGS} "C;CXX" "")
endif()
force_compiler_flags_to_cache("C;CXX")

#
# Extra runtime libraries...
#

# Locate a Windows sockets library (required!)
if("${CMAKE_SIZEOF_VOID_P}" STREQUAL 8)
  set(winsock_suffix_dir "x64")
else()
  set(winsock_suffix_dir "x86")
endif()
set(winkitlibdir "$ENV{ProgramFiles\(x86\)}/Windows Kits/10/Lib")
file(GLOB winkitdirs "${winkitlibdir}/*/um/${winsock_suffix_dir}")
foreach(lib ws2_32;wsock32;winsock32;mswsock32)
  if(NOT Lib_win_winsock)
    find_library(
      winsock_lib_${lib}
      NAMES ${lib}
      HINTS C:/Windows/System32;${winkitdirs})
    if(EXISTS "${winsock_lib_${lib}}")
      set(Lib_win_winsock
          "${winsock_lib_${lib}}"
          CACHE FILEPATH "Windows sockets library.")
    endif()
  endif()
endforeach()
unset(winkitdirs)
unset(winkitlibdir)

# winsock is a required dependency.
if(NOT Lib_win_winsock)
  message(FATAL_ERROR "Could not find library wsock32, mswsock32 or ws2_32!")
endif()

# Extra logic to ensure correct winsock is found.
if("${Lib_win_winsock}" MATCHES "um/x86" AND CMAKE_CL_64)
  message(FATAL_ERROR "Found 32-bit winsock (${Lib_win_winsock} but targeting x64 architecture. "
                      "Ensure that cmake is run from the x64 Visual Studio Command Prompt.")
elseif("${Lib_win_winsock}" MATCHES "um/x64" AND NOT CMAKE_CL_64)
  message(FATAL_ERROR "Found 64-bit winsock (${Lib_win_winsock} but targeting x86 architecture. "
                      "Ensure that cmake is run from the x86 Visual Studio Command Prompt.")
endif()

# ------------------------------------------------------------------------------------------------ #
# End windows-msvc.cmake
# ------------------------------------------------------------------------------------------------ #
