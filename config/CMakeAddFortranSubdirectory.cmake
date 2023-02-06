#.rst:
# CMakeAddFortranSubdirectory
# ---------------------------
#
# Use a version of gfortran that is not available from within the current project.  For example, use
# MinGW gfortran from Visual Studio if a Fortran compiler is not found, or use GNU gfortran from a
# XCode/clang build project.
#
# The 'add_fortran_subdirectory' function adds a subdirectory to a project that contains a Fortran
# only sub-project.  The module will check the current compiler and see if it can support Fortran.
# If no Fortran compiler is found and the compiler is MSVC or if the Generator is XCode, then this
# module will try to find a gfortran compiler in local environment (e.g.: MinGW gfortran).  It will
# then use an external project to build with alternate (MinGW/Unix) tools.  It will also create
# imported targets for the libraries created.
#
# For visual studio, this will only work if the Fortran code is built into a dll, so
# BUILD_SHARED_LIBS is turned on in the project. In addition the CMAKE_GNUtoMS option is set to on,
# so that the MS .lib files are created.
#
# Usage is as follows:
#
# ::
#
#   cmake_add_fortran_subdirectory(
#    <subdir>                # name of subdirectory
#    PROJECT <cafs_project_name>  # project name in subdirectories's top CMakeLists.txt
#                            # recommendation: use the same project name as listed in
#                            # <subdir>/CMakeLists.txt
#    ARCHIVE_DIR <dir>       # directory where project places .lib files
#    RUNTIME_DIR <dir>       # directory where project places .dll files
#    LIBRARIES <lib>...      # names of library targets to import
#    TARGET_NAMES <string>...# target names assigned to the libraries listed above available in the
#                            # primary project.
#    LINK_LIBRARIES          # link interface libraries for LIBRARIES
#     [LINK_LIBS <lib> <dep>...]...
#    DEPENDS                 # Register dependencies external for this AFSD project.
#    CMAKE_COMMAND_LINE ...  # extra command line flags to pass to cmake
#    NO_EXTERNAL_INSTALL     # skip installation of external project
#    )
#
# Relative paths in ARCHIVE_DIR and RUNTIME_DIR are interpreted with respect to the build directory
# corresponding to the source directory in which the function is invoked.

# ================================================================================================ #
# This is a heavily modified version of CMakeAddFortranSubdirectory.cmake that is distributed with
# CMake

# Copyright (C) 2014-2023 Triad National Security, LLC., All rights reserved.

# cmake-lint: disable=R0912,R0915,W0106

set(_CAFS_CURRENT_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})
include(CheckLanguage)
include(ExternalProject)

# ------------------------------------------------------------------------------------------------ #
# Find gfortran and check/setup x86/x64 information.
# ------------------------------------------------------------------------------------------------ #
function(_setup_cafs_config_and_build source_dir build_dir)

  # Try to find a Fortran compiler (use MinGW gfortran for MSVC).
  find_program(
    CAFS_Fortran_COMPILER
    NAMES ${CAFS_Fortran_COMPILER} gfortran
    PATHS c:/MinGW/bin c:/msys64/mingw64/bin)
  if(NOT EXISTS ${CAFS_Fortran_COMPILER})
    message(
      FATAL_ERROR
        "A Fortran compiler was not found.  Please set CAFS_Fortran_COMPILER to the full path of a "
        "working Fortran compiler. For Windows platforms, you need to install MinGW with the "
        "gfortran option.")
  endif()

  # Validate flavor/architecture of specified gfortran
  if(MSVC)
    # MinGW gfortran under MSVS.
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
      set(_cafs_fortran_target_arch "Target:.*64.*mingw")
    else()
      set(_cafs_fortran_target_arch "Target:.*mingw32")
    endif()
  elseif(APPLE)
    # GNU gfortran under XCode.
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
      set(_cafs_fortran_target_arch "Target:.*64-apple*")
    else()
      set(_cafs_fortran_target_arch "Target:.*86-apple*")
    endif()
  else()
    # GNU gfortran with Ninja generator or clang CXX compiler.
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
      set(_cafs_fortran_target_arch "Target: x86_64*|ppc64")
    else()
      set(_cafs_fortran_target_arch "Target:.*86*")
    endif()
  endif() # MSVC
  execute_process(
    COMMAND "${CAFS_Fortran_COMPILER}" -v
    ERROR_VARIABLE out
    ERROR_STRIP_TRAILING_WHITESPACE)
  if(NOT "${out}" MATCHES "${_cafs_fortran_target_arch}")
    string(REPLACE "\n" "\n  " out "  ${out}")
    message(
      FATAL_ERROR
        "CAFS_Fortran_COMPILER is set to\n${CAFS_Fortran_COMPILER}\n, which is not a valid Fortran "
        "compiler for this architecture.  The output from '${CAFS_Fortran_COMPILER} -v' does not "
        "match '${_cafs_fortran_target_arch}':\n${out}\nSet CAFS_Fortran_COMPILER to a compatible "
        "Fortran compiler for this architecture.")
  endif()

  # Configure scripts to run Fortran tools with the proper PATH.
  get_filename_component(CAFS_Fortran_COMPILER_PATH ${CAFS_Fortran_COMPILER} PATH)
  file(TO_NATIVE_PATH "${CAFS_Fortran_COMPILER_PATH}" CAFS_Fortran_COMPILER_PATH)
  string(REPLACE "\\" "\\\\" CAFS_Fortran_COMPILER_PATH "${CAFS_Fortran_COMPILER_PATH}")
  # Generator type
  if(MSVC)
    set(GENERATOR_TYPE "-GMinGW Makefiles")
    set(CAFS_GNUtoMS "-DCMAKE_GNUtoMS=ON")
  else() # Unix Makefiles, Xcode or Ninja.
    set(GENERATOR_TYPE "-GUnix Makefiles")
  endif()

  # Generate the config_cafs_proj.cmake command file:
  if(ARGS_VERBOSE)
    message(
      "
    configure_file(
    ${_CAFS_CURRENT_SOURCE_DIR}/CMakeAddFortranSubdirectory/config_cafs_proj.cmake.in
    ${build_dir}/config_cafs_proj.cmake
    @ONLY)
  ")
  endif()
  configure_file(${_CAFS_CURRENT_SOURCE_DIR}/CMakeAddFortranSubdirectory/config_cafs_proj.cmake.in
                 ${build_dir}/config_cafs_proj.cmake @ONLY)

  # Generate the build_cafs_proj.cmake command file:
  include(ProcessorCount)
  processorcount(NumProcs)
  set(build_command_args "--build . -j ${NumProcs}")
  if(NOT ARGS_NO_EXTERNAL_INSTALL)
    string(APPEND build_command_args " --target install")
  endif()
  set(build_cafs_proj_command
      "# This file generated by CMakeAddFortranSubdirectory.cmake
#------------------------------------------------------------------------
set(ENV{PATH} \"${CAFS_Fortran_COMPILER_PATH}\;\$ENV{PATH}\")
set(VERBOSE ${ARGS_VERBOSE})
if(VERBOSE)
    message(\"${CMAKE_COMMAND} ${build_command_args}\")
endif()
execute_process(COMMAND \"${CMAKE_COMMAND}\" ${build_command_args})

# end build_cafs_proj.cmake
#------------------------------------------------------------------------
")
  if(ARGS_VERBOSE)
    message("Generating ${build_dir}/build_cafs_proj.cmake")
  endif()
  file(WRITE "${build_dir}/build_cafs_proj.cmake" ${build_cafs_proj_command})

endfunction()

# ------------------------------------------------------------------------------------------------ #
# _add_fortran_library_link_interface
# ------------------------------------------------------------------------------------------------ #
function(_add_fortran_library_link_interface library depend_library)
  set_target_properties(
    ${library}
    PROPERTIES
    IMPORTED_LINK_INTERFACE_LIBRARIES_NOCONFIG "${depend_library}")
  if(ARGS_VERBOSE)
    message(
      "
  set_target_properties(
    ${library}
    PROPERTIES
    IMPORTED_LINK_INTERFACE_LIBRARIES_NOCONFIG \"${depend_library}\")
")
  endif()
endfunction()

# ------------------------------------------------------------------------------------------------ #
# This is the main function.  This generates the required external_project pieces that will be run
# under a different generator (MinGW Makefiles).
# ------------------------------------------------------------------------------------------------ #
function(cmake_add_fortran_subdirectory subdir)

  # Parse arguments to function
  set(options NO_EXTERNAL_INSTALL VERBOSE)
  set(oneValueArgs PROJECT ARCHIVE_DIR RUNTIME_DIR)
  set(multiValueArgs LIBRARIES TARGET_NAMES LINK_LIBRARIES DEPENDS CMAKE_COMMAND_LINE)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(NOT ARGS_NO_EXTERNAL_INSTALL AND ARGS_VERBOSE)
    message(
      "
-- The external_project ${ARGS_PROJECT} will be installed to the location specified by
   CMAKE_INSTALL_PREFIX. This install location should be set via the values provided by
   cmake_add_fortran_subdirectory's CMAKE_COMMAND_LINE parameter.
   ")
  endif()

  # If the current generator/system already supports Fortran, then simply add the requested
  # directory to the project.
  if(ARGS_VERBOSE)
    if(MSVC)
      message("MSVC = TRUE")
    endif()
    message(
      "
    _LANGUAGES_            = ${_LANGUAGES_}
    CMAKE_Fortran_COMPILER = ${CMAKE_Fortran_COMPILER}
    ")
  endif()
  if(_LANGUAGES_ MATCHES Fortran OR (MSVC AND "${CMAKE_Fortran_COMPILER}" MATCHES ifort))
    add_subdirectory(${subdir})
    return()
  endif()

  # Setup external projects to build with alternate Fortran:
  set(source_dir "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}")
  set(cafs_project_name "${ARGS_PROJECT}")
  set(library_dir "${ARGS_ARCHIVE_DIR}")
  set(binary_dir "${ARGS_RUNTIME_DIR}")
  set(libraries ${ARGS_LIBRARIES})
  set(target_names "${ARGS_TARGET_NAMES}")
  list(LENGTH libraries numlibs)
  list(LENGTH target_names numtgtnames)
  if(ARGS_VERBOSE)
    message(
      "
    Preparing CAFS external project:
    - source_dir        = \"${CMAKE_CURRENT_SOURCE_DIR}/${subdir}\"
    - cafs_project_name = \"${ARGS_PROJECT}\"
    - library_dir       = \"${ARGS_ARCHIVE_DIR}\"
    - binary_dir        = \"${ARGS_RUNTIME_DIR}\"
    - libraries         = \"${ARGS_LIBRARIES}\"
    - target_names      = \"${ARGS_TARGET_NAMES}\"
    - numlibs           = ${numlibs}
    - numtgtnames       = ${numtgtnames}
    ")
  endif()
  if(${numtgtnames} STREQUAL 0)
    set(target_names ${libraries})
    set(numtgtnames ${numlibs})
  endif()
  if(NOT ${numlibs} STREQUAL ${numtgtnames})
    message(FATAL_ERROR "If TARGET_NAMES are provided, you must provide an equal number of entries"
                        " for both TARGET_NAMES and LIBRARIES.")
  endif()
  # use the same directory that add_subdirectory would have used
  set(build_dir "${CMAKE_CURRENT_BINARY_DIR}/${subdir}")
  foreach(dir_var library_dir binary_dir)
    if(NOT IS_ABSOLUTE "${${dir_var}}")
      get_filename_component(${dir_var} "${CMAKE_CURRENT_BINARY_DIR}/${${dir_var}}" ABSOLUTE)
    endif()
  endforeach()
  # create build and configure wrapper scripts
  _setup_cafs_config_and_build("${source_dir}" "${build_dir}")
  # If the current build tool has multiple configurations, use the generator expression $<CONFIG> to
  # drive the build type for the Fortran subproject. Otherwise, force the Fortran subproject to use
  # the same build type as the main project.
  if(CMAKE_CONFIGURATION_TYPES)
    set(ep_build_type "$<CONFIG>")
  else()
    set(ep_build_type "${CMAKE_BUILD_TYPE}")
  endif()
  # create the external project
  if(ARGS_VERBOSE)
    message(
      "
    externalproject_add(${cafs_project_name}_build
    DEPENDS           ${ARGS_DEPENDS}
    SOURCE_DIR        ${source_dir}
    BINARY_DIR        ${build_dir}
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=${ep_build_type}
                      -P ${build_dir}/config_cafs_proj.cmake
    BUILD_COMMAND     ${CMAKE_COMMAND} -P ${build_dir}/build_cafs_proj.cmake
    BUILD_ALWAYS 1
    INSTALL_COMMAND   \"\"
    )
  ")
  endif()
  ExternalProject_Add(
    ${cafs_project_name}_build
    DEPENDS ${ARGS_DEPENDS}
    SOURCE_DIR ${source_dir}
    BINARY_DIR ${build_dir}
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=${ep_build_type} -P
                      ${build_dir}/config_cafs_proj.cmake
    BUILD_COMMAND ${CMAKE_COMMAND} -P ${build_dir}/build_cafs_proj.cmake
    BUILD_ALWAYS 1
    INSTALL_COMMAND "")
  # create imported targets for all libraries
  set(idx 0)
  foreach(lib ${libraries})
    list(GET target_names ${idx} tgt)
    if(ARGS_VERBOSE)
      message("    add_library(${tgt} SHARED IMPORTED GLOBAL)")
    endif()
    add_library(${tgt} SHARED IMPORTED GLOBAL)
    if(CMAKE_RUNTIME_OUTPUT_DIRECTORY)
      set(crod "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
      if(ARGS_VERBOSE)
        message(
          "    set_target_properties(
        ${tgt}
        PROPERTIES
          IMPORTED_LOCATION \"${tmp}/${ep_build_type}/lib${lib}${CMAKE_SHARED_LIBRARY_SUFFIX}\"
          IMPORTED_LINK_INTERFACE_LANGUAGES \"Fortran\"
          IMPORTED_LINK_INTERFACE_LIBRARIES \"${ARGS_DEPENDS}\"
        )    ")
      endif()
      set_target_properties(
        ${tgt} PROPERTIES
        IMPORTED_LOCATION "${crod}/${ep_build_type}/lib${lib}${CMAKE_SHARED_LIBRARY_SUFFIX}"
        IMPORTED_LINK_INTERFACE_LANGUAGES "Fortran"
        IMPORTED_LINK_INTERFACE_LIBRARIES "${ARGS_DEPENDS}")
      unset(crod)
    else()
      if(ARGS_VERBOSE)
        message(
          "    set_target_properties(
        ${tgt}
        PROPERTIES
          IMPORTED_LOCATION \"${binary_dir}/lib${lib}${CMAKE_SHARED_LIBRARY_SUFFIX}\"
          IMPORTED_LINK_INTERFACE_LANGUAGES \"Fortran\"
          IMPORTED_LINK_INTERFACE_LIBRARIES \"${ARGS_DEPENDS}\"
        )")
      endif()
      set_target_properties(
        ${tgt} PROPERTIES
        IMPORTED_LOCATION "${binary_dir}/lib${lib}${CMAKE_SHARED_LIBRARY_SUFFIX}"
        IMPORTED_LINK_INTERFACE_LANGUAGES "Fortran"
        IMPORTED_LINK_INTERFACE_LIBRARIES "${ARGS_DEPENDS}")
    endif()
    if(WIN32)
      if(ARGS_VERBOSE)
        message("    set_target_properties(${tgt} PROPERTIES
        IMPORTED_IMPLIB
          \"${library_dir}/lib${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}\" )")
      endif()
      set_target_properties(
        ${tgt} PROPERTIES IMPORTED_IMPLIB "${library_dir}/lib${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif()
    # [2015-01-29 KT/Wollaber: We don't understand why this is needed, but adding
    # IMPORTED_LOCATION_DEBUG to the target_properties fixes a missing RPATH problem for Xcode
    # builds.  Possibly, this is related to the fact that the Fortran project is always built in
    # Debug mode.
    if(APPLE)
      set_target_properties(
        ${tgt} PROPERTIES
        IMPORTED_LOCATION_DEBUG "${binary_dir}/lib${lib}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    endif()
    add_dependencies(${tgt} ${cafs_project_name}_build)

    # The Ninja Generator appears to want to find the imported library
    # ${binary_dir}/lib${lib}${CMAKE_SHARED_LIBRARY_SUFFIX or a rule to generate this target before
    # it runs any build commands.  Since this library will not exist until the external project is
    # built, we need to trick Ninja by creating a place-holder file to satisfy Ninja's dependency
    # checker.  This library will be overwritten during the actual build.
    if(${CMAKE_GENERATOR} MATCHES Ninja)
      # artificially create some targets to help Ninja resolve dependencies.
      execute_process(COMMAND ${CMAKE_COMMAND} -E touch
                              "${binary_dir}/lib${lib}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    endif()

    if(ARGS_VERBOSE)
      message(
        "
cmake_add_fortran_subdirectory
   Project    : ${cafs_project_name}
   Directory  : ${source_dir}
   Target name: ${tgt}
   Library    : ${binary_dir}/lib${lib}${CMAKE_SHARED_LIBRARY_SUFFIX}
   Target deps: ${cafs_project_name}_build --> ${ARGS_DEPENDS}
   Extra args : ${ARGS_CMAKE_COMMAND_LINE}
      ")
      include(print_target_properties)
      print_targets_properties(${tgt})
    endif()
    math(EXPR idx "${idx} + 1")
  endforeach()

  # now setup link libraries for targets
  set(start FALSE)
  set(target)
  foreach(lib ${ARGS_LINK_LIBRARIES})
    if("${lib}" STREQUAL "LINK_LIBS")
      set(start TRUE)
    else()
      if(start)
        if(DEFINED target)
          # process current target and target_libs
          _add_fortran_library_link_interface(${target} "${target_libs}")
          # zero out target and target_libs
          set(target)
          set(target_libs)
        endif()
        # save the current target and set start to FALSE
        set(target ${lib})
        set(start FALSE)
      else()
        # append the lib to target_libs
        list(APPEND target_libs "${lib}")
      endif()
    endif()
  endforeach()
  # process anything that is left in target and target_libs
  if(DEFINED target)
    _add_fortran_library_link_interface(${target} "${target_libs}")
  endif()

  # If we are installing this target, then create a string that can be saved to the installed
  # <project>-config.cmake file for import by clients. This information will be saved in the
  # variable CAFS_EXPORT_TARGET_PROPERTIES
  if(NOT ARGS_NO_EXTERNAL_INSTALL AND CMAKE_CONFIGURATION_TYPES)

    string(
      APPEND
      CAFS_EXPORT_DEFINE_IMPORT_PREFIX
      "
# Compute the installation prefix relative to this file.
# (generated by CMakeAddFortranSubdirectory.cmake)
if( NOT DEFINED _IMPORT_PREFIX )
  get_filename_component(_IMPORT_PREFIX \"\${CMAKE_CURRENT_LIST_FILE}\" PATH)
  get_filename_component(_IMPORT_PREFIX \"\${_IMPORT_PREFIX}\" PATH)
  if(_IMPORT_PREFIX STREQUAL \"/\")
    set(_IMPORT_PREFIX \"\")
  endif()
endif()
    ")

    string(
      APPEND
      CAFS_EXPORT_TARGET_PROPERTIES
      "
# Create imported target ${ARGS_TARGET_NAMES}
# (generated by CMakeAddFortranSubdirectory.cmake)
add_library(${ARGS_TARGET_NAMES} SHARED IMPORTED)
foreach( build_type ${CMAKE_CONFIGURATION_TYPES} )
  string( TOUPPER \${build_type} build_type_allcaps )
  if( EXISTS \"\${_IMPORT_PREFIX}/\${build_type}/lib/lib${ARGS_LIBRARIES}.lib\" AND
      EXISTS \"\${_IMPORT_PREFIX}/\${build_type}/bin/lib${ARGS_LIBRARIES}.dll\" AND
      EXISTS \"\${_IMPORT_PREFIX}/\${build_type}/include/${subdir}\" )
    set_property(TARGET ${ARGS_TARGET_NAMES} APPEND PROPERTY IMPORTED_CONFIGURATIONS \${build_type})
    set_target_properties(${ARGS_TARGET_NAMES} PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES \"\${_IMPORT_PREFIX}/\${build_type}/include/${subdir}\"
      INTERFACE_LINK_LIBRARIES \"${ARGS_DEPENDS}\"
      IMPORTED_IMPLIB_\${build_type_allcaps}
        \"\${_IMPORT_PREFIX}/\${build_type}/lib/lib${ARGS_LIBRARIES}.lib\"
      IMPORTED_LOCATION_\${build_type_allcaps}
        \"\${_IMPORT_PREFIX}/\${build_type}/bin/lib${ARGS_LIBRARIES}.dll\" )
  endif()
endforeach()
    ")
    set(CAFS_EXPORT_TARGET_PROPERTIES
        "${CAFS_EXPORT_TARGET_PROPERTIES}"
        PARENT_SCOPE)
    set(CAFS_EXPORT_DEFINE_IMPORT_PREFIX
        "${CAFS_EXPORT_DEFINE_IMPORT_PREFIX}"
        PARENT_SCOPE)
  endif()

endfunction()

# ------------------------------------------------------------------------------------------------ #
# When generating Visual Studio projects and using MSYS/MinGW gfortran for CAFS -based Fortran
# subprojects, you must link to an MSYS/MinGW style '.a' library instead of msmpi.lib
# (implementation portion of msmpi.dll). This macro provides this fix (call this from the
# CMakeLists.txt found in the CAFS subproject.
#
# Don't link to the C++ MS-MPI library when compiling with MinGW gfortran. Instead, link to
# libmsmpi.a that was created via gendef.exe and dlltool.exe from msmpi.dll.  Refs:
#
# * https://github.com/KineticTheory/Linux-HPC-Env/wiki/Setup-Win32-development-environment,
# * http://www.geuz.org/pipermail/getdp/2012/001519.html
#
# There are also issues with MS-MPI's mpif.h when using gfortran's '-frange-check' compiler flag. If
# this compiler option is enabled, remove it.
# ------------------------------------------------------------------------------------------------ #
function(cafs_fix_mpi_library)

  set(verbose false)

  # MS-MPI and gfortran do not play nice together...
  if(WIN32 AND "${DRACO_C4}" STREQUAL "MPI")
    if(verbose)
      message("CAFS: MPI_Fortran_LIBRARIES= ${MPI_Fortran_LIBRARIES}")
    endif()
    if(NOT MPI_Fortran_LIBRARIES OR "${MPI_Fortran_LIBRARIES}" MATCHES "msmpi.lib")
      # msmpi.lib should be located one of:
      #
      # * C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\[x86|x64]
      # * %vcpkg%/installed/x64-windows/[debug]/lib
      if(MPI_msmpi_LIBRARY)
        # FindMPI.cmake should have set $MPI_msmpi_LIBRARY.
        get_filename_component(MSMPI_SDK_DIR ${MPI_msmpi_LIBRARY} DIRECTORY)
        if(verbose)
          message("CAFS: MSMPI_SDK_DIR = ${MSMPI_SDK_DIR}")
        endif()
        find_library(
          MPI_gfortran_LIBRARIES
          NAMES "libmsmpi.a"
          HINTS ${MSMPI_SDK_DIR})
        unset(MSMPI_SDK_DIR)
      else()
        find_file(MPI_gfortran_LIBRARIES NAMES "libmsmpi.a")
      endif()
    elseif("${MPI_Fortran_LIBRARIES}" MATCHES "libmsmpi.a")
      foreach(lib ${MPI_Fortran_LIBRARIES})
        if("${lib}" MATCHES "libmsmpi.a")
          set(MPI_gfortran_LIBRARIES "${lib}")
        endif()
      endforeach()
    endif()

    # Sanity check
    if(MPI_gfortran_LIBRARIES)
      set(MPI_gfortran_LIBRARIES
          ${MPI_gfortran_LIBRARIES}
          CACHE FILEPATH "msmpi for gfortran" FORCE)
    else()
      message(
        FATAL_ERROR
          "Unable to find libmsmpi.a. This library must be created from msmpi.dll and saved as a "
          "MinGW library. For more help see "
          "https://github.com/KineticTheory/Linux-HPC-Env/wiki/Setup-Win32-development-environment")
    endif()

    # 1. Create a new imported target to represent MPI::CAFS
    add_library(MPI::MPI_cafs SHARED IMPORTED)
    set_target_properties(
      MPI::MPI_cafs PROPERTIES
      IMPORTED_CONFIGURATIONS RELEASE
      INTERFACE_INCLUDE_DIRECTORIES "${MPI_C_HEADER_DIR};${MPI_mpifptr_INCLUDE_DIR}"
      IMPORTED_IMPLIB_RELEASE "${MPI_gfortran_LIBRARIES}"
      IMPORTED_LOCATION_RELEASE "${MPI_gfortran_LIBRARIES}")
    # 1. Strip MPI deps from Lib_c4
    set_target_properties(Lib_c4 PROPERTIES INTERFACE_LINK_LIBRARIES "Lib_dsxx;MPI::MPI_cafs")

    # Remove '-fno-range-check' gfortran compiler flag (remove it from Debug flags)
    if("${CMAKE_Fortran_FLAGS_DEBUG}" MATCHES "frange-check")
      string(REPLACE "range-check" "no-range-check" CMAKE_Fortran_FLAGS_DEBUG
                     ${CMAKE_Fortran_FLAGS_DEBUG})
      set(CMAKE_Fortran_FLAGS_DEBUG
          "${CMAKE_Fortran_FLAGS_DEBUG}"
          PARENT_SCOPE)
    endif()
    string(APPEND CMAKE_Fortran_FLAGS " -fno-range-check")
    set(CMAKE_Fortran_FLAGS
        "${CMAKE_Fortran_FLAGS}"
        PARENT_SCOPE)
    if(verbose)
      message("CAFS: MPI_gfortran_LIBRARIES= ${MPI_gfortran_LIBRARIES}"
              "\n      CMAKE_Fortran_FLAGS_DEBUG = ${CMAKE_Fortran_FLAGS_DEBUG}")
    endif()
  endif()

endfunction(cafs_fix_mpi_library)

# ------------------------------------------------------------------------------------------------ #
# Create imported libraries owned by the Visual Studio project to be used in the CAFS subproject.
# ------------------------------------------------------------------------------------------------ #
function(cafs_create_imported_targets targetName libName targetPath linkLang)

  set(verbose_target_name "Lib_foo")
  if(targetName STREQUAL verbose_target_name)
    message(
      "
    targetName = ${targetName}
    libName    = ${libName}
    targetPath = ${targetPath}
    linkLang   = ${linkLang}")
  endif()

  get_filename_component(pkgloc "${targetPath}" ABSOLUTE)
  if(targetName STREQUAL verbose_target_name)
    message(
      "
    find_library(
      lib
      NAMES ${libName}
      PATHS ${pkgloc}
      PATH_SUFFIXES Release Debug )
  ")
  endif()
  find_library(
    lib
    NAMES ${libName}
    PATHS ${pkgloc}
    PATH_SUFFIXES Release Debug)
  get_filename_component(libloc ${lib} DIRECTORY)
  if(targetName STREQUAL verbose_target_name)
    message("      libloc = ${libloc}")
  endif()
  if("${libloc}x" STREQUAL "x")
    message(
      FATAL_ERROR
        "cafs_create_imported_targets :: Did not find library ${lib} at location ${pkgloc} when "
        "trying to create import target ${targetName}.\n" "libName    = ${libName}\n"
        "targetPath = ${targetPath}")
  endif()

  # Debug case?
  find_library(
    lib_debug
    NAMES ${libName}
    PATHS ${pkgloc}/Debug)
  get_filename_component(libloc_debug ${lib_debug} DIRECTORY)

  #
  # Generate the imported library target and set properties...
  #
  if(DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(dll_loc "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
    set(dll_loc_debug "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
  else()
    set(dll_loc "${libloc}")
    set(dll_loc_debug "${libloc_debug}")
  endif()

  add_library(${targetName} SHARED IMPORTED GLOBAL)
  set_target_properties(
    ${targetName} PROPERTIES
    IMPORTED_LOCATION
      "${dll_loc}/${CMAKE_SHARED_LIBRARY_PREFIX}/${libName}${CMAKE_SHARED_LIBRARY_SUFFIX}"
    IMPORTED_LINK_INTERFACE_LANGUAGES ${linkLang})
  if(targetName STREQUAL verbose_target_name)
    message(
      "
  add_library( ${targetName} SHARED IMPORTED GLOBAL )
  set_target_properties(
    ${targetName} PROPERTIES
    IMPORTED_LOCATION
      ${dll_loc}/${CMAKE_SHARED_LIBRARY_PREFIX}/${libName}${CMAKE_SHARED_LIBRARY_SUFFIX}
    IMPORTED_LINK_INTERFACE_LANGUAGES ${linkLang}")
  endif()
  if(lib_debug)
    set_target_properties(
      ${targetName} PROPERTIES
      IMPORTED_LOCATION_DEBUG
        "${dll_loc_debug}/${CMAKE_SHARED_LIBRARY_PREFIX}${libName}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    if(targetName STREQUAL verbose_target_name)
      message("    IMPORTED_LOCATION_DEBUG
      ${dll_loc_debug}/${CMAKE_SHARED_LIBRARY_PREFIX}${libName}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    endif()
  endif()

  # Fortran mod files
  if(linkLang STREQUAL Fortran)
    set_target_properties(${targetName} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${targetPath})
    if(targetName STREQUAL verbose_target_name)
      message("    INTERFACE_INCLUDE_DIRECTORIES ${targetPath}")
    endif()
  endif()

  # platform specific properties
  if(APPLE)
    set_target_properties(${targetName} PROPERTIES MACOSX_RPATH TRUE)
  elseif(WIN32)
    if(CMAKE_GNUtoMS)
      set(CMAKE_IMPORT_LIBRARY_PREFIX "")
      set(CMAKE_IMPORT_LIBRARY_SUFFIX ".lib")
    endif()
    set_target_properties(
      ${targetName} PROPERTIES
      IMPORTED_IMPLIB
        "${libloc}/${CMAKE_IMPORT_LIBRARY_PREFIX}${libName}${CMAKE_IMPORT_LIBRARY_SUFFIX}")
    if(targetName STREQUAL verbose_target_name)
      message("    IMPORTED_IMPLIB
      ${libloc}/${CMAKE_IMPORT_LIBRARY_PREFIX}${libName}${CMAKE_IMPORT_LIBRARY_SUFFIX}")
    endif()
    if(lib_debug)
      set_target_properties(
        ${targetName} PROPERTIES
        IMPORTED_IMPLIB_DEBUG
          "${libloc_debug}/${CMAKE_IMPORT_LIBRARY_PREFIX}${libName}${CMAKE_IMPORT_LIBRARY_SUFFIX}")
      if(targetName STREQUAL verbose_target_name)
        message("    IMPORTED_IMPLIB_DEBUG
      ${libloc_debug}/${CMAKE_IMPORT_LIBRARY_PREFIX}${libName}${CMAKE_IMPORT_LIBRARY_SUFFIX}")
      endif()
    endif()
  endif()
  if(targetName STREQUAL verbose_target_name)
    message("  )")
  endif()
  unset(lib CACHE)
  unset(lib_debug CACHE)
endfunction()

# ------------------------------------------------------------------------------------------------ #
# Generate a set of default variables for the CAFS cmake command line. This macro sets these cmake
# variables:
#
# * build_system_state
#
# Optional arguments:
#
# * PROJECT "Draco"
#
# ------------------------------------------------------------------------------------------------ #
function(init_build_system_state)

  # If the current generator/system already supports Fortran, then simply add the requested
  # directory to the project.
  if(_LANGUAGES_ MATCHES Fortran OR (MSVC AND "${CMAKE_Fortran_COMPILER}" MATCHES ifort))
    return()
  endif()

  # Parse arguments to function
  set(options)
  set(oneValueArgs PROJECT)
  set(multiValueArgs)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(NOT ARGS_PROJECT)
    set(ARGS_PROJECT "Draco")
    set(draco_DIR "${Draco_SOURCE_DIR}/config")
  endif()

  # The alternate build system (Makefiles if we are Apple/OSX or Linux/Ninja) will need some of the
  # current build system parameters:
  set(build_system_state "-DDRACO_C4=${DRACO_C4}" "-DDRACO_LIBRARY_TYPE=${DRACO_LIBRARY_TYPE}"
                         "-Ddraco_DIR=${draco_DIR}")
  if(${DRACO_C4} MATCHES "MPI")
    list(APPEND build_system_state "-DMPI_C_LIBRARIES=${MPI_C_LIBRARIES}"
         "-DMPI_C_INCLUDE_DIRS=${MPI_C_INCLUDE_DIRS}")
  endif()
  if(WIN32)
    # For Win32 builds, DLL and applications are built in the directory specified by
    # CMAKE_RUNTIME_OUTPUT_DIRECTORY.
    list(APPEND build_system_state
         "-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=${${ARGS_PROJECT}_BINARY_DIR}/\${CMAKE_BUILD_TYPE}")
    if(CMAKE_TOOLCHAIN_FILE)
      list(APPEND build_system_state "-DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}")
    endif()
  else()
    list(APPEND build_system_state "-DHAVE_GPU=${HAVE_GPU}")
  endif()
  set(build_system_state
      "${build_system_state}"
      PARENT_SCOPE)
  set(draco_DIR
      "${draco_DIR}"
      PARENT_SCOPE)

endfunction()

# ------------------------------------------------------------------------------------------------ #
# Capture boilerplate setup for Fortran-only directories built with CAFS
#
# 1. Ensure that draco/config is listed in CMAKE_MODULE_PATH
# 2. Include basic build system setup routines that define helper macros used by the main Draco bld
#    system (do platform_checks, compilerEnv, vendors, etc.)
# 3. Define draco_BINARY_DIR and create local import targets for common dependencies like Lib_dsxx
#    and Lib_c4.
# 4. Include extra directories to find header files.
#
# ------------------------------------------------------------------------------------------------ #
macro(CAFS_Fortran_dir_boilerplate_setup)

  # Parse arguments to function
  set(options)
  set(oneValueArgs PROJECT)
  set(multiValueArgs)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(NOT ARGS_PROJECT)
    set(ARGS_PROJECT "Draco")
  endif()

  if(NOT TARGET Lib_dsxx)

    # Build system configuration files are located here.
    if(NOT DEFINED draco_DIR OR NOT EXISTS ${draco_DIR})
      message(FATAL_ERROR "can't find draco/config directory at draco_DIR = \"${draco_DIR}\"")
    endif()

    if(NOT ARGS_PROJECT STREQUAL "Draco")
      set(CROD "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
      find_package(draco REQUIRED CONFIG)
      dbs_basic_setup()
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CROD}")
      unset(CROD)
    else()
      get_filename_component(draco_BINARY_DIR ${PROJECT_BINARY_DIR}/../../.. ABSOLUTE)
      cafs_create_imported_targets(Lib_dsxx "rtt_ds++" "${draco_BINARY_DIR}/src/ds++" CXX)
      cafs_create_imported_targets(Lib_c4 "rtt_c4" "${draco_BINARY_DIR}/src/c4" CXX)

      # If we get here, we also need to use the Draco scripts to setup compiler flags and MPI
      # options.
      include(buildEnv)
      dbssetdefaults()
      # On Win32, set the default top level output directory:
      if(WIN32)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
            ${draco_BINARY_DIR}/${CMAKE_BUILD_TYPE}
            CACHE PATH "Build runtime objects at this location.")
      endif()
      include(platform_checks)
      include(compilerEnv)
      dbssetupfortran()
      include(component_macros)
      include(vendor_libraries)
    endif()

    # include_directories( "${cafs_Draco_TPL_INCLUDE_DIRS}" )
  endif(NOT TARGET Lib_dsxx)
endmacro(CAFS_Fortran_dir_boilerplate_setup)

# ------------------------------------------------------------------------------------------------ #
# Capture MPI-specific boilerplate setup for Fortran-only directories built with CAFS
#
# 1. Call Draco's setupMPILibraries() to discover and configure MPI for use.
# 2. For MSVC+MSYS-gfortran, do some extra MPI setup to help this project find MPI's headers and
#    libraries.
# 3. Sets and returns CAFS_MPI_DEPS
#
# ------------------------------------------------------------------------------------------------ #
macro(CAFS_Fortran_dir_MPI_setup)

  if("${DRACO_C4}" STREQUAL "MPI")
    # CAFS setup unique to this directory; vendor discovery
    setupmpilibraries()

    # Link to msys64 formatted msmpi.a instead of VS formatted msmpi.lib/dll, and ensure that the
    # compiler option '-frange-check' is disabled.
    cafs_fix_mpi_library()

    add_definitions(-DC4_MPI)
  endif()

endmacro(CAFS_Fortran_dir_MPI_setup)

# ------------------------------------------------------------------------------------------------ #
# End of CMakeAddFortranSubdirectory.cmake
# ------------------------------------------------------------------------------------------------ #
