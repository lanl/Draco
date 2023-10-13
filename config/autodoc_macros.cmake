# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/autodoc_macros.cmake
# author Kelly G. Thompson, kgt@lanl.gov
# date   Wednesday, Nov 14, 2018, 19:01 pm
# brief  Provide extra macros to simplify CMakeLists.txt for autodoc directories.
# note   Copyright (C) 2018-2023 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)
set(draco_config_dir
    ${CMAKE_CURRENT_LIST_DIR}
    CACHE INTERNAL "")

# ------------------------------------------------------------------------------------------------ #
# Set the target directory where the html files wll be created.
# ------------------------------------------------------------------------------------------------ #
function(set_autodocdir)
  # if AUTODOCDIR is set in environment (or make command line), create a CMake variable with this
  # value.
  if(DEFINED ENV{AUTODOCDIR})
    set(AUTODOCDIR "$ENV{AUTODOCDIR}")
  endif()
  # if AUTODOCDIR is set, use it, otherwise, default to CMAKE_INSTALL_PREFIX.
  if(DEFINED AUTODOCDIR)
    set(DOXYGEN_OUTPUT_DIR
        "${AUTODOCDIR}"
        PARENT_SCOPE)
  else()
    set(DOXYGEN_OUTPUT_DIR
        ${CMAKE_INSTALL_PREFIX}/autodoc
        PARENT_SCOPE)
  endif()
endfunction()

# ------------------------------------------------------------------------------------------------ #
# Build a list of directories that include sources that doxygen should examine.
# ------------------------------------------------------------------------------------------------ #
function(set_doxygen_input)

  set(DOXYGEN_INPUT "${PROJECT_SOURCE_DIR}/autodoc" "${PROJECT_BINARY_DIR}/autodoc")
  set(DOXYGEN_EXAMPLE_PATH "")

  # BUG: Move this list generation into component_macros.cmake so that inactive packages are not
  # included in this list.
  file(GLOB package_list ${PROJECT_SOURCE_DIR}/src/*)
  foreach(package ${package_list})
    if(EXISTS ${package}/CMakeLists.txt)
      list(APPEND DOXYGEN_INPUT "${package}")
    endif()
    if(EXISTS ${package}/test/CMakeLists.txt)
      list(APPEND DOXYGEN_INPUT "${package}/test")
      list(APPEND DOXYGEN_EXAMPLE_PATH "${package}/test")
    endif()
    if(EXISTS ${package}/autodoc)
      list(APPEND DOXYGEN_INPUT "${package}/autodoc")
    endif()
  endforeach()

  # Also look for files that have been configured (.in files) and placed in the BINARY_DIR.
  file(GLOB package_list ${PROJECT_BINARY_DIR}/src/*)
  foreach(package ${package_list})
    # pick up processed .dcc files
    if(EXISTS ${package}/autodoc)
      list(APPEND DOXYGEN_INPUT "${package}/autodoc")
    endif()
  endforeach()

  # convert list of directories into a space delimited string
  unset(temp)
  foreach(dir ${DOXYGEN_INPUT})
    set(temp "${temp} ${dir}")
  endforeach()
  set(DOXYGEN_INPUT
      "${temp}"
      PARENT_SCOPE)

  unset(temp)
  foreach(dir ${DOXYGEN_EXAMPLE_PATH})
    set(temp "${temp} ${dir}")
  endforeach()
  set(DOXYGEN_EXAMPLE_PATH
      "${temp}"
      PARENT_SCOPE)
  unset(temp)

  if(${DRACO_DBC_LEVEL} GREATER 0)
    set(DOXYGEN_ENABLED_SECTIONS
        "REMEMBER_ON"
        PARENT_SCOPE)
  endif()

  # Tell doxygen where Draco's include files are: set( DOXYGEN_INCLUDE_PATH
  # "${CMAKE_INSTALL_PREFIX}/include" PARENT_SCOPE)

endfunction()

# ------------------------------------------------------------------------------------------------ #
# Build a list of directories that include images that doxygen should be able to find.
# ------------------------------------------------------------------------------------------------ #
function(set_doxygen_image_path)
  # Tell doxygen where image files are located so they can be copied to the output directory.
  #
  # The list of source files (this variable also set by
  # comonent_macros.cmake::process_autodoc_pages()
  list(APPEND DOXYGEN_IMAGE_PATH "${PROJECT_SOURCE_DIR}/autodoc"
       "${PROJECT_SOURCE_DIR}/autodoc/html")
  list(REMOVE_DUPLICATES DOXYGEN_IMAGE_PATH)
  # convert list of image directories into a space delimited string
  unset(temp)
  foreach(image_dir ${DOXYGEN_IMAGE_PATH})
    set(temp "${temp} ${image_dir}")
  endforeach()
  set(DOXYGEN_IMAGE_PATH
      "${temp}"
      PARENT_SCOPE)
  unset(temp)
endfunction()

# ------------------------------------------------------------------------------------------------ #
# Set the number of cpu threads to use when generating the documentation.
# ------------------------------------------------------------------------------------------------ #
function(set_doxygen_dot_num_threads)
  # Doxygen only allows 32 threads max
  if(NOT DEFINED MPIEXEC_MAX_NUMPROCS OR "${MPIEXEC_MAX_NUMPROCS}empty" STREQUAL "empty")
    set(DOXYGEN_DOT_NUM_THREADS
        1
        PARENT_SCOPE)
  elseif("${MPIEXEC_MAX_NUMPROCS}" GREATER 32)
    set(DOXYGEN_DOT_NUM_THREADS
        32
        PARENT_SCOPE)
  else()
    set(DOXYGEN_DOT_NUM_THREADS
        "${MPIEXEC_MAX_NUMPROCS}"
        PARENT_SCOPE)
  endif()
  # hack in a couple of other settings based on the version of doxygen discovered.
  if(${DOXYGEN_VERSION} VERSION_GREATER 1.8.14)
    set(DOXYGEN_HTML_DYNAMIC_MENUS
        "HTML_DYNAMIC_MENUS = YES"
        PARENT_SCOPE)
  endif()
  if(${DOXYGEN_VERSION} VERSION_LESS 1.8.17)
    set(PERL_PATH
        "PERL_PATH = /usr/bin/perl"
        PARENT_SCOPE)
  endif()

  # Escalate doxygen warnings into errors for CI builds
  if(DEFINED ENV{CI})
    set(DOXYGEN_WARN_AS_ERROR
        "YES"
        PARENT_SCOPE)
  else()
    set(DOXYGEN_WARN_AS_ERROR
        "NO"
        PARENT_SCOPE)
  endif()
endfunction()

# ------------------------------------------------------------------------------------------------ #
# Generate and install HTML support files
#
# Requires the following variables to be set:
#
# * PROJECT_SOURCE_DIR - Always provided by CMake (but should point to the to top level source
#   directory.
# * DOXYGEN_OUTPUT_DIR - Directory where html code will be written by doxygen. Actual location of
#   HTML files will be ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}.
# * DOXYGEN_HTML_OUTPUT - The project name in all lowercase.
#
# ------------------------------------------------------------------------------------------------ #
macro(doxygen_provide_support_files)

  add_custom_command(
    OUTPUT "${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/footer.html"
    COMMAND
      "${CMAKE_COMMAND}" -DINFILE="${PROJECT_SOURCE_DIR}/autodoc/html/footer.html.in"
      -DOUTFILE="${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/footer.html" -P
      "${draco_config_dir}/configureFileOnMake.cmake"
    DEPENDS "${PROJECT_SOURCE_DIR}/autodoc/html/footer.html.in"
    COMMENT "Generating footer.html")
  add_custom_command(
    OUTPUT "${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/header.html"
    COMMAND
      "${CMAKE_COMMAND}" -DINFILE="${PROJECT_SOURCE_DIR}/autodoc/html/header.html.in"
      -DOUTFILE="${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/header.html" -P
      "${draco_config_dir}/configureFileOnMake.cmake"
    DEPENDS "${PROJECT_SOURCE_DIR}/autodoc/html/header.html.in"
    COMMENT "Generating header.html")

  if(EXISTS "${draco_config_dir}/doxygen.css")
    # use Draco's version of the style sheet
    set(doxygen_ccs_file "${draco_config_dir}/doxygen.css")
  elseif(EXISTS "${PROJECT_SOURCE_DIR}/autodoc/html/doxygen.css")
    # use Draco's version of the style sheet
    set(doxygen_ccs_file "${PROJECT_SOURCE_DIR}/autodoc/html/doxygen.css")
  else()
    message(
      FATAL_ERROR
        "I can't find a style sheet to install for autodoc. Expected to find doxygen.css at either:
- ${draco_config_dir}/, or
- ${PROJECT_SOURCE_DIR}/autodoc/html/")
  endif()
  add_custom_command(
    OUTPUT "${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/doxygen.css"
    COMMAND
      "${CMAKE_COMMAND}" -DINFILE="${doxygen_ccs_file}"
      -DOUTFILE="${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/doxygen.css" -P
      "${draco_config_dir}/configureFileOnMake.cmake"
    DEPENDS "${doxygen_css_file}"
    COMMENT "Generating doxygen.css")
endmacro()

# ------------------------------------------------------------------------------------------------ #
# Create a string to locate Draco.tag
# ------------------------------------------------------------------------------------------------ #
function(set_doxygen_tagfiles)
  # Create links to Draco autodoc installation.
  unset(DRACO_TAG_FILE CACHE)
  find_file(DRACO_TAG_FILE Draco.tag HINTS ${DOXYGEN_OUTPUT_DIR} ${DOXYGEN_OUTPUT_DIR}/..)
  get_filename_component(DRACO_AUTODOC_DIR ${DRACO_TAG_FILE} PATH)
  file(RELATIVE_PATH DRACO_AUTODOC_HTML_DIR ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}
       ${DRACO_AUTODOC_DIR}/draco)
  set(TAGFILES
      "${DRACO_TAG_FILE}=${DRACO_AUTODOC_HTML_DIR}"
      PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------------------------ #
# End config/autodoc_macros.cmake
# ------------------------------------------------------------------------------------------------ #
