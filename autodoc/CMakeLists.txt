#--------------------------------------------*-cmake-*---------------------------------------------#
# file   draco/autodoc/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Oct 14
# brief  Instructions for building Doxygen documentation
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

# 2018-03-10 (KT): As of cmake-3.9.0, we have the option of using the cmake function
# doxygen_add_docs() instead of the complex logic below.  We would lose a few features, but wouldn't
# need to maintain this file.  Ref. https://cmake.org/cmake/help/v3.10/module/FindDoxygen.html

# config/vendor_libraries.cmake is responsible for detecting if doxygen is available.

if(NOT DOXYGEN_FOUND)
  return()
endif()

include(autodoc_macros)
set_autodocdir() # --> AUTODOCDIR

# Source directories that should be examined by doxygen Sets ${DOXYGEN_INPUT}.  This value is
# embedded in the doxygen_config file below.
set_doxygen_input() # --> DOXYGEN_INPUT

# Tell doxygen where image files are located so they can be copied to the output directory.
set_doxygen_image_path() # --> DOXYGEN_IMAGE_PATH

# Doxygen only allows 32 threads max
set_doxygen_dot_num_threads() # --> DOXYGEN_DOT_NUM_THREADS --> DOXYGEN_HTML_DYNAMIC_MENUS (v
                              # 1.8.14+) --> DOXYGEN_WARN_AS_ERROR --> PERL_PATH

set(DOXYGEN_HTML_OUTPUT draco)
string(CONCAT project_brief "Object-oriented components for numerically intensive, parallel, "
              "radiation transport applications.")
string(REPLACE " " "___" project_brief ${project_brief})
string(REPLACE " " "___" DOXYGEN_HTML_DYNAMIC_MENUS "${DOXYGEN_HTML_DYNAMIC_MENUS}")
string(REPLACE " " "___" PERL_PATH "${PERL_PATH}")

# doxygen_config.in needs the following symbols:
add_custom_command(
  OUTPUT "${PROJECT_BINARY_DIR}/autodoc/doxygen_config"
  COMMAND
    "${CMAKE_COMMAND}" -DINFILE="${PROJECT_SOURCE_DIR}/config/doxygen_config.in"
    -DOUTFILE="${PROJECT_BINARY_DIR}/autodoc/doxygen_config" -DPROJECT_NAME=${PROJECT_NAME}
    -Dpackage=${PROJECT_NAME} -DPROJECT_NUMBER="${Draco_VERSION_FULL}"
    -Dnumber="${Draco_VERSION_FULL}" -DOUTPUT_DIRECTORY="${DOXYGEN_OUTPUT_DIR}"
    -Ddoxygen_output_top="${DOXYGEN_OUTPUT_DIR}" -DINPUT=${DOXYGEN_INPUT}
    -Ddoxygen_input=${DOXYGEN_INPUT} -DDOXYGEN_EXAMPLE_PATH=${DOXYGEN_EXAMPLE_PATH}
    -DHTML_OUTPUT=${DOXYGEN_HTML_OUTPUT} -Ddoxygen_html_output=${DOXYGEN_HTML_OUTPUT}
    -DDOTFILE_DIRS="" -Ddotpath="" -Dproject_brief=\"${project_brief}\" -Dlatex_yes_no=NO
    -Ddoxygen_latex_output=${DOXYGEN_HTML_OUTPUT}/../latex
    -Dhtmlheader="${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/header.html"
    -Dhtmlfooter="${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/footer.html"
    -Dproject_logo="${PROJECT_SOURCE_DIR}/autodoc/html/dragon.jpg"
    -Dnumthreads=${DOXYGEN_DOT_NUM_THREADS} -DDOXYGEN_IMAGE_PATH=${DOXYGEN_IMAGE_PATH}
    -DDOXYGEN_HTML_DYNAMIC_MENUS="${DOXYGEN_HTML_DYNAMIC_MENUS}"
    -DDOXYGEN_WARN_AS_ERROR="${DOXYGEN_WARN_AS_ERROR}"
    -DDOXYGEN_ENABLED_SECTIONS="${DOXYGEN_ENABLED_SECTIONS}"
    -DDOXYGEN_INCLUDE_PATH="${DOXYGEN_INCLUDE_PATH}" -P
    "${PROJECT_SOURCE_DIR}/config/configureFileOnMake.cmake"
  DEPENDS "${PROJECT_SOURCE_DIR}/config/doxygen_config.in"
          "${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/header.html"
          "${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/footer.html"
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/autodoc
  COMMENT "")

# Create header.html, footer.html and doxygen.css.
doxygen_provide_support_files()

add_custom_command(
  OUTPUT "${PROJECT_BINARY_DIR}/autodoc/mainpage.dcc"
  COMMAND
    "${CMAKE_COMMAND}" -DINFILE="${PROJECT_SOURCE_DIR}/autodoc/mainpage.dcc.in"
    -DOUTFILE="${PROJECT_BINARY_DIR}/autodoc/mainpage.dcc" -DCOMP_LINKS=${COMP_LINKS}
    -DPACKAGE_LINKS=${PACKAGE_LINKS} -DDraco_VERSION_FULL="${Draco_VERSION_FULL}"
    -DDRACO_INFO="$<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>"
    -DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}" -P
    "${PROJECT_SOURCE_DIR}/autodoc//generate_mainpage_dcc.cmake"
  DEPENDS "${PROJECT_SOURCE_DIR}/autodoc/mainpage.dcc.in" Exe_draco_info
  COMMENT "")
add_custom_command(
  OUTPUT "${PROJECT_BINARY_DIR}/autodoc/build_system.dcc"
  COMMAND
    "${CMAKE_COMMAND}" -DINFILE="${PROJECT_SOURCE_DIR}/autodoc/build_system.dcc.in"
    -DOUTFILE="${PROJECT_BINARY_DIR}/autodoc/build_system.dcc" -DCOMP_LINKS=${COMP_LINKS}
    -DPACKAGE_LINKS=${PACKAGE_LINKS} -P "${PROJECT_SOURCE_DIR}/config/configureFileOnMake.cmake"
  DEPENDS "${PROJECT_SOURCE_DIR}/autodoc/build_system.dcc.in"
  COMMENT "")

# The actual doxygen command that will be run
add_custom_command(
  OUTPUT ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/always_out_of_date.html
  BYPRODUCTS ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/index.html
  COMMAND "${DOXYGEN_EXECUTABLE}" ${PROJECT_BINARY_DIR}/autodoc/doxygen_config
  DEPENDS ${PROJECT_BINARY_DIR}/autodoc/doxygen_config ${PROJECT_BINARY_DIR}/autodoc/mainpage.dcc
          ${PROJECT_BINARY_DIR}/autodoc/build_system.dcc
          ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/doxygen.css
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/autodoc
  COMMENT "")

# Target for building doxygen files.
add_custom_target(
  autodoc
  DEPENDS ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/always_out_of_date.html
  COMMENT
    "Building Doxygen based autodocs (HTML)...

*** Autodoc HTML pages are at
  file://${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/index.html

")

# To Do: 1. How to 'make clean'?

# ------------------------------------------------------------------------------------------------ #
# End autodoc/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
