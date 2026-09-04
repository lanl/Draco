#--------------------------------------------*-cmake-*---------------------------------------------#
# file   config/lcov-summary.cmake
# brief  From coverage targets, run lcov to generate log files.
# note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

# Example use: cmake -DLCOV=`which lcov` -Doutput_file=build/linux-debug-coverage/coverage.txt \
# -DCOVERAGE_FILE=build/linux-debug-coverage/coverage.info -DMODE=summary \
# -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR} \ -Pconfig/lcov-summary.cmake

set(COVERAGE_FILE ${COVERAGE_FILE})
set(IGNORE_LIST ${IGNORE_LIST})
set(LCOV ${LCOV})
set(MODE ${MODE}) # summary or list
set(output_file ${output_file})
set(PROJECT_SOURCE_DIR ${PROJECT_SOURCE_DIR})

if(${MODE} STREQUAL "list")
  execute_process(
    COMMAND "${LCOV}" --list "${COVERAGE_FILE}" --quiet --ignore-errors unused,unused --include
            ${PROJECT_SOURCE_DIR}/src/* ${IGNORE_LIST}
    OUTPUT_FILE "${output_file}"
    OUTPUT_QUIET)
else()
  execute_process(
    COMMAND "${LCOV}" --summary "${COVERAGE_FILE}" --include ${PROJECT_SOURCE_DIR}/src/*
            ${IGNORE_LIST}
    OUTPUT_FILE "${output_file}"
    OUTPUT_QUIET)
endif()

# ------------------------------------------------------------------------------------------------ #
# End
# ------------------------------------------------------------------------------------------------ #
