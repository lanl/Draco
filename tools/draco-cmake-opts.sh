#!/bin/bash
#--------------------------------------------------------------------------------------------------#
# CMake options for Draco releases
#--------------------------------------------------------------------------------------------------#

export OPTIMIZE_ON="-DCMAKE_BUILD_TYPE=Release -DDRACO_LIBRARY_TYPE=SHARED -DDRACO_ROUNDOFF_MODE=FAST"
export OPTIMIZE_OFF="-DCMAKE_BUILD_TYPE=Debug  -DDRACO_LIBRARY_TYPE=SHARED -DDRACO_ROUNDOFF_MODE=ACCURATE"
export OPTIMIZE_RWDI="-DCMAKE_BUILD_TYPE=RELWITHDEBINFO -DDRACO_LIBRARY_TYPE=SHARED -DDRACO_DBC_LEVEL=15 -DDRACO_ROUNDOFF_MODE=ACCURATE"

export LOGGING_ON="-DDRACO_DIAGNOSTICS=7 -DDRACO_TIMING=1"
export LOGGING_OFF="-DDRACO_DIAGNOSTICS=0 -DDRACO_TIMING=0"
