#!/bin/bash
## Compilation script for SQGMU
# Compiles the SQGMU project for a licence-free execution using the 
# Matlab Component runtime (MCR). Requires mcc to be in the $PATH.
#
# Usage: sh compile.sh [PATH_TO_MCC]
#
# This script requires a Matlab Compiler token, which is awarded for 30 min.
# This sucks.
#
# Written by P. DERIAN 2017-02-11.
echo 'Beginning compilation of SQGMU...'

# Note: -I is not recursive, so every directory has to be listed here.
mcc -m                                  \
    -I ..                               \
    -I ../functions                     \
    -I ../functions/advection           \
    -I ../functions/initial_fields      \
    -I ../functions/output              \
    -I ../functions/sigma/misc          \
    -I ../functions/sigma/spectrum      \
    -I ../functions/sigma/SVDnoise      \
    -I ../functions/SQG_sto             \
    -a ../functions/output/BuYlRd.mat   \
    main -o sqgmu
