#!/bin/bash
###
#OAR -l /nodes=1/core=8,walltime=06:00:00
#OAR -O /temp_dd/igrida-fs1/pderian/matlab.%jobid%.output
#OAR -E /temp_dd/igrida-fs1/pderian/matlab.%jobid%.error
###
# path to sqgmu root directory
SQGMU_ROOT=~/matlab/sqgmu 

echo "*** Log for job #$OAR_JOB_ID ***"
echo "Running on:"
cat $OAR_NODEFILE
echo "Beginning job - " $(date)

# make script aware of "module"
source /etc/profile.d/modules.sh

# install MCR
echo "Installing MCR..."
source ~/matlab/setup_MCR.sh
echo "MCR installed in "$MCR_ROOT

# load MATLAB module
module load matlab/$MATLAB_RELEASE
# compile code
# [TODO] destination in SCRATCH?
echo "Compiling M-code..."
cd "$SQGMU_ROOT"/compile
matlab -nojvm -r compile
# unload module
module unload matlab/$MATLAB_RELEASE

# start
echo "Launching M-code..."
cd $SCRATCHDIR
cd "$SQGMU_ROOT"/compile/run_sqgmu.sh $MCR_RUN

# clean MCR
echo "Removing MCR..."
rm -rf $MCR_ROOT

echo
echo "*** job $OAR_JOB_ID complete ***"

