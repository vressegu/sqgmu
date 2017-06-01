#!/bin/bash
#
# This script install the MCR (Matlab Component Runtime) on the current node
# in /srv/local/MCR_<HOSTNAME>, available in $MCR_ROOT after successful execution.
#
# M-code must be compiled with Matlab release corresponding to the MCR's.
# This release name (version number) is available in $MATLAB_RELEASE ($MATLAB_VERSION)
# after execution of this script.
#
# MCR executables are installed in $MCR_RUN, this is the path to be supplied
# to Matlab's automatically generated run_<exec_name>.sh
#
# see https://fr.mathworks.com/help/compiler/install-the-matlab-runtime.html#bun27bq-1
# and https://fr.mathworks.com/products/compiler/mcr/
#
# Written by P. DERIAN - 2016-11.
#
###

HOST=$(hostname -s)
MCR_ROOT=/srv/local/MCR_"$HOST"
MCR_INSTALL=$MCR_ROOT/installer
MCR_ARCHIVE=/soft/matlab/toolbox/compiler/deploy/glnxa64/MCRInstaller.zip #note: a64 for AMD64 architectures
MCR_LOG=$MCR_ROOT/log_install.txt

# create root and installation directory
echo "Creating $MCR_ROOT" 
mkdir -p $MCR_ROOT

# unzip installation program
unzip $MCR_ARCHIVE -d $MCR_INSTALL > $MCR_LOG

# start the installation
echo "Beginning installation"
$MCR_INSTALL/install -mode silent -agreeToLicense yes -destinationFolder $MCR_ROOT >> $MCR_LOG

# upon completion, check the MCR version
MCR_VERSION="v"$(grep "Installing Product: MATLAB Compiler Runtime - Core" $MCR_LOG | sed 's|.*Core \([0-9]*\)\.\([0-9]*\)|\1\2|')
MCR_RUN=$MCR_ROOT/$MCR_VERSION # this is to be passed to the running script

# and the corresponding matlab version
readme=$MCR_INSTALL/readme.txt
MATLAB_VERSION=$(grep -o "MATLAB [0-9]\(.[0-9]\?\)" $readme)
MATLAB_RELEASE=$(grep -m 1 -o "R20[0-9]\{2\}[a-z]" $readme) #look for the first occurence of R20XXy, with XX the year and y the release letter

# display info
echo
echo "Remember to compile MATLAB code with the appropriate version!"
echo "See $readme for details."
echo "This MCR $MCR_VERSION is valid for $MATLAB_VERSION - release $MATLAB_RELEASE."
echo

