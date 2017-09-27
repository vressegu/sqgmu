%% Compilation script for SQGMU
% Compiles the SQGMU project for a licence-free execution using the
% Matlab Component runtime (MCR).
%
% Note: this script is meant to be started in command line, so as to exit the session
% at the end and release the Matlab Compiler (MCC) token. E.g. with:
%   matlab -nodesktop -nosplash -r compile
% It requires a Matlab license token as well as a Compiler token.
%
% Note: when compiling remotely, the script may fail if no screen is active.
%  (known MATLAB bug) On Os X, a workaround is to run: 
%       caffeinate -u -t 1
%   which wakes up the screen.
%   See https://fr.mathworks.com/matlabcentral/answers/202908-how-can-i-run-a-batch-program-on-mac-osx-that-produces-high-quality-figures-when-there-is-no-display
%
% Written by P. DERIAN 2016-10-12.
fprintf(1, 'Beginning compilation of SQGMU...\n');

% Add the project functions to the path
addpath( genpath('../functions') );

% Notes on compilation:
% -I to add path to the compiler
% -a <some_file> to add specific files (e.g. data) to the archive
% -o <output_file> for the name of the output archive
% for igrida:
% -R -nojvm to tell the Matlab runtime to disable the jvm (no graphics!)
% -R -singleCompThread to tell runtime not to use multithreading
%
% However: the jvm is mandatory for the Parallel Toolbox:
% "Java must be initialized in order to use the Parallel Computing Toolbox."

% Notes on run:
% use ./run_sqgmu <MCR_PATH/vxx> <SQGMU_args>
% see also the readme.txt generated automatically

% Notes on parallel:
% see http://fr.mathworks.com/help/compiler/use-the-parallel-computing-toolbox.html
% for the parallel computing stuff.
%
% So: to use a pool, the JVM is mandatory.
mcc -m -I .. ...
       -I /Applications/MATLAB_R2015b.app/toolbox/Wavelab850/Orthogonal ...
       -a ../functions/output/BuYlRd.mat main -o ...
       sqgmu

   
% Job's done
fprintf(1, 'Compilation complete, exiting.\n');
fprintf(1, 'Try: ./run_sqgmu.sh <PATH_TO_MCR>\n');
exit;

