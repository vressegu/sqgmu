%%%%%%%%%%%%%%%%%%
% Initialisation
%%%%%%%%%%%%%%%%%%
%
% Modified by P. DERIAN 2016-10-10
%   - added "isdeployed" for compilation support
% Modified by P. DERIAN 2017-03-10
%   - removed "clear all"

% Cleaning
close all;clc;

% Note: isdeployed necessary for compiled code (addpath won't work)
if ~isdeployed
    % Debug
    dbstop if error;

    % Paths
    fct = genpath([ pwd '/functions' ]);
    addpath(pwd)
    addpath(fct)
    clear fct
end
    
home;