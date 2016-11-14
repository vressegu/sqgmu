%%%%%%%%%%%%%%%%%%
% Initialisation
%%%%%%%%%%%%%%%%%%
%
% Modified by P. DERIAN 2016-10-10
%   - added isdeployed for compilation support

% Cleaning
clear all;
close all;clc;

% Debug
dbstop if error;

% Paths
if ~isdeployed
    % Note: isdeployed necessary for compiled code (addpath won't work)
    fct = genpath([ pwd '/functions' ]);
    addpath(pwd)
    addpath(fct)
    clear fct
end
    
home;