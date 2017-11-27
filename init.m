%%%%%%%%%%%%%%%%%%
% Initialisation
%%%%%%%%%%%%%%%%%%

% Cleaning
clear;
% clear all;
close all;clc;

% Debug
dbstop if error;

% Paths
fct = genpath([ pwd '/functions' ]);
addpath(pwd)
addpath(fct)

% Cleaning
clear fct
home;
