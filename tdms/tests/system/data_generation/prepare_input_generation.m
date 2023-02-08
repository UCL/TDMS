%% Setup steps for generating the input file to a TDMS system test

% Relative path to the matlab functions that generate the input data
addpath('../../../../matlab');

% Create directory to place input files in, if it doesn't exist already
if ~exist('in', 'dir')
    mkdir('in');
end

% clear variables that may polute workspace
clear; close all;
