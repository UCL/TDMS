%% Setup steps for generating the input file to a TDMS system test

% Relative path to the shared matlab functions that generate the input data
% NOTE: If there are function-name conflicts, MATLAB uses the version of the function whose file is located closest to the current working directory. This means that function files in the data_generation folders for each test serve as overrides where necessary for the methods that are defined in the shared functionality folder.
addpath('../../../../matlab');

% Create directory to place input files in, if it doesn't exist already
if ~exist('in', 'dir')
    mkdir('in');
end

% clear variables that may polute workspace
clear; close all;
