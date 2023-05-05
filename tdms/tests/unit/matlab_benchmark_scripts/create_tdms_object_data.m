%% This script creates the class_data.mat file for use in the hdf5 tests for reading in TDMS objects.
close all;
clear;

%% interface
% interface is a struct with fields {I, J, K} {0, 1}, corresponding to the planes at which a source field is introduced.
% Each field is a 1 - by - 2 array of doubles; the first element being the
% {I, J, K} index of the Yee cell in which the particular plane lies.
% The second element is cast to a bool and indicates whether any boundary
% conditions are to be applied on said plane.
interface = struct();
interface.I0 = [1 0]; % I0 in the I = 1 plane, no source condition
interface.I1 = [4 1]; % I1 in the I = 4 plane, source condition applied
interface.J0 = [2 0]; % J0 in the J = 2 plane, no source condition
interface.J1 = [5 0]; % J1 in the J = 5 plane, no source condition
interface.K0 = [3 1]; % K0 in the K = 3 plane, source condition applied
interface.K1 = [6 1]; %K1 in the K=6 plane, source condition applied

%% f_vec & f_vec_bad
% f_vec is a struct with two fields, fx_vec and fy_vec.
% These are each 1D vectors of doubles, usually they have the same length
% however for our testing purposes we will make them different lengths.
f_vec = struct();
f_vec.fx_vec = [0.25 0.5 0.75 1.];          % Row vector w/ 4 elements
f_vec.fy_vec = [-0.25; -0.5; -0.75; -1.];   % Col vector w/ 4 elements

% f_vec_bad will be used as the failure case in our tests.
f_vec_bad = struct();
f_vec_bad.fx_vec = zeros(2,2);      % Non-2D elements are not expected!
f_vec_bad.fy_vec = zeros(2,1);      % This is permitted, but the former field should error

%% phasorsurface
% This array is read into the Cuboid class. It is just an array of 6 integers (stored as doubles OFC) that correspond to Yee cell indices in the various axial directions
phasorsurface = [1 4 2 5 3 6];

%% save variables to the file we need
% Save the files to the expected filename for the unit tests to read the
% data back in.

save("matlab_data/class_data.mat", "-v7.3");
