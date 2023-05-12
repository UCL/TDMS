%% This script creates the bad_class_data.mat file for use in the hdf5 tests for reading in TDMS objects.
% These variables are used to check failure-cases for when the reader methods should detect an error.
close all;
clear;

%% f_vec
% f_vec is a struct with two fields, fx_vec and fy_vec.
f_vec = struct();
f_vec.fx_vec = zeros(2,2);      % Non-2D elements are not expected!
f_vec.fy_vec = zeros(2,1);      % This is permitted, but the former field should error

%% phasorsurface
% This array is read into the Cuboid class. It is just an array of 6 integers (stored as doubles OFC) that correspond to Yee cell indices in the various axial directions
phasorsurface = [1; 4; 2; 5; 3; 6; 7;]; % 7 elements should throw an error

%% save variables to the file we need
% Save the files to the expected filename for the unit tests to read the
% data back in.

save("unit_test_data/bad_class_data.mat", "-v7.3");
