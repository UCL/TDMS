%% This script creates the structure_array.mat file for use in the hdf5 tests for reading in structure arrays.
close all;
clear;

%% An example to test how different matlab datatypes are stored in HDF5 files
example_struct = struct();

example_struct.double_no_decimal = 1.;
example_struct.double_half = 0.5;
example_struct.string = 'tdms';
example_struct.boolean = true;
example_struct.uint_345 = uint8(ones(3, 4, 5));
example_struct.double_22 = [0.25, 0.5; 0.75, 1.];
example_struct.complex_22 = [0., -1.i; 1.i, 0.];

%% An example to check that buffers are read in correctly from .mat files
read_in_test = struct();
read_in_test.vector = int32(0:11);
read_in_test.matrix = reshape(0:11, 2, 6);
read_in_test.tensor = reshape(0:11, 2, 3, 2);

%% The two forms of the empty struct / matrix that we might come across
empty_struct = struct([]);
empty_array = [];

%% A non-empty standalone array that we can attempt to read from
nonempty_dataset = 0:4;

%% An example to check that XYZVectors are read in correctly from .mat files
XYZVector = struct();
XYZVector.xyz_x = [0.1 0.2 0.3];
XYZVector.xyz_y = [0.4 0.5 0.6];
XYZVector.xyz_z = [0.7 0.8 0.9];

%% save variables to the file we need
% Save the files to the expected filename for the unit tests to read the
% data back in.

save("unit_test_data/structure_array.mat", "-v7.3");
