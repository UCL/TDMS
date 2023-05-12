%% This script creates the structure_array.mat file for use in the hdf5 tests for reading in structure arrays.
close all;
clear;

example_struct = struct();

example_struct.double_no_decimal = 1.;
example_struct.double_half = 0.5;
example_struct.string = 'tdms';
example_struct.boolean = true;
example_struct.uint_345 = uint8(ones(3, 4, 5));
example_struct.double_22 = [0.25, 0.5; 0.75, 1.];
example_struct.complex_22 = [0., -1.i; 1.i, 0.];

%% save variables to the file we need
% Save the files to the expected filename for the unit tests to read the
% data back in.

save("structure_array.mat", "example_struct", "-v7.3");
