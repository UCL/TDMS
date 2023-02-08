%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This m-file is a test case for steady state operation
%
%This section generates the files used as input to the executeable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../');
prepare_input_generation;

%% Generate the input file
% By default, intmethod is set to 1 if not present in the input file.
% Ergo, we are creating the cubic interpolation input file.

%start by defining the coordinates of the computational grid
[x,y,z,lambda] = fdtd_bounds('pstd_input_file.m');

%refractive index of scatterer
refind = 1.42;

%insert a sphere at the origin
[X,Y,Z] = ndgrid(x,y,z);

%generate scattering matrix with a single scatterer at the origin
I = zeros(size(X));
I( ( (X==0) & (Y==0) & (Z==0) ) ) = 1;

inds = find(I(:));
[ii,jj,kk] = ind2sub(size(I), inds);
composition_matrix = [ii jj kk ones(size(ii))];
material_matrix = [1 refind^2 1 0 0 0     0     0     0 0 0];
save('gridfile_sc', 'composition_matrix', 'material_matrix');

%setup free space matrix
composition_matrix = [];
save('gridfile_fs', 'composition_matrix', 'material_matrix');

%generate tdms executable files
iteratefdtd_matrix('pstd_input_file.m','filesetup','in/pstd_sc_cubic','gridfile_sc.mat','');
iteratefdtd_matrix('pstd_input_file.m','filesetup','in/pstd_fs_cubic','gridfile_fs.mat','');

%% Band-limited interpolation input file
% This can be generated from the cubic input file that we just created.
% The setup is identical, we just need to change the intmethod flag to 2.
% Thus, we can simply clear (all variables), load the file we just created,
% set intmethod to 2, then save the workspace as the new input file.

clear; close all;

set_interpolation_method('in/pstd_sc_cubic', 'in/pstd_sc_bli', 2);
set_interpolation_method('in/pstd_fs_cubic', 'in/pstd_fs_bli', 2);

%% Cleanup when done
clear;
