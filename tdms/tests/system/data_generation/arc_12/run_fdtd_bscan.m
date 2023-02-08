%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section generates the files used as input to the executeable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../');
prepare_input_generation;

%% Generate the input file
% By default, intmethod is set to 1 if not present in the input file.
% Ergo, we are creating the cubic interpolation input file.

%start by defining the coordinates of the computational grid
[x,y,z,lambda] = fdtd_bounds('fdtd_input_file.m');

%15 micron radius cylinder
rad = 1e-6;
%refractive index of cylinder
refind = 1.42;

%insert a sphere at the origin
[X,Y,Z] = ndgrid(x,y,z);

%generate scattering matrix
I = zeros(size(X));
%set all Yee cells within the cylinder to have index of 1
I( (X.^2 + Y.^2 + Z.^2) < rad^2 ) = 1;
I( (end-3):end,1,:) = 0;
I( :, 1, (end-3):end) = 0;
inds = find(I(:));
[ii,jj,kk] = ind2sub(size(I), inds);
composition_matrix = [ii jj kk ones(size(ii))];
material_matrix = [1 refind^2 1 0 0 0     0     0     0 0 0];
save('gridfile_sph', 'composition_matrix', 'material_matrix');

%setup free space matrix
composition_matrix = [];
save('gridfile_fs', 'composition_matrix', 'material_matrix');

%generate the illumination file
iteratefdtd_matrix('illumination_input_file.m','illsetup','illfile','gridfile_sph.mat','');

%generate tdms executable files
iteratefdtd_matrix('fdtd_input_file.m','filesetup','in/fdtd_sph_cubic','gridfile_sph.mat','illfile');
iteratefdtd_matrix('fdtd_input_file.m','filesetup','in/fdtd_fs_cubic','gridfile_fs.mat','illfile');

%% Band-limited interpolation input file
% This can be generated from the cubic input file that we just created.
% The setup is identical, we just need to change the intmethod flag to 2.
% Thus, we can simply clear (all variables), load the file we just created,
% set intmethod to 2, then save the workspace as the new input file.

clear; close all;

set_interpolation_method('in/fdtd_sph_cubic', 'in/fdtd_sph_bli', 2);
set_interpolation_method('in/fdtd_fs_cubic', 'in/fdtd_fs_bli', 2);

%% Cleanup when done
clear;
