%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates the files used as input to the executeable

%Ensure that you are running MATLAB in the examples/arc_01 directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The name of the input file to be read
input_filename = 'arc_01_example_input.m';

% Add the tdms/tests/system/data/input_generation/matlab directory to the
% current MATLABPATH.
% Functions such as iteratefdtd_matrix (and its dependents) must be on the
% MATLABPATH in order to allow them to be called.
addpath('../../tdms/tests/system/data/input_generation/matlab/');


%% 4. Material Specification
% We now setup the quantities specified in section 4 of the PDF
% documentation.
% These can be setup in an existing .mat or .m file if the user desires.
% In which case, it is necessary to include its name as the material_file input to iteratefdtd_matrix.
% In this example, we will use the functions in
% tdms/tests/system/data/input_generation/matlab/ to generate the material
% inputs at runtime (IE, now).
% For more information about the functions called here, see the
% documentation in the source files located in the aforementioned
% directory.

% The shape of the scattering object(s) that will be present in the simulation
% 'cyl': cylinder, 'sph': sphere, 'fs': freespace, 'sc': point-source
obstacle_cyl = 'cyl';
obstacle_fs = 'fs';
% The radius of the scattering obstacle (if it has one) in microns
radius_cyl = 15e-6;
% The (relative) refractive index of the scattering obstacle
refind = 1.42;

% See S4 of the PDF documentation for explanations of the composition_matrix and material_matrix.
% See the docstring in composition_matrix_builder.m for an explanation of
% the function.
composition_matrix_cyl = composition_matrix_builder(input_filename, obstacle_cyl, radius_cyl);
material_matrix_cyl = [1 refind^2 1 0 0 0     0     0     0 0 0];

composition_matrix_fs = composition_matrix_builder(input_filename, obstacle_fs, 0.);
material_matrix_fs = [1 refind^2 1 0 0 0     0     0     0 0 0];

% Save the material files for the two scattering objects, so we can pass
% them to iteratefdtd_matrix.
% We could have created these files via a separate script, or re-used
% material files from another simulation if suitable.
gridfile_cyl = 'gridfile_cyl.mat';
save(gridfile_cyl, 'composition_matrix_cyl', 'material_matrix_cyl', '-v7.3');
gridfile_fs = 'gridfile_fs.mat';
save(gridfile_fs, 'composition_matrix_fs', 'material_matrix_fs', '-v7.3');

%% 3.2.2 Source file specification
% If our input_filename did not specify that we wish to use a matlab
% function to compute the incident fields on the interface planes, we could
% provide a pre-made .mat file containing this information to
% iteratefdtd_matrix by populating the illfile variable.
% Since we do not wish to do this, we leave the name empty and
% iteratefdtd_matrix will infer from the input file that we want to compute
% the source fields from scratch.
illfile = '';

%% Use iteratefdtd_matrix to setup the input file(s).
% Create the input to the simulation that models scattering from a
% cylindrical obstacle.
iteratefdtd_matrix(input_filename,'filesetup','cylinder_input.mat',gridfile_cyl,illfile);
% Create the input to the simulation that has no obstacle, and only models
% free space.
iteratefdtd_matrix(input_filename,'filesetup','freespace_input.mat',gridfile_fs,illfile);

%% Run the tdms executable
% These commands can be carried out via the CLI.

system('!tdms in_pstd_fs.mat out_pstd_fs.mat');
system('!tdms in_pstd_cyl.mat out_pstd_cyl.mat');

%plot the data
dat_cyl = load('out_pstd_cyl');
dat_fs = load('out_pstd_fs');

figure(1);clf;
subplot(2,1,1);
imagesc(dat_fs.z_i,dat_fs.x_i,abs(squeeze(dat_fs.Ex_i)));
axis equal;
title('Focussed beam in free space');

subplot(2,1,2);
imagesc(dat_fs.z_i,dat_fs.x_i,abs(squeeze(dat_cyl.Ex_i)));
title('Focussed beam with scattering cylinder');
axis equal;

%%
figure;
imagesc(dat_cyl.z_i, dat_cyl.x_i, abs(squeeze(dat_cyl.Hz_out)));
axis equal;
title('Normalised Hz component of the scattered field from a cylinder');
