%{
An example script showing how to use TDMS to simulate the scattering of a
beam of light due to the presence of a cylindrical scatterer.

This script should be run via MATLAB from the examples/arc_01 directory.

In this script, we setup two simulations over the same grid and with the same source, a pulsed beam of light incident from the K0-plane:
- One simulation has no scatterer in the path of the beam (freespace
obstacle)
- The other simulation has a cylindrical (or since this is 2D, circular)
scattering object in the path of the beam (cylindrical obstacle).
The simulations are 2D, and differ only in the specification of the
material file.
Because of this, the same input file (arc_01_example_input)
can be used for both simulations to setup the grid, timestep, simulation
options, etc.

This script will clear your MATLAB workspace of variables and figures via
the clear and close commands, so please ensure that you have
extracted/saved your workspace before running this script in an existing
instance.
%}
% Clean the environment
clear; close all;

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
% When defining your own obstacles of custom shape, you will need to define
% the composition matrix accordingly.
composition_matrix = composition_matrix_builder(input_filename, obstacle_cyl, radius_cyl);
material_matrix = [1 refind^2 1 0 0 0     0     0     0 0 0];
% Due to the way iteratefdtd_matrix is coded, the material_file must
% contain variables under the names composition_matrix and material_matrix.
% As such, we now save these for the cylinder now, then overwrite them and
% save to a different file for the freespace simulation
gridfile_cyl = 'gridfile_cyl.mat';
save(gridfile_cyl, 'composition_matrix', 'material_matrix', '-v7.3');

% Now produce the material file for the freespace simulation
composition_matrix = composition_matrix_builder(input_filename, obstacle_fs, 0.);
material_matrix = [1 refind^2 1 0 0 0     0     0     0 0 0];
gridfile_fs = 'gridfile_fs.mat';
save(gridfile_fs, 'composition_matrix', 'material_matrix', '-v7.3');


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
% Write the created input file to the cylinder_input.mat file.
iteratefdtd_matrix(input_filename,'filesetup','cylinder_input.mat',gridfile_cyl,illfile);
% Create the input to the simulation that has no obstacle, and only models
% free space.
% Write the created input file to the freespace_input.mat file.
iteratefdtd_matrix(input_filename,'filesetup','freespace_input.mat',gridfile_fs,illfile);

%% Run the tdms executable
% These commands can be carried out via the CLI
% Because MATLAB ships with versions of libstdc++ and prepends these to the
% path when calling system, when running tdms from MATLAB it is necessary
% to overwrite the search path for the libraries we need.
% TDMS is built in such a way that it does not rely on environment
% variables being present to run.
% ON UNIX: You require the 'LD_LIBRARY_PATH', '' option
% ON MAC: You require the 'DYLD_LIBRARY_PATH', '' option
% ON WINDOWS: You should not require any additional arguments

% Run the simulation for the freespace setup
% On the command line, in the examples/arc_01 folder, run
% tdms freespace_input.mat freespace_output.mat
system('tdms freespace_input.mat freespace_output.mat', 'LD_LIBRARY_PATH', '', 'DYLD_LIBRARY_PATH', '');
% Run the simulation for the cylindrical scatterer
% On the command line, in the examples/arc_01 folder, run
% tdms cylinder_input.mat cylinder_output.mat
system('tdms cylinder_input.mat cylinder_output.mat', 'LD_LIBRARY_PATH', '', 'DYLD_LIBRARY_PATH', '');

%% View the results
% Clear the variables that we have created during the setup and run phases
% above.
clear; close all;

% Load the output data from the files
data_cylinder = load('cylinder_output.mat');
data_freespace = load('freespace_output.mat');

% Plot the profiles of the incident beams in freespace and when the
% scattering cylinder is present
beam_figure = figure(1);

subplot(2,1,1);
imagesc(data_freespace.x_i,data_freespace.z_i,abs(squeeze(data_freespace.Ex_i)));
axis square;
title('Focussed beam ($E_{x}$) in free space', 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$z$', 'Interpreter', 'latex');

subplot(2,1,2);
imagesc(data_cylinder.x_i,data_cylinder.z_i,abs(squeeze(data_cylinder.Ex_i)));
axis square;
title('Focussed beam ($E_{x}$) in with scattering cylinder', 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$z$', 'Interpreter', 'latex');

% Plot a familiar TDMS image
tdms_figure = figure(2);
imagesc(data_cylinder.x_i, data_cylinder.z_i, abs(squeeze(data_cylinder.Hz_out)));
axis square;
title('Normalised $H_{z}$ component of the scattered field from a cylinder', 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$z$', 'Interpreter', 'latex');

exit;
