function [] = run_bscan(test_directory, input_filename, non_fs_obstacle, illfile_extra_file, obstacle_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function generates the files used as input to the executeable.
% It is an overarching hope that the system tests will be converted to Python after MATLAB header dependency removal from the tdms source code. In which case, this function will become part of the Python BScanArguments class, and will not need the long list of input arguments that it currently depends on.

% test_directory : Path to directory into which to place generated input files
% input_filename : Path to the input file, defining run-specific dimensions, functions, etc
% non_fs_obstacle: String, either 'sph', 'cyl', defining the shape of the obstacle present in the non-free-space simulation
% illfile_extra_file: If present, we need to call iteratefdtd_matrix twice, once to setup the illumination and again to setup the .mat inputs. input_filename must be passed when iteratefdtd_matrix is in illsetup mode, and this file must be passed when it is in filesetup mode. If this variable contains an empty string, we simply need to pass input_filename to iteratefdtd_matrix in filesetup mode as usual.
% obstacle_radius: The radius in microns of the obstacle (radius of the circular face for a cyl, radius of the sphere for sph)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create directory into which to place the input files, if it doesn't exist already
dir_to_place_input_mats = test_directory;%strcat(test_directory,'/in')
if ~exist(dir_to_place_input_mats, 'dir')
    mkdir(dir_to_place_input_mats);
end

% Refractive index of obstacle
refind = 1.42;

% Insert obstacle, typically at the origin, by introducing the scattering matrix
I = scattering_matrix(input_filename, obstacle_radius, non_fs_obstacle);

% Generate additional matrices
inds = find(I(:));
[ii,jj,kk] = ind2sub(size(I), inds);
composition_matrix = [ii jj kk ones(size(ii))];
material_matrix = [1 refind^2 1 0 0 0     0     0     0 0 0];

% Save obstacle matrices
obstacle_gridfile = sprintf('gridfile_%s.mat',non_fs_obstacle);
save(obstacle_gridfile, 'composition_matrix', 'material_matrix');
% Setup & save freespace matrix
composition_matrix = [];
save('gridfile_fs', 'composition_matrix', 'material_matrix');

%% Generate TDMS executable input files

% This is the file that should be passed to iteratefdtd_matrix in filesetup mode
filesetup_input_file = '';
% This is the illfile that should be passed to iteratefdtd_matrix in filesetup mode. Empty string implies no illfile is needed, IE we didn't run in illsetup beforehand
illfile_produced = '';
% Setup illumination file if necessary
if strcmp(illfile_extra_file, "")
    % This is empty, so we just call iteratefdtd_matrix immediately using input_filename
    filesetup_input_file = input_filename;
else
    % We need a call to iteratefdtd_matrix in illsetup mode first
    illfile_produced = 'illfile';
    iteratefdtd_matrix(input_filename,'illsetup',illfile_produced,obstacle_gridfile,'');
    % Then call iteratefdtd_matrix in filesetup mode using illfile_extra_file as the input file
    filesetup_input_file = illfile_extra_file;
end

% Names to save .mat files under
obstacle_output_file = sprintf('%s/pstd_%s_input',dir_to_place_input_mats, non_fs_obstacle);
freespace_output_file = sprintf('%s/pstd_fs_input',dir_to_place_input_mats);

% Call iteratefdtd_matrix in filesetup mode to generate input .mat files for tdms executable
iteratefdtd_matrix(filesetup_input_file,'filesetup',obstacle_output_file,obstacle_gridfile,illfile_produced);
iteratefdtd_matrix(filesetup_input_file,'filesetup',freespace_output_file,'gridfile_fs.mat',illfile_produced);

end
