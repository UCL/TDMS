function [] = run_bscan(test_directory, input_filename, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function generates the files used as input to the executeable.
% It is an overarching hope that the system tests will be converted to Python after MATLAB header dependency removal from the tdms source code. In which case, this function will become part of the Python BScanArguments class, and will not need the long list of input arguments that it currently depends on.

% test_directory    : Path to directory into which to place generated input files.
% input_filename    : Path to the input file, defining run-specific dimensions, functions, etc.
% options           : 1-by-1 struct whose fieldnames contain the options for this regeneration. Details below:
% --------------------------------------------------------------
% Fieldname         | Description
% obstacle          | String, either 'fs', 'sph', 'cyl', or 'sc', defining the shape of the obstacle present.
% illsetup          | Bool, if False we do not need to call iteratefdtd_matrix in illsetup mode prior to filesetup mode. If True, we must do this, and the ill_filesetup field is populated.
% ill_filesetup     | String, not used if illsetup is false. If illsetup is True, points to the input file iteratefdtd_matrix needs to read in filesetup mode. input_filename is the file that needs to be read in illsetup mode.
% obstacle_radius   | Float, radius in microns of the obstacle (radius of the circular face for a cyl, radius of the sphere for sph).
% calc_tdfield      | Bool, whether calc_field_tdfield needs to be run prior to setting the inputs.
% refind            | Float, refractive index of the non-freespace portion of the medium
% output_name       | String, the name under which to save the .mat file containing the generated tdms input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create directory into which to place the input files, if it doesn't exist already
dir_to_place_input_mats = test_directory;
if ~exist(dir_to_place_input_mats, 'dir')
    mkdir(dir_to_place_input_mats);
end

%% Run BScan to generate input files

% Run calc_field_tdfield if we need to
if options.calc_tdfield
    tdfield_saved_to = calc_field_tdfield(input_filename, strcat(test_directory,'/eivars.mat'));
end

% Insert obstacle, typically at the origin, by introducing the composition matrix
composition_matrix = composition_matrix(input_filename, options.obstacle, options.obstacle_radius);
% Create the material matrix from the refractive index provided
material_matrix = [1 options.refind^2 1 0 0 0     0     0     0 0 0];

% Save matrices TSTK the freespace runs used to use the obstacle gridfiles though!?!?! Does this still work?!!?!?
gridfile = sprintf('gridfile_%s.mat',options.obstacle);
save(gridfile, 'composition_matrix', 'material_matrix');

%% Generate TDMS executable input files

% Work around illsetup and calc_tdfield being needed
if options.illsetup
    % Need to run in illsetup mode first
    illfile_mat = 'illfile.mat';
    iteratefdtd_matrix(input_filename, 'illsetup', illfile_mat, gridfile, '');
    % Now prepare for filesetup mode
    filesetup_input_file = options.ill_filesetup;
else
    if options.calc_field_tdfield
        % We don't need to call in illsetup mode, but do need to pass an illumination file
        illfile_mat = tdfield_saved_to;
    else
        illfile_mat = '';
    end
    % Pass in the original input file
    filesetup_input_file = input_filename;
end

% Call iteratefdtd_matrix in filesetup mode to generate input .mat files for tdms executable
iteratefdtd_matrix(filesetup_input_file, 'filesetup', options.output_name, gridfile, illfile_mat);
end
