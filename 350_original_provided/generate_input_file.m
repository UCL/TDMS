% Adapted gen_input_file code so that we can run this on an arbitrary machine, though
% default values are setup for Will's laptop

function exitcode = generate_input_file(input_file_name, output_file_name, tdms_root)
if ~exist('input_file_name', 'var')
    input_file_name = 'input_file_ml.m';
end % if if ~exist('output_file_name', 'var') output_file_name = 'oct_ml_in';
end % if if ~exist('tdms_root',
                   'var') tdms_root = '/home/ccaegra/Documents/TDMS';
end % if

        % Add additional input generation paths matlab_supporting_fn_path =
        sprintf('%s/tdms/tests/system/data/input_generation/matlab', tdms_root);
addpath(matlab_supporting_fn_path);

composition_matrix = [];
material_matrix = [];
save scatfile_fs material_matrix composition_matrix;

iteratefdtd_matrix(input_file_name, 'filesetup', output_file_name,
                   'scatfile_fs.mat','');

% Remove additional paths rmpath(matlab_supporting_fn_path);

exitcode = 0;
end
