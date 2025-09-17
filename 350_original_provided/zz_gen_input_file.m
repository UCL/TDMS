% Original file as downloaded from the issue.
% Edited into generate_input_file.m to allow for testing on a different machine,
% and from the command-line.

tdms_root = '/usr/local/TDMS';


addpath(sprintf('%s/tdms/tests/system/data/input_generation/matlab',
                tdms_root));


composition_matrix = [];
material_matrix = [];
save scatfile_fs material_matrix composition_matrix;


iteratefdtd_matrix('input_file_ml.m', 'filesetup', 'oct_ml_in',
                   'scatfile_fs.mat','');

rmpath(sprintf('%s/tdms/tests/system/data/input_generation/matlab', tdms_root));
