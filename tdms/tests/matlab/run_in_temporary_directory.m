function run_in_temporary_directory(testCase, func, exception)
    [this_file_location, ~, ~] = fileparts(mfilename('fullpath'));
    path_to_iteratefdtdmatrix = strcat(this_file_location, '/../system/data/input_generation/matlab');
    path_to_utils = strcat(this_file_location, '/utils/');
    path_to_data_folder = strcat(this_file_location, '/data/');

    addpath(path_to_iteratefdtdmatrix, path_to_utils, path_to_data_folder);

    oldFolder = cd(createTemporaryFolder(testCase));
    if (strlength(exception) > 0)
        verifyError(testCase, func, exception);
    else
        func();
    end
    cd(oldFolder);
end
