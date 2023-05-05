function run_in_temporary_directory(testCase, func, exception)
    %% Runs the testCase provided in a temporary directory, adding the necessary additional paths so that the unit tests can locate the matlab functions that need to be tested.
    % func      : Function handle that performs the workflow of the unit test, including throwing any expected exceptions.
    % exception : char array corresponding to the error type that is expected to be thrown. If empty, then we assume that func should not raise an error when called.

    %% Add necessary paths for unit tests to run
    [this_file_location, ~, ~] = fileparts(mfilename('fullpath'));
    path_to_iteratefdtdmatrix = strcat(this_file_location, '/../system/data/input_generation/matlab');
    path_to_utils = strcat(this_file_location, '/utils/');
    path_to_data_folder = strcat(this_file_location, '/data/');

    addpath(path_to_iteratefdtdmatrix, path_to_utils, path_to_data_folder);

    %% Create a temporary folder and run func in it, checking for exception.
    oldFolder = cd(createTemporaryFolder(testCase));
    if (strlength(exception) > 0)
        verifyError(testCase, func, exception);
    else
        func();
    end

    %% Return to the location we started in prior to the unit test
    cd(oldFolder);
end
