% Define the tests as all the locally defined functions
function tests = test_iteratefdtd_matrix
    tests = functiontests(localfunctions);
end

% These tests check cases where iteratefdtd_matrix recieves a valid combination of input arguments from the input file, but those values themselves are incorrect. This covers array shapes, field members not being present, and other such things that assume _existence_ of the variable itself.

function testIJKSource_not_defined(testCase)
    %% When one of Isource, Jsource or Ksource is not defined in a source file that should contain all three, an error should be raised.
    run_in_temporary_directory(testCase, ...
    @()create_source_missing_K('pstd_input_file_2D.m'), ...
    'TDMSException:InvalidIlluminationFile');
end

function create_source_missing_K(input_file)
    set_compactsource(input_file, 1);
    set_efname(input_file, 1);

    gridfile = 'gridfile.mat';
    create_gridfile(gridfile);

    Isource = zeros(size(1));
    Jsource = zeros(size(1));
    % Illumnation file must also have a Ksource tensor to be valid, deliberately don't include here.
    save('invalid_illumnation_file', 'Isource', 'Jsource', 'Jsource');

    iteratefdtd_matrix(input_file,'filesetup','tmp_input',gridfile,'invalid_illumnation_file.mat');
end

function testValid_inputs(testCase)
    %% Test that input files for 2D and 3D simulations that we expect to pass, do indeed pass.
    two_dimensions = 'pstd_input_file_2D.m';
    run_in_temporary_directory(testCase, @()create_valid_input(two_dimensions), '');

    three_dimensions = 'pstd_input_file_3D.m';
    run_in_temporary_directory(testCase, @()create_valid_input(three_dimensions), '');
end

function create_valid_input(input_file)
    set_compactsource(input_file, 1);
    set_usecd(input_file, 0);
    set_efname(input_file, 1);

    gridfile = 'gridfile.mat';
    illfile = 'illumination.mat';
    create_gridfile(gridfile);
    create_illumination_file(illfile, input_file, gridfile);

    iteratefdtd_matrix(input_file,'filesetup','tmp_input',gridfile,illfile);
end

function testInvalid_illumination_dimensions(testCase)
    %% Test that if the dimensionality of the grid does not match the illumination arrays provided, an error is raised.
    two_dimensions = 'pstd_input_file_2D.m';
    run_in_temporary_directory(testCase, ...
            @()create_souce_with_wrong_dims(two_dimensions, 2), ...
            'TDMSException:InvalidIlluminationDimensions');

    three_dimensions = 'pstd_input_file_3D.m';
    run_in_temporary_directory(testCase, ...
            @()create_souce_with_wrong_dims(three_dimensions, 3), ...
            'TDMSException:InvalidIlluminationDimensions');
end

function create_souce_with_wrong_dims(input_file, n_dimensions)
    set_compactsource(input_file, 1);
    set_usecd(input_file, 0);
    set_efname(input_file, 1);

    gridfile = 'gridfile.mat';
    create_gridfile(gridfile);
    illfile = 'illumination.mat';
    create_illumination_file(illfile, input_file, gridfile);

    % Change the grid size specified in the input file,
    % so the source dimensions are now invalid
    if n_dimensions == 2
        replace_in_file(input_file, 'I = 256;', 'I = 264;');
    elseif n_dimensions == 3
        replace_in_file(input_file, 'I = 128;', 'I = 138;');
    end
    iteratefdtd_matrix(input_file,'filesetup','tmp_input',gridfile,illfile);
end

function testInvalid_tdfield_dimensions(testCase)
    %% Test that passing in a pre-computed time-domain field with the incorrect dimensions throws an error.
    run_in_temporary_directory(testCase, ...
            @()create_tdfield_with_wrong_dims('pstd_input_file_2D.m'), ...
            'TDMSException:InvalidIlluminationDimensions');
end

function create_tdfield_with_wrong_dims(input_file)
    set_compactsource(input_file, 1);

    gridfile = 'gridfile.mat';
    create_gridfile(gridfile);
    tdfield = 'tdfield.mat';
    % In 2D, we are expecting TDfields of size 227 x 1 x 500.
    % This will create fields of size 1 x 1 x 500
    create_tdfield_data(tdfield, 1);

    iteratefdtd_matrix(input_file,'filesetup','tmp_input',gridfile, tdfield);
end
