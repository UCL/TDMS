% Define the tests as all the locally defined functions
function tests = test_iteratefdtd_matrix_function
    tests = functiontests(localfunctions);
end

% These tests check cases where iteratefdtd_matrix recieves a valid combination of input arguments from the input file, but those values themselves are incorrect. This covers array shapes, field members not being present, and other such things that assume _existence_ of the variable itself.

function testInvalidIlluminationSourceIJKErrors(testCase)
    input_file = 'pstd_input_file_2D.m';
    run_in_temporary_directory(testCase, ...
    @()createInputWithInvalidIJKSource(input_file), ...
    'TDMSException:InvalidIlluminationFile');
end

function createInputWithInvalidIJKSource(input_file)

    gridfile = 'gridfile.mat';
    create_gridfile(gridfile);
    set_compactsource(input_file, 1);
    set_efname(input_file, 1);

    Isource = zeros(size(1));
    Jsource = zeros(size(1));
    % Illumnation file must also have a Ksource tensor
    save(sprintf('invalid_illumnation_file'), 'Isource', 'Jsource', 'Jsource');

    iteratefdtd_matrix(input_file,'filesetup',...
    'tmp_input',gridfile,'invalid_illumnation_file.mat');
end

function testFileSetupValidIlluminationFile2DSource(testCase)
    run_in_temporary_directory(testCase, ...
            @()createInputWithValidIlluminationSource2D('pstd_input_file_2D.m'), ...
            '');
end

function createInputWithValidIlluminationSource2D(input_file)
    set_compactsource(input_file, 1);
    set_usecd(input_file, 0);
    set_efname(input_file, 1);

    gridfile = 'gridfile.mat';
    illfile = 'illumination.mat';
    create_gridfile(gridfile);
    create_illumination_file(illfile, input_file, gridfile);

    iteratefdtd_matrix(input_file,'filesetup','tmp_input',gridfile,illfile);
end

function testFileSetupInvalidIlluminationFile2D(testCase)
    run_in_temporary_directory(testCase, ...
            @()createInputIlluminationSourceWithInvalidDimensions2D('pstd_input_file_2D.m'), ...
            'TDMSException:InvalidIlluminationDimensions');
end

function createInputIlluminationSourceWithInvalidDimensions2D(input_file)
    set_compactsource(input_file, 1);
    set_usecd(input_file, 0);
    set_efname(input_file, 1);

    gridfile = 'gridfile.mat';
    create_gridfile(gridfile);
    illfile = 'illumination.mat';
    create_illumination_file(illfile, input_file, gridfile);

    % Cannot generate the input with an invalid dimension size
    replace_in_file(input_file, 'I = 256;', 'I = 264;');
    iteratefdtd_matrix(input_file,'filesetup','tmp_input',gridfile,illfile);
end

% refactor with the 2d thing above, since it's now the same test functionally just parameterised by the 2d/3d ness
function testFileSetupValidIlluminationFile3D(testCase)
    run_in_temporary_directory(testCase, ...
            @()createInputWithValidIlluminationSource3D('pstd_input_file_3D.m'), ...
            '');
end

function createInputWithValidIlluminationSource3D(input_file)
    set_compactsource(input_file, 1);
    set_usecd(input_file, 0);
    set_efname(input_file, 1);

    gridfile = 'gridfile.mat';
    create_gridfile(gridfile);
    illfile = 'illumination.mat';
    create_illumination_file(illfile, input_file, gridfile);

    iteratefdtd_matrix(input_file,'filesetup','tmp_input',gridfile, illfile);
end

% Refactor with it's 2d counterpart above (use same input file and just change the IJK values, since that's like the only difference!)
function testFileSetupInvalidIlluminationFile3D(testCase)
    run_in_temporary_directory(testCase, ...
            @()createInputIlluminationSourceWithInvalidDimensions3D('pstd_input_file_3D.m'), ...
            'TDMSException:InvalidIlluminationDimensions');
end

function createInputIlluminationSourceWithInvalidDimensions3D(input_file)
    set_compactsource(input_file, 1);
    set_usecd(input_file, 0);
    set_efname(input_file, 1);

    gridfile = 'gridfile.mat';
    create_gridfile(gridfile);
    illfile = 'illumination.mat';
    create_illumination_file(illfile, input_file, gridfile);

    % Cannot generate the input with an invalid dimension size
    replace_in_file(input_file, 'I = 128;', 'I = 138;');
    iteratefdtd_matrix(input_file,'filesetup','tmp_input',gridfile,illfile);
end

function testFileSetupInvalidIlluminationFile2DExi(testCase)
    run_in_temporary_directory(testCase, ...
            @()createInputIlluminationSourceWithInvalidExiDimensions2D('pstd_input_file_2D.m'), ...
            'TDMSException:InvalidIlluminationDimensions');
end

function createInputIlluminationSourceWithInvalidExiDimensions2D(input_file)
    set_compactsource(input_file, 1);

    gridfile = 'gridfile.mat';
    create_gridfile(gridfile);
    tdfield = 'tdfield.mat';
    % In 2D, we are expecting TDfields of size 227 x 1 x 500.
    % This will create fields of size 1 x 1 x 500
    create_tdfield_data(tdfield, 1);

    iteratefdtd_matrix(input_file,'filesetup','tmp_input',gridfile, tdfield);
end
