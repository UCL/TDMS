% Define the tests as all the locally defined functions
function tests = test_iteratefdtd_matrix_function
    tests = functiontests(localfunctions);
end

% These tests check cases where iteratefdtd_matrix recieves a valid combination of input arguments from the input file, but those values themselves are incorrect. This covers array shapes, field members not being present, and other such things that assume _existence_ of the variable itself.

function testInvalidIlluminationSourceIJKErrors(testCase)
    input_file = 'pstd_input_file_2D.m';
    run_in_temporary_directory(testCase, ...
    @()createInputWithInvalidIlluminationSourceIJ(input_file), ...
    'TDMSException:InvalidIlluminationFile');
end

function createInputWithInvalidIlluminationSourceIJ(input_file)

    gridfile = 'gridfile.mat';
    create_gridfile(gridfile);
    set_compactsource(input_file, 1);

    Isource = zeros(size(1));
    Jsource = zeros(size(1));
    % Illumnation file must also have a Ksource tensor
    save(sprintf('invalid_illumnation_file'), 'Isource', 'Jsource', 'Jsource');

    iteratefdtd_matrix(input_file,'filesetup',...
    'tmp_input',gridfile,'invalid_illumnation_file.mat');
end

function testInvalidIlluminationSourceExiEyiErrors(testCase)
    run_in_temporary_directory(testCase, ...
        @()createInputWithInvalidIlluminationSourceExi(), ...
        'TDMSException:InvalidIlluminationFile');
end

function testInvalidStateWithIlluminationAndEField(testCase)
    run_in_temporary_directory(testCase, ...
        @()createInputWithIlluminationAndEField(), ...
        'TDMSException:IncompatibleInput');
end

function testFileSetupValidIlluminationFile2DSource(testCase)
    run_in_temporary_directory(testCase, ...
            @()createInputWithValidIlluminationSource2D(), ...
            '');
end

function testFileSetupInvalidIlluminationFile2D(testCase)
    run_in_temporary_directory(testCase, ...
            @()createInputIlluminationSourceWithInvalidDimensions2D(), ...
            'TDMSException:InvalidIlluminationDimensions');
end

function testFileSetupValidIlluminationFile3D(testCase)
    run_in_temporary_directory(testCase, ...
            @()createInputWithValidIlluminationSource3D(), ...
            '');
end

function testFileSetupInvalidIlluminationFile3D(testCase)
    run_in_temporary_directory(testCase, ...
            @()createInputIlluminationSourceWithInvalidDimensions3D(), ...
            'TDMSException:InvalidIlluminationDimensions');
end

function testFileSetupInvalidIlluminationFile2DExi(testCase)
    run_in_temporary_directory(testCase, ...
            @()createInputIlluminationSourceWithInvalidExiDimensions2D(), ...
            'TDMSException:InvalidIlluminationDimensions');
end

function createInputWithValidIlluminationSource2D()
    gridfile = 'gridfile.mat';
    create_gridfile(gridfile);
    createIlluminaionFileFrom('pstd_input_file_2D.m');

    [~] = iteratefdtd_matrix('pstd_input_file_2D.m','filesetup','input_file',gridfile,'illfile.mat');
end

function createInputWithValidIlluminationSource3D()

    create_gridfile();
    createIlluminaionFileFrom('pstd_input_file_3D.m');

    [~] = iteratefdtd_matrix('pstd_input_file_3D.m','filesetup','input_file','gridfile.mat','illfile.mat');
end

function createInputIlluminationSourceWithInvalidDimensions2D()

    create_gridfile();
    createIlluminaionFileFrom('pstd_input_file_2D.m');

    % Cannot generate the input with an invalid dimension size
    replace_in_file('pstd_input_file_2D.m', 'I = 256;', 'I = 264;');
    [~] = iteratefdtd_matrix('pstd_input_file_2D.m','filesetup','input_file','gridfile.mat','illfile.mat');
end

function createInputIlluminationSourceWithInvalidDimensions3D()

    create_gridfile();
    createIlluminaionFileFrom('pstd_input_file_3D.m');

    % Cannot generate the input with an invalid dimension size
    replace_in_file('pstd_input_file_3D.m', 'I = 128;', 'I = 138;');
    [~] = iteratefdtd_matrix('pstd_input_file_3D.m','filesetup','input_file','gridfile.mat','illfile.mat');
end


function createInputWithInvalidIlluminationSourceExi()

    create_gridfile();
    exi = zeros(277, 277, 500);
    % Illumnation file must both exi and eyi
    save(sprintf('invalid_illumnation_file'), 'exi', 'exi');

    [~] = iteratefdtd_matrix('pstd_input_file_2D.m','filesetup',...
    'tmp_input','gridfile.mat','invalid_illumnation_file.mat');
end

function createInputWithIlluminationAndEField()

    create_gridfile();
    replace_in_file('pstd_input_file_2D.m', 'efname = ''''', 'efname = ''tmp''')

    % creates a valid illumnation file
    exi = zeros(277, 277, 500);  % I_tot + 1, J_tot + 1, Nt
    eyi = zeros(277, 277, 500);  % I_tot + 1, J_tot + 1, Nt
    save(sprintf('illumnation_file'), 'exi', 'eyi');

    [~] = iteratefdtd_matrix('pstd_input_file_2D.m','filesetup',...
    'tmp_input','gridfile.mat','illumnation_file.mat');
end

function createInputIlluminationSourceWithInvalidExiDimensions2D()

    create_gridfile();

    % creates an invalid illumnation file
    exi = zeros(1, 277, 500);  % I_tot + 1, J_tot + 1, Nt
    eyi = zeros(277, 277, 500);  % I_tot + 1, J_tot + 1, Nt
    save(sprintf('invalid_illumnation_file'), 'exi', 'eyi');

    [~] = iteratefdtd_matrix('pstd_input_file_2D.m','filesetup',...
    'tmp_input','gridfile.mat','invalid_illumnation_file.mat');
end

% Create a valid illumination file
function createIlluminaionFileFrom(filename)
    defineEfnameHfnameIn(filename);
    [~] = iteratefdtd_matrix(filename,'illsetup','illfile','gridfile.mat','');
    removeEfnameHfnameFrom(filename);
end

% Edit an input file to define efname and hfname
function defineEfnameHfnameIn(filename)
    replace_in_file(filename, 'efname = ''''', 'efname = ''efield_gauss_base''');
    replace_in_file(filename, 'hfname = ''''', 'hfname = ''hfield_focused_equiv''');
end

% Edit an input file to remove the definitions of efname and hfname
function removeEfnameHfnameFrom(filename)
    replace_in_file(filename, 'efname = ''efield_gauss_base''', 'efname = ''''');
    replace_in_file(filename, 'hfname = ''hfield_focused_equiv''', 'hfname = ''''');
end
