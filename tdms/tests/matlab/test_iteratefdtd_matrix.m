
% Define the tests as all the locally defined functions
function tests = test_iteratefdtd_matrix_function
    tests = functiontests(localfunctions);
end


function testInvalidIlluminationSourceIJKErrors(testCase)
    runInTempoaryDirectory(testCase, ...
    @()createInputWithInvalidIlluminationSourceIJ(), ...
    'TDMSException:InvalidIlluminationFile');
end

function testInvalidIlluminationSourceExiEyiErrors(testCase)
    runInTempoaryDirectory(testCase, ...
        @()createInputWithInvalidIlluminationSourceExi(), ...
        'TDMSException:InvalidIlluminationFile');
end

function testInvalidStateWithIlluminationAndEField(testCase)
    runInTempoaryDirectory(testCase, ...
        @()createInputWithIlluminationAndEField(), ...
        'TDMSException:IncompatibleInput');
end

function testFileSetupValidIlluminationFile2D(testCase)
    runInTempoaryDirectory(testCase, ...
            @()createInputWithValidIlluminationSource2D(), ...
            '');
end

function testFileSetupInvalidIlluminationFile2D(testCase)
    runInTempoaryDirectory(testCase, ...
            @()createInputIlluminationSourceWithInvalidDimensions2D(), ...
            'TDMSException:InvalidIlluminationFileDimensions');
end

function testFileSetupValidIlluminationFile3D(testCase)
    runInTempoaryDirectory(testCase, ...
            @()createInputWithValidIlluminationSource3D(), ...
            '');
end


function runInTempoaryDirectory(testCase, func, exception)
    addpath('../../matlab/', 'data/');

    oldFolder = cd(createTemporaryFolder(testCase));
    if (strlength(exception) > 0)
        verifyError(testCase, func, exception);
    else;
        func();
    end
    cd(oldFolder);
end

function createInputWithValidIlluminationSource2D()

    createGridFile();
    createIlluminaionFileFrom('pstd_input_file_2D.m');

    [~] = iteratefdtd_matrix('pstd_input_file_2D.m','filesetup','input_file','gridfile.mat','illfile.mat');
end

function createInputWithValidIlluminationSource3D()

    createGridFile();
    createIlluminaionFileFrom('pstd_input_file_3D.m');

    [~] = iteratefdtd_matrix('pstd_input_file_3D.m','filesetup','input_file','gridfile.mat','illfile.mat');
end

function createInputIlluminationSourceWithInvalidDimensions2D()

    createGridFile();
    createIlluminaionFileFrom('pstd_input_file_2D.m');

    % Cannot generate the input with an invalid dimension size
    replaceInFile('pstd_input_file_2D.m', 'I = 256;', 'I = 264;');
    [~] = iteratefdtd_matrix('pstd_input_file_2D.m','filesetup','input_file','gridfile.mat','illfile.mat');
end

function createInputWithInvalidIlluminationSourceIJ()

    createGridFile();
    Isource = zeros(size(1));
    Jsource = zeros(size(1));
    % Illumnation file must also have a Ksource tensor
    save(sprintf('invalid_illumnation_file'), 'Isource', 'Jsource', 'Jsource');

    [~] = iteratefdtd_matrix('pstd_input_file_2D.m','filesetup',...
    'tmp_input','gridfile.mat','invalid_illumnation_file.mat');
end

function createInputWithInvalidIlluminationSourceExi()

    createGridFile();
    exi = zeros(size(1));
    % Illumnation file must both exi and eyi
    save(sprintf('invalid_illumnation_file'), 'exi', 'exi');

    [~] = iteratefdtd_matrix('pstd_input_file_2D.m','filesetup',...
    'tmp_input','gridfile.mat','invalid_illumnation_file.mat');
end

function createInputWithIlluminationAndEField()

    createGridFile();
    replaceInFile('pstd_input_file_2D.m', 'efname = ''''', 'efname = ''tmp''')

    % creates a valid illumnation file
    exi = zeros(size(1));
    eyi = zeros(size(1));
    save(sprintf('illumnation_file'), 'exi', 'eyi');

    [~] = iteratefdtd_matrix('pstd_input_file_2D.m','filesetup',...
    'tmp_input','gridfile.mat','illumnation_file.mat');
end

% Create a valid grid file
function createGridFile()

    composition_matrix = [];
    material_matrix = [1 0 1 0 0 0 0 0 0 0 0];
    save(sprintf('gridfile'), 'composition_matrix', 'material_matrix');
end

% Create a valid illumination file
function createIlluminaionFileFrom(filename)
    defineEfnameHfnameIn(filename);
    [~] = iteratefdtd_matrix(filename,'illsetup','illfile','gridfile.mat','');
    removeEfnameHfnameFrom(filename);
end

function defineEfnameHfnameIn(filename)
    replaceInFile(filename, 'efname = ''''', 'efname = ''efield_gauss''');
    replaceInFile(filename, 'hfname = ''''', 'hfname = ''hfield_focused_equiv''');
end

function removeEfnameHfnameFrom(filename)
    replaceInFile(filename, 'efname = ''efield_gauss''', 'efname = ''''');
    replaceInFile(filename, 'hfname = ''hfield_focused_equiv''', 'hfname = ''''');
end

% Replace a string present in a file with another
function replaceInFile(filename, oldString, newString)

    file  = fopen(filename, 'r');
    lines = fread(file,'*char')';
    fclose(file);

    file  = fopen(filename,'w');
    lines = strrep(lines, oldString, newString);
    fprintf(file,'%s',lines);
    fclose(file);

end
