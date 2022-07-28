
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

function runInTempoaryDirectory(testCase, func, exception)
    addpath('../../matlab/', 'data/');

    oldFolder = cd(createTemporaryFolder(testCase));
    verifyError(testCase, func, exception);
    cd(oldFolder);
end

function createInputWithInvalidIlluminationSourceIJ()

    createGridFile();
    Isource = zeros(size(1));
    Jsource = zeros(size(1));
    % Illumnation file must also have a Ksource tensor
    save(sprintf('invalid_illumnation_file'), 'Isource', 'Jsource', 'Jsource');

    [~] = iteratefdtd_matrix('pstd_input_file.m','filesetup',...
    'tmp_input','gridfile.mat','invalid_illumnation_file.mat');
end

function createInputWithInvalidIlluminationSourceExi()

    createGridFile();
    exi = zeros(size(1));
    % Illumnation file must both exi and eyi
    save(sprintf('invalid_illumnation_file'), 'exi', 'exi');

    [~] = iteratefdtd_matrix('pstd_input_file.m','filesetup',...
    'tmp_input','gridfile.mat','invalid_illumnation_file.mat');
end

function createInputWithIlluminationAndEField()

    createGridFile();
    replaceInFile('pstd_input_file.m', 'efname = ''''', 'efname = ''tmp''')

    % creates a valid illumnation file
    exi = zeros(size(1));
    eyi = zeros(size(1));
    save(sprintf('illumnation_file'), 'exi', 'eyi');

    [~] = iteratefdtd_matrix('pstd_input_file.m','filesetup',...
    'tmp_input','gridfile.mat','illumnation_file.mat');
end

function createGridFile()

    lambda = 1300e-9;
    illorigin = [128 0 128]

    x = ((1:256) - illorigin(1))*lambda/4;
    z = ((1:256) - illorigin(3))*lambda/4;

    [X,Y,Z] = ndgrid(x,0,z);
    I = zeros(size(X));
    inds = find(I(:));
    [ii,jj,kk] = ind2sub(size(I), inds);
    composition_matrix = [ii jj kk ones(size(ii))];
    material_matrix = [1 0 1 0 0 0     0     0     0 0 0];
    save(sprintf('gridfile'), 'composition_matrix', 'material_matrix');
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
