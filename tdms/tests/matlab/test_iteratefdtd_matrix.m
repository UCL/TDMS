
function tests = test_iteratefdtd_matrix_function
tests = functiontests(localfunctions);
end

function testInvalidIlluminationSourceIJKErrors(testCase)
    verifyError(testCase, @()createInputWithInvalidIlluminationSourceIJ(), 'TDMSException:InvalidIlluminationFile');
    delete *.mat
end

function testInvalidIlluminationSourceExiEyiErrors(testCase)
    verifyError(testCase, @()createInputWithInvalidIlluminationSourceExi(), 'TDMSException:InvalidIlluminationFile');
    delete *.mat
end

function createInputWithInvalidIlluminationSourceIJ()

    addPathAndCreateGridFile();
    Isource = zeros(size(1));
    Jsource = zeros(size(1));
    % Illumnation file must also have a Ksource tensor
    save(sprintf('invalid_illumnation_file'), 'Isource', 'Jsource', 'Jsource');

    [~] = iteratefdtd_matrix('pstd_input_file.m','filesetup',...
    'tmp_input','gridfile_cyl.mat','invalid_illumnation_file.mat');
end

function createInputWithInvalidIlluminationSourceExi()

    addPathAndCreateGridFile();
    exi = zeros(size(1));
    % Illumnation file must both exi and eyi
    save(sprintf('invalid_illumnation_file'), 'exi', 'exi');

    [~] = iteratefdtd_matrix('pstd_input_file.m','filesetup',...
    'tmp_input','gridfile_cyl.mat','invalid_illumnation_file.mat');
end

function addPathAndCreateGridFile()
    addpath('../../matlab/', 'data/');

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
    save(sprintf('gridfile_cyl'), 'composition_matrix', 'material_matrix');
end
