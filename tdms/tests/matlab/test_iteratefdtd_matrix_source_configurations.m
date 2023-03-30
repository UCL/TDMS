% Define the tests to be all locally defined functions
function tests = test_iteratefdtd_matrix_source_configurations
    tests = functiontests(localfunctions);
end

%% The source field can be produced by iteratefdtd_matrix in a number of ways, but is governed by a combination of 5 variables that come from either an illumination file, or the input file itself:
%%  a time-domain field (exi, eyi), usecd, compactsource, efname, hfname.
% Certain combinations of these inputs are contradictory, whilst other combinations are not sufficient to create the source data.
% The allowable combinations (with 1 indicating that the variable is present and non-empty, 0 otherwise) are given in the iteratefdtd_matrix_allowable_source_configurations.md.
% This file tests that combinations that we expect to fail raise errors, and those that we don't exit without errors.

function testIteratefdtdMatrixSourceCombinations(testcase)
    % The 2*2*2*2*2 array below is indexed as (TD-field, usecd, compactsource, efname, hfname), using index 1 -> that variable takes the value logical 0 and 2->1 in similar fashion. Entries are strings: empty indicate no error should be raised in that case, otherwise the string indicates an error should be raised with the appropriate message.
    cases_to_errors = zeros(2, 2, 2, 2, 2);
    cases_to_errors(2, 2, 2, 2, 2) = 1;
    cases_to_errors(2, 2, 2, 1, 2) = 1;
    cases_to_errors(2, 2, 1, 2, 1) = 1;
    cases_to_errors(2, 2, 1, 1, 2) = 1;
    cases_to_errors(2, 1, 2, 2, 2) = 1;
    cases_to_errors(2, 1, 2, 1, 2) = 1;
    cases_to_errors(2, 1, 1, 2, 2) = 1;
    cases_to_errors(2, 1, 1, 2, 1) = 1;
    cases_to_errors(2, 1, 1, 1, 2) = 1;
    cases_to_errors(2, 1, 1, 1, 1) = 1;
    cases_to_errors(1, 2, 2, 2, 2) = 1;
    cases_to_errors(1, 2, 2, 1, 2) = 1;
    cases_to_errors(1, 2, 2, 1, 1) = 1;
    cases_to_errors(1, 2, 1, 2, 1) = 1;
    cases_to_errors(1, 2, 1, 1, 2) = 1;
    cases_to_errors(1, 2, 1, 1, 1) = 1;
    cases_to_errors(1, 1, 2, 2, 2) = 1;
    cases_to_errors(1, 1, 2, 1, 2) = 1;
    cases_to_errors(1, 1, 2, 1, 1) = 1;
    cases_to_errors(1, 1, 1, 2, 2) = 1;
    cases_to_errors(1, 1, 1, 2, 1) = 1;
    cases_to_errors(1, 1, 1, 1, 2) = 1;
    cases_to_errors(1, 1, 1, 1, 1) = 1;

    for td_field=0:1
        for usecd=0:1
            for compactsource=0:1
                for efname_populated=0:1
                    for hfname_populated=0:1
                        % If we are expecting an error, run as such
                        if cases_to_errors(td_field+1, ...
                                            usecd+1, ...
                                            compactsource+1, ...
                                            efname_populated+1, ...
                                            hfname_populated+1)
                            error_type = 'TDMSException:IncompatibleSourceInput';
                        else
                            error_type = '';
                        end
                        % Run the test case, expecting either the error above or success
                        runInTempoaryDirectory(testcase, ...
                            @() call_to_iteratefdtd_matrix('source_terms_input.m', ...
                                                            td_field, ...
                                                            usecd, ...
                                                            compactsource, ...
                                                            efname_populated, ...
                                                            hfname_populated), ...
                            error_type);
                    end
                end
            end
        end
    end
end

function runInTempoaryDirectory(testCase, func, exception)
    [this_file_location, ~, ~] = fileparts(mfilename('fullpath'));
    path_to_iteratefdtdmatrix = strcat(this_file_location, '/../system/data/input_generation/matlab');

    addpath(path_to_iteratefdtdmatrix, 'data/');

    oldFolder = cd(createTemporaryFolder(testCase));
    if (strlength(exception) > 0)
        verifyError(testCase, func, exception);
    else
        func();
    end
    cd(oldFolder);
end

function call_to_iteratefdtd_matrix(filename, td_field, usecd, compactsource, efname_nonempty, hfname_nonempty)
    gridfile = 'gridfile.mat';
    illfile = '';

    create_gridfile(gridfile);
    if td_field
        illfile = 'illfile.mat';
        create_tdfield_data(illfile);
    end
    set_usecd(filename, usecd);
    set_compactsource(filename, compactsource);
    set_efname(filename, efname_nonempty);
    set_hfname(filename, hfname_nonempty);

    iteratefdtd_matrix(filename, 'filesetup', 'tmp_input', gridfile, illfile);
end

function replaceInFile(filename, oldString, newString)
    %% Replace a particular string present in a file with another

    % Open original file and read all the lines
    file  = fopen(filename, 'r');
    lines = fread(file,'*char')';
    fclose(file);
    % Open the file again in write mode and overwrite line-by-line, replacing the oldString with the newString
    file  = fopen(filename,'w');
    lines = strrep(lines, oldString, newString);
    fprintf(file,'%s',lines);
    fclose(file);
end

function set_usecd(filename, value)
    %% Overwrites usecd with value, in the input file provided
    if value == 0
        replaceInFile(filename, 'usecd = 1', 'usecd = 0');
    elseif value == 1
        replaceInFile(filename, 'usecd = 0', 'usecd = 1');
    end
end

function set_compactsource(filename, value)
    %% Overwrites compactsource with value, in the input file provided
    if value == 0
        replaceInFile(filename, 'compactsource = 1', 'compactsource = 0');
    elseif value == 1
        replaceInFile(filename, 'compactsource = 0', 'compactsource = 1');
    end
end

function set_efname(filename, nonempty)
    %% Overwrites the value of efname in the file provided.
    % If nonempty is 1, efname will be populated (and set to 'efield_gauss_base'),
    % otherwise efname will be set to the empty string.
    if nonempty
        replaceInFile(filename, 'efname = ''''', 'efname = ''efield_gauss_base''');
    else
        replaceInFile(filename, 'efname = ''efield_gauss_base''', 'efname = ''''');
    end
end

function set_hfname(filename, nonempty)
    %% Overwrites the value of hfname in the file provided.
    % If nonempty is 1, hfname will be populated (and set to 'hfield_focused_equiv'),
    % otherwise hfname will be set to the empty string.
    if nonempty
        replaceInFile(filename, 'hfname = ''''', 'hfname = ''hfield_focused_equiv''');
    else
        replaceInFile(filename, 'hfname = ''hfield_focused_equiv''', 'hfname = ''''');
    end
end

function create_tdfield_data(output_name, I_tot, J_tot, Nt)
    %% Creates valid time-domain field data and saves it to a .mat folder.
    % The input file we are using as a basis, pstd_input_file_2D.m, will accept td-fields of dimensions
    % 277-277-500 as valid inputs, which are the defaults set here.
    if ~exist('I_tot', 'var')
        I_tot = 256;
    end
    if ~exist('J_tot', 'var')
        J_tot = 256;
    end
    if ~exist('Nt', 'var')
        Nt = 499;
    end
    exi = zeros(I_tot + 1, J_tot + 1, Nt);  % I_tot + 1, J_tot + 1, Nt
    eyi = zeros(I_tot + 1, J_tot + 1, Nt);  % I_tot + 1, J_tot + 1, Nt
    save(output_name, 'exi', 'eyi');
end

function create_gridfile(output_name)
    %% Create a valid gridfile to pass into iteratefdtd_matrix
    composition_matrix = [];
    material_matrix = [1 0 1 0 0 0 0 0 0 0 0];
    save(output_name, 'composition_matrix', 'material_matrix');
end
