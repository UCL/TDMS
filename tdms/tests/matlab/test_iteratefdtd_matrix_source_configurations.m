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
                        run_in_temporary_directory(testcase, ...
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
