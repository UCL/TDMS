% Define the tests to be all locally defined functions
function tests = test_iteratefdtd_matrix_source_configurations
    tests = functiontests(localfunctions);
end

function testIteratefdtdMatrixSourceCombinations(testcase)
    %% The source field can be produced by iteratefdtd_matrix in a number of ways, but is governed by a combination of 5 variables that come from either an illumination file, or the input file itself:
    %%  a time-domain field (exi, eyi), using FDTD, compactsource, efname, hfname.
    % Certain combinations of these inputs are contradictory, whilst other combinations are not sufficient to create the source data.
    % This file tests that combinations that we expect to fail raise errors, and those that we don't exit without errors. This unit test was initially named arc_19, however since it does not contain any TDMS calls it has been absorbed into the unit test framework.

    case_errors = zeros(2, 2, 2, 2, 2);
    % The 2*2*2*2*2 array case_errors is indexed as (TD-field, using FDTD, compactsource, efname, hfname), using index 1 -> that variable takes the value logical 0 and 2->1 in similar fashion. Entries are also boolean flags: zeros error should be raised for the particular combination of (TD-field, using FDTD, compactsource, efname, hfname) that makes up the index. 1s obviously indicate the converse - an IncompatibleSourceInput error should be raised.
    case_errors(2, 2, 2, 2, 2) = 1;
    case_errors(2, 2, 2, 1, 2) = 1;
    case_errors(2, 2, 1, 2, 1) = 1;
    case_errors(2, 2, 1, 1, 2) = 1;
    case_errors(2, 1, 2, 2, 2) = 1;
    case_errors(2, 1, 2, 1, 2) = 1;
    case_errors(2, 1, 1, 2, 2) = 1;
    case_errors(2, 1, 1, 2, 1) = 1;
    case_errors(2, 1, 1, 1, 2) = 1;
    case_errors(2, 1, 1, 1, 1) = 1;
    case_errors(1, 2, 2, 2, 2) = 1;
    case_errors(1, 2, 2, 1, 2) = 1;
    case_errors(1, 2, 2, 1, 1) = 1;
    case_errors(1, 2, 1, 2, 1) = 1;
    case_errors(1, 2, 1, 1, 2) = 1;
    case_errors(1, 2, 1, 1, 1) = 1;
    case_errors(1, 1, 2, 2, 2) = 1;
    case_errors(1, 1, 2, 1, 2) = 1;
    case_errors(1, 1, 2, 1, 1) = 1;
    case_errors(1, 1, 1, 2, 2) = 1;
    case_errors(1, 1, 1, 2, 1) = 1;
    case_errors(1, 1, 1, 1, 2) = 1;
    case_errors(1, 1, 1, 1, 1) = 1;

    for td_field=0:1
        for using_fdtd=0:1
            for compactsource=0:1
                for efname_populated=0:1
                    for hfname_populated=0:1
                        % If we are expecting an error, run as such
                        if case_errors(td_field+1, ...
                                        using_fdtd+1, ...
                                        compactsource+1, ...
                                        efname_populated+1, ...
                                        hfname_populated+1)
                            error_type = 'TDMSException:IncompatibleSourceInput';
                        else
                            error_type = '';
                        end
                        % Run the test case, expecting either the error above or success
                        fprintf('%d %d (%d) %d %d %d\n',td_field, ...
                                                            using_fdtd, ~using_fdtd, ...
                                                            compactsource, ...
                                                            efname_populated, ...
                                                            hfname_populated);
                        run_in_temporary_directory(testcase, ...
                            @() call_to_iteratefdtd_matrix('pstd_input_file_2D.m', ...
                                                            td_field, ...
                                                            using_fdtd, ...
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

function call_to_iteratefdtd_matrix(filename, td_field, using_fdtd, compactsource, efname_nonempty, hfname_nonempty)
    gridfile = 'gridfile.mat';
    illfile = '';

    create_gridfile(gridfile);
    if td_field
        illfile = 'illfile.mat';
        create_tdfield_data(illfile);
    end
    set_use_pstd(filename, ~using_fdtd); % The input file flag is whether to swap off the FDTD method, hence negation
    set_compactsource(filename, compactsource);
    set_efname(filename, efname_nonempty);
    set_hfname(filename, hfname_nonempty);

    iteratefdtd_matrix(filename, 'filesetup', 'tmp_input', gridfile, illfile);
end
