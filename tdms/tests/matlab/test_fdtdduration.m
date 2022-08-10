
% Define the tests as all the locally defined functions
function tests = test_iteratefdtd_matrix_function
    addpath('../../matlab/', 'data/');
    tests = functiontests(localfunctions);
end


function testGaussianPulseParameters(testCase)

    % ? and Half width at half maximum (HWHM)
    [to_actual hwhm_actual] = fdtdduration('pstd_input_file_2D.m');

    % Check that these are the expected values
    verifyEqual(testCase, 3.035046893364502e-14, hwhm_actual, "AbsTol", 1e-10);
    verifyEqual(testCase, 7.349254840150354e-14, to_actual, "AbsTol", 1e-10);
end


function testMinStepsFDTD(testCase)

    n_expected = 1483;
    % This also tests fdtdts, which computes the upper limite of the timestep
    n_acual = minsteps_fdtd('pstd_input_file_2D.m');

    verifyEqual(testCase, n_expected, n_acual);
end


function testMinStepsPSTD(testCase)

    n_expected = 891;
    n_acual = minsteps_pstd('pstd_input_file_2D.m');

    verifyEqual(testCase, n_expected, n_acual);
end
