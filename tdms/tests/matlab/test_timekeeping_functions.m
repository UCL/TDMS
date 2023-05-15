% Tests are defined to be all locally-defined functions.
function tests = test_timekeeping_functions
    addpath('../system/data/input_generation/matlab')
    tests = functiontests(localfunctions);
end

function testGaussianPulseParameters(testCase)
    %% Test that fdtdduration successfully returns the expected Gaussian pulse parameters from the input file.
    % Parameters to be read are ? and the half-width at half-maximum (HWHM)

    [to_actual, hwhm_actual] = fdtdduration('data/pstd_input_file_2D.m');
    verifyEqual(testCase, 3.035046893364502e-14, hwhm_actual, "AbsTol", 1e-10);
    verifyEqual(testCase, 7.349254840150354e-14, to_actual, "AbsTol", 1e-10);
end

function testMinStepsFDTD(testCase)
    %% Verify that the minimum number of timesteps required by a FDTD simulation in the input file deduced by the minstep_fdtd function matches value returned by an analytic formula.
    % This also tests fdtdts, which computes the upper limit of the number of timesteps

    n_expected = 1483;
    n_acual = minsteps_fdtd('data/pstd_input_file_2D.m');
    verifyEqual(testCase, n_expected, n_acual);
end

function testMinStepsPSTD(testCase)
    %% Verify that the minimum number of timesteps required by a PSTD simulation in the input file deduced by the minstep_fdtd function matches value returned by an analytic formula.
    % This also tests fdtdts, which computes the upper limit of the number of timesteps

    n_expected = 891;
    n_acual = minsteps_pstd('data/pstd_input_file_2D.m');
    verifyEqual(testCase, n_expected, n_acual);
end
