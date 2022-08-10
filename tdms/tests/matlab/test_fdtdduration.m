
% Define the tests as all the locally defined functions
function tests = test_iteratefdtd_matrix_function
    tests = functiontests(localfunctions);
end


function testGaussianPulseParameters(testCase)
    addpath('../../matlab/', 'data/');

    [to_actual hwhm_actual] = fdtdduration('pstd_input_file_2D.m');

    pstd_input_file_2D;      % Run the input file to set variables
    epsilon = 1.35;          % Refractive index
    c = 299792458;           % speed of light in m/s

    lambda = c / (f_an * epsilon);

    % Half width at half maximum (HWHM)
    hwhm_expected = (lambda^2 ...
                     / ((c/epsilon) * wavelengthwidth)) ...
                     * 2 * sqrt(log(2)/pi);

    to_expected = hwhm_expected * sqrt(log(1e8)/pi);

    % Check that these are the expected values
    verifyEqual(testCase, hwhm_expected, hwhm_actual, "AbsTol", 1e-10);
    verifyEqual(testCase, to_expected, to_actual, "AbsTol", 1e-10);
end
