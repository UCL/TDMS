%{
    BENCHMARK ERROR VALUES FOR BAND-LIMITED INTERPOLATION

Determines the error in MATLAB's BLi scheme for the test functions used in
the BLi unit tests in test_interpolation_functions.cpp:
- The constant function
- sin(2\pi x)
- A pulse function
- A complex-valued function of the form sin(2\pi x) + pulse(x)*i

These errors are printed to the screen upon execution of the script.
%}

clear;
close all;

nSamples = 100; % number of sample points to use
r = 2; % interpolate to midpoints of samples
N = 4; % use 2*4=8 sample points in interpolation
x = linspace(0,1,r*nSamples); % sample points

% Constant function interpolation
const_fn_data = const_fn(x);
const_fn_interp = interp(const_fn_data(1:r:end),r,N);
const_fn_err = max(abs( const_fn_interp(2:r:end) - const_fn_data(2:r:end) ));

% sin(2\pi x) interpolation
s2pi_data = s2pi(x);
s2pi_interp = interp(s2pi_data(1:r:end),r,N);
s2pi_err = max(abs( s2pi_interp(2:r:end) - s2pi_data(2:r:end) ));

% pulse function interpolation
pulse_data = pulse_fn(x);
pulse_interp = interp(pulse_data(1:r:end),r,N);
pulse_err = max(abs( pulse_interp(2:r:end) - pulse_data(2:r:end) ));

% complex: s2pi(x) + i*pulse_fn(x) interpolation
complex_data = complex_test(x);
complex_interp = interp(complex_data(1:r:end),r,N);
complex_err = max(abs( complex_interp(2:r:end) - complex_data(2:r:end) ));

fprintf("Constant function error: %.8e \n", const_fn_err);
fprintf("sin function error:      %.8e \n", s2pi_err);
fprintf("pulse function error:    %.8e \n", pulse_err);
fprintf("complex fn error:        %.8e \n", complex_err);

%% constant function
function [out] = const_fn(x)
    out = 1. + x*0.;
end

%% sin(2\pi x)
function [out] = s2pi(x)
    out = sin(2*pi*x);
end

%% compact pulse function
% f(x) = exp( -1 / (1-3|2x-1|)^2 )  when 3|2x-1|<1,
% f(x) = 0                          otherwise.
function [out] = pulse_fn(x)
    xhat = 3 * (2*x - 1);
    out = zeros(size(xhat));
    out(abs(xhat)>=1) = 0.;
    out(abs(xhat)<1) = exp( -1 ./ ( 1 - abs( xhat(abs(xhat)<1) ).^2 ));
end

%% complex-valued function
% Real part is s2pi(x)
% Imag part is pulse_fn(x)
function [out] = complex_test(x)
    out = s2pi(x) + pulse_fn(x)*1.i;
end