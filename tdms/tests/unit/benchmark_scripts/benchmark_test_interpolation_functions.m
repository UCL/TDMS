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
xi = linspace(0,1,nSamples);
xi5 = linspace(0,1,nSamples) + (0.5/(nSamples-1));
xi5 = xi5(1:end-1);

% Constant function interpolation
const_fn_data = const_fn(xi);
const_fn_exact = const_fn(xi5);
const_fn_interp = interp(const_fn_data,r,N);
const_fn_interp = const_fn_interp(2:r:end-1);
const_fn_ptwise_errors = const_fn_exact - const_fn_interp;
const_fn_max_ptwise_err = max(abs( const_fn_ptwise_errors ));
const_fn_norm_err = norm(const_fn_ptwise_errors);

% sin(2\pi x) interpolation
s2pi_data = s2pi(xi);
s2pi_exact = s2pi(xi5);
s2pi_interp = interp(s2pi_data,r,N);
s2pi_interp = s2pi_interp(2:r:end-1);
s2pi_ptwise_errors = s2pi_exact - s2pi_interp;
s2pi_max_ptwise_err = max(abs( s2pi_ptwise_errors ));
s2pi_norm_err = norm(s2pi_ptwise_errors);

% pulse function interpolation
pulse_data = pulse_fn(xi);
pulse_exact = pulse_fn(xi5);
pulse_interp = interp(pulse_data,r,N);
pulse_interp = pulse_interp(2:r:end-1);
pulse_ptwise_errors = pulse_exact - pulse_interp;
pulse_max_ptwise_err = max(abs( pulse_ptwise_errors ));
pulse_norm_err = norm(pulse_ptwise_errors);

% complex: s2pi(x) + i*pulse_fn(x) interpolation
complex_data = complex_test(xi);
complex_exact = complex_test(xi5);
complex_interp = interp(complex_data,r,N);
complex_interp = complex_interp(2:r:end-1);
complex_ptwise_errors = complex_interp - complex_exact;
complex_max_ptwise_err = max(abs( complex_ptwise_errors ));
complex_norm_err = norm(complex_ptwise_errors);

fprintf("Constant function max ptwise error: %.8e \n", const_fn_max_ptwise_err);
fprintf("sin function max ptwise error:      %.8e \n", s2pi_max_ptwise_err);
fprintf("pulse function max ptwise error:    %.8e \n", pulse_max_ptwise_err);
fprintf("complex fn max ptwise error:        %.8e \n", complex_max_ptwise_err);

fprintf("constant function norm error: %.8e \n", const_fn_norm_err);
fprintf("sin function norm error:      %.8e \n", s2pi_norm_err);
fprintf("pulse function norm error:    %.8e \n", pulse_norm_err);
fprintf("complex fn norm error:        %.8e \n", complex_norm_err);

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
