function [n] = minsteps_fdtd(input_file)
    %% Estimates the minmum number of steps required in a fdtd simulation.
    %% Assumes the time step used is 0.95 * maximum time step, and that the trailing edge of the guassian pulse travels at the speed of light (dispersion free).
    %% Calculates the edge of the pulse travelling from the interface to the end wall and back.
    % input_file    : Config file to read values from
    %
    % n             : Minimum number of timesteps that will need to be performed in a fdtd simulation

%% Fetch the configuration information for this test
[delta, K, interface, epsr] = get_from_input_file(input_file, struct(), ...
												'delta','K','interface','epsr');

%% Load / compute parameters required
[~, ~, c] = import_constants;
[t0, hwhm] = fdtdduration(input_file);
dt_upper = fdtdts(input_file);

% Adjust so that we have the pulse close to 0 at the interface
t = 2*t0;

%% Compute timestep estimate
% This pulse has to propagate from the interface, to the other edge of the grid, and then out again
% This is the time that will be ellapsed
t = t + (K-interface.K0(1))*delta.z/(c/sqrt(max(epsr))) + K*delta.z/(c/sqrt(max(epsr)));
% Hence from the ellapsed time, we can compute the number of timesteps we need for the pulse to travel back & forth
n = ceil(t/(0.95*dt_upper));
end
