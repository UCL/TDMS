function [n] = minsteps_pstd(input_file)
    %% Estimates the minmum number of steps required in a pstd simulation.
    %% Assumes the time step used is 0.95 * maximum time step, and that the trailing edge of the guassian pulse travels at the speed of light (dispersion free).
    %% Calculates the edge of the pulse travelling from the interface to the end wall and back.
    % input_file    : Config file to read values from
    %
    % n             : Minimum number of timesteps that will need to be performed in a fdtd simulation

%% Fetch the configuration information for this test

% Check that the input file can be found on the path
if isfile(input_file)
	% Run input file as a script to import variables into the workspace
	run(input_file);
else
	% Throw error - config file cannot be found
    error(sprintf('File %s could not be opened for reading',input_file));
end

% Check required variables have been set in the config file, and loaded
% Note that lambda will be returned once it is imported from the input file
required_variables = {'delta','K','interface','f_an','dt','epsr'};
n_required_variables = length(required_variables);
i = 1;
% Search through all required variables and check they exist in the workspace
% Throw error (and break loop early) if they are not found
while i<=n_required_variables
	if ~exist(required_variables{i}, 'var')
		error(sprintf('Required variable %s not present in %s', required_variables{i}, input_file));
		break
	end
	i = i + 1;
end

%% Load / compute parameters required
[~, ~, c] = import_constants;
[t0, hwhm] = fdtdduration(inputfile);

% Adjust so that we have the pulse close to 0 at the interface
t = 2*t0;

%% Compute timestep estimate
% This pulse has to propagate from the interface, to the other edge of the grid, and then out again
% This is the time that will be ellapsed
t = t + (K-interface.K1(1))*delta.z/(c/sqrt(epsr(1))) + K*delta.z/(c/sqrt(epsr(1)));
% Hence from the ellapsed time, we can compute the number of timesteps we need for the pulse to travel back & forth
n = ceil(t/(0.95*dt));
end
