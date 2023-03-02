function [dt_upper] = fdtdts(input_file)
    %% Calculates the maximum allowable time step, subject to the stability criterion
    %% This corresponds to:
    %% dt_upper = ( c * ( dx^{-2} + dy^{-2} + dx^{-2} ) )^{-1}.
    % input_file    : Configuration file to read parameters from
    %
    % dt_upper      : Maximum allowable timestep

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
required_variables = {'delta'};
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

%% Define internal constants
[~, ~, c] = import_constants;

%% Compute maximum permissable timestep
dt_upper = 1/(c*sqrt(1/delta.x^2 + 1/delta.y^2 + 1/delta.z^2));
end
