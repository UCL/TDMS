function [dt_upper] = fdtdts(input_file)
    %% Calculates the maximum allowable time step, subject to the stability criterion
    %% This corresponds to:
    %% dt_upper = ( c * ( dx^{-2} + dy^{-2} + dx^{-2} ) )^{-1}.
    % input_file    : Configuration file to read parameters from
    %
    % dt_upper      : Maximum allowable timestep

%% Fetch the configuration information for this test
delta = get_from_input_file(input_file, struct(), 'delta');

%% Define internal constants
[~, ~, c] = import_constants;

%% Compute maximum permissable timestep
dt_upper = 1/(c*sqrt(1/delta.x^2 + 1/delta.y^2 + 1/delta.z^2));
end
