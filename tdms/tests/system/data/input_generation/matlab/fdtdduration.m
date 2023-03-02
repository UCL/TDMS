function [t0, hwhm] = fdtdduration(input_file)
    %% Calculates the parameters for the Gaussian pulse employed in the fdtd code according to:
    %% G(t) = \exp( -\pi(\frac{t-t_0}{hwhm})^2 ).
    %% hwhm is the half width at half maximum of the pulse.
    %% t_0 is the time delay. This delay ensures that when the pulse starts to be introduced, the envelope of the pulse has a magnitude of 10^{-8}. In turn, this ensures that there aren't any erroneous spectral components introduced as a result of having a sharp discontinuity in time.
    %
    % input_file : The configuration file to read parameters from
    % t0

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
required_variables = {'wavelengthwidth', 'f_an', 'epsr'};
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
refractive_index = sqrt(epsr(1));
omega_an         = 2*pi*f_an;
lambda_an        = c/(f_an*refractive_index);

%% Produce output values
hwhm = lambda_an^2/((c/refractive_index)*wavelengthwidth)*2*sqrt(log(2)/pi);
t0   = hwhm*sqrt(log(1e8)/pi);
end
