function [t0, hwhm] = fdtdduration(input_file)
    %% Calculates the parameters for the Gaussian pulse employed in the fdtd code according to:
    %% G(t) = \exp( -\pi(\frac{t-t_0}{hwhm})^2 ).
    %% hwhm is the half width at half maximum of the pulse.
    %% t_0 is the time delay. This delay ensures that when the pulse starts to be introduced, the envelope of the pulse has a magnitude of 10^{-8}. In turn, this ensures that there aren't any erroneous spectral components introduced as a result of having a sharp discontinuity in time.
    %
    % input_file : The configuration file to read parameters from
    % t0

%% Fetch the configuration information for this test
[wavelengthwidth, f_an, epsr] = get_from_input_file(input_file, struct(), 'wavelengthwidth', 'f_an', 'epsr');

%% Define internal constants
[~, ~, c] = import_constants;
refractive_index = sqrt(epsr(1));
omega_an         = 2*pi*f_an;
lambda_an        = c/(f_an*refractive_index);

%% Produce output values
hwhm = lambda_an^2/((c/refractive_index)*wavelengthwidth)*2*sqrt(log(2)/pi);
t0   = hwhm*sqrt(log(1e8)/pi);
end
