%function [to hwhm] = fdtdduration(inputfile)
%
%Calculates the parameters for the gaussian pulse employed in the
%fdtd code according to:
%
%G(t) = exp( -pi((t-to)/hwhm)^2 )
% where hwhm is the half width at half maximum of the pulse and to (t_0) is the
% time delay. The latter ensures that when the pulse starts to be introduced
% the envelope of the pulse has a magnitude of 10^(-8). This ensures that there
% aren't any erroneous spectral components introduced as a result of having a
% sharp discontinuity in time.
function [to, hwhm] = fdtdduration(inputfile)

[fid_input,message] = fopen(inputfile,'r');

%check if file was opened successfully
if fid_input== -1
    error(sprintf('File %s could not be opened for reading',inputfile));
end

%proceed to_l read in config information
current_line = fgets(fid_input);

while current_line ~= -1
    eval(current_line);
    current_line = fgets(fid_input);
end

[epso muo c] = import_constants;

if exist('wavelengthwidth') ~= 1
    error('wavelengthwidth is not defined - cannot determine hwhm or to');
end

if exist('f_an') ~= 1
    error('f_an is not defined - cannot determine hwhm or to');
end

refractive_index = sqrt(epsr(1));
omega_an         = 2*pi*f_an;
lambda_an        = c/(f_an*refractive_index);

hwhm = lambda_an^2/((c/refractive_index)*wavelengthwidth)*2*sqrt(log(2)/pi);
to   = hwhm*sqrt(log(1e8)/pi);
