function [ex_coords, ey_coords, tvec_E, fvec_E, f_an, hwhm, to_l] = getsourcecoords(input_file)
	%% Computes the spatial coordinates of the source field
	% input_file	: Config file to read parameters from
	%
	% ex_coords		: Spatial coordinates of the Ex component of the source field
	% ey_coords		: Spatial coordinates of the Ey component of the source field
	% tvec_E		: Time "coordinates" as a 1D array
	% fvec_E		: Frequency "coordinates" as a 1D array
	% f_an			:
	% hwhm			:
	% to_l			:

%% Fetch the configuration information for this test
% Optionals
optional_args = struct();
optional_args.z_launch = 0;
optional_args.multilayer = [];
% Fetch variables that we need
[delta, I, J, K, Dxl, ...
Dxu, Dyl, Dyu, Dzl, Dzu, ...
epsr, f_an, Nt, interface, wavelengthwidth, ...
illorigin, sourcemode, dt, z_launch, multilayer] = get_from_input_file(input_file, optional_args, ...
'delta', 'I', 'J', 'K', 'Dxl', ...
'Dxu', 'Dyl', 'Dyu', 'Dzl', 'Dzu', ...
'epsr', 'f_an', 'Nt', 'interface', 'wavelengthwidth', ...
'illorigin', 'sourcemode', 'dt', 'z_launch', 'multilayer');

%% Variable checks - interface
% Check that interface has been fully specified.
% If not, un-specified interface conditions will be set to 0 (IE no interface condition will be implemented)
interface_field = {'I0','I1','J0','J1','K0','K1'};
for lvar = 1:length(interface_field)
    if ~isfield(interface,interface_field{lvar})
	fprintf('interface.%s has not been defined, setting it to [0 0]\n', interface_field{lvar});
	interface.(interface_field{lvar}) = [0 0];
    end
end

% Check that the entries for I0...K1 are valid
if interface.I0(2) && interface.I1(2)
    if interface.I1(1) < interface.I0(1)
		error('Must have interface.I1>=interface.I0');
    end
    if interface.I1(1) > I
		error('Must have interface.I1 <= I');
    end
    if interface.I0(1) < 1
		error('Must have interface.I0 >= 1');
    end
end
if interface.J0(2) && interface.J1(2)
    if interface.J1(1) < interface.J0(1)
		error('Must have interface.J1>=interface.J0');
    end
    if interface.J1(1) > J
		error('Must have interface.J1 <= J');
    end
    if interface.J0(1) < 1
		error('Must have interface.J0 >= 1');
    end
end
if interface.K0(2) && interface.K1(2)
    if interface.K1(1) < interface.K0(1)
		error('Must have interface.K1>=interface.K0');
    end
    if interface.K1(1) > K
		error('Must have interface.K1 <= K');
    end
    if interface.K0(1) < 1
		error('Must have interface.K0 >= 1');
    end
end

%% Variable checks - other variables
% sourcemode
if ~(strncmp(sourcemode,'pulsed',6) | strncmp(sourcemode,'steadystate',11))
    error(sprintf('sourcemode ''%s'' invalid',sourcemode));
end

% multilayer and epsr
if isempty(multilayer)
    if numel(epsr)~=1
	error('epsr should have only a single element in the case of a single layer');
    end
else
    if numel(multilayer)~=(numel(epsr)-1)
	error('epsr should have one more member than multilayer');
    end
end

% Extent of the grid - including the PML
I_tot = I + Dxl + Dxu;
J_tot = J + Dyl + Dyu;
K_tot = K + Dzl + Dzu;

fprintf('Allocating grid...');
fdtdgrid = initialisesplitgrid(I,J,K,Dxl,Dxu,Dyl,Dyu,Dzl,Dzu);
fprintf('Done\n');

%% Impliment incident planewave
% We now implement an incident planewave using a total/scattered field formulation.
% For k<K1 we have scattered and for k>=K1 we have total. The field is coupled into system by means of update equations to maintain consistency.

if strcmp(sourcemode, 'pulsed')
	interface.I0(1) = 1;
    interface.J0(1) = 1;
    interface.I1(1) = I_tot+1;
    interface.J1(1) = J_tot+1;
    interface.I0(2) = 0;
    interface.I1(2) = 0;
    interface.J0(2) = 0;
    interface.J1(2) = 0;
    interface.K1(2) = 0;
    if (interface.K0(2)==0) & K~=0
		error('Running in pulsed mode with k0[0]=0, there is no point running');
    end
else
	% Adjust interface to use global coordinates
	interface.I0(1) = interface.I0(1) + Dxl;
	interface.I1(1) = interface.I1(1) + Dxl;
	interface.J0(1) = interface.J0(1) + Dyl;
	interface.J1(1) = interface.J1(1) + Dyl;
	interface.K0(1) = interface.K0(1) + Dzl;
	interface.K1(1) = interface.K1(1) + Dzl;
end
% Specify origin of the source
illorigin = illorigin + [Dxl Dyl Dzl];

% Set up the Ksource field.
% This has to be defined on a 2d array over the range (I0,I1)x(J0,J1). We calculate the field values assuming an origin for the illumination
i_source = (interface.I0(1):interface.I1(1)) - illorigin(1);
j_source = (interface.J0(1):interface.J1(1)) - illorigin(2);
k_source = interface.K0(1) - illorigin(3);

if interface.K0(2)
    %Ex, K0
    [x,y,z] = yeeposition(i_source, j_source, k_source, delta, 'Ex');
    z = z + z_launch;
    ex_coords.x = x;
    ex_coords.y = y;
    ex_coords.z = z;
    %Ey, K0
    [x,y,z] = yeeposition(i_source, j_source, k_source, delta, 'Ey');
    z = z + z_launch;
    ey_coords.x = x;
    ey_coords.y = y;
    ey_coords.z = z;
else
    ex_coords = [];
    ey_coords = [];
end

tvec_E = (1:Nt)*dt;

fvec_E = (0:(numel(tvec_E)-1))/numel(tvec_E)/diff(tvec_E(1:2));
fvec_E = fftshift(fvec_E);
ind0 = find(fvec_E==0);
fvec_E(1:(ind0-1)) = fvec_E(1:(ind0-1)) -fvec_E(ind0-1) - fvec_E(ind0+1);
fvec_E = ifftshift(fvec_E);

[~, ~, c] = import_constants;
refractive_index = sqrt(real(epsr(1)));
lambda_an = c/(f_an*refractive_index);
hwhm = lambda_an^2/((c/refractive_index)*wavelengthwidth)*2*sqrt(log(2)/pi);
to_l = hwhm*sqrt(log(1e8)/pi);
end
