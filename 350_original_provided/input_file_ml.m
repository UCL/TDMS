% PSTD simulation use_pstd = 1;

% Use linear interpolation use_bli = 0;

% Not a compact source condition compactsource = 1;

% these are not involved in the formal input file spec lambda = 1300e-9;

% size of Yee cell in metres delta.x = lambda / 6;
delta.y = lambda / 6;
delta.z = lambda / 6;

% define the grid size I = 20;
J = 20;
K = 40;

% initialise grid

        % order of the PML conductivity profile curve n = 4;

% maximum reflection at PML R0 = 1e-7;

% number of PML cells in each direction,
        this layer absorbs light % but does not reflect it Dxl = 10;
Dxu = 10;
Dyl = 10;
Dyu = 10;
Dzl = 10;
Dzu = 10;

% time step - subject to restriction dt =
        2 / sqrt(3) / pi * delta.x / (3e8 / 1.0) * .95;

% define the number of time steps Nt = 3000;

% define the dielectric propertries of the interior epsr = [1. ^ 2 1.1 ^ 2];
% epsr = (refractive index) ^ 2 mur = 1;
kappa_max = [1 1];
multilayer = [20];

% frequency in Hz f_an = 2.997924580105029e+08 / lambda;

interface.I0 = [5 0];
interface.I1 = [I - 5 0];
interface.J0 = [5 0];
interface.J1 = [J - 5 0];
interface.K0 = [10 1];
interface.K1 = [K - 5 0];

% not used as run mode is complete outputs_array = {};

% these are the function names used to generate the field efname =
        'efield_plane';
hfname = '';

% the wavelength width(in m).This corresponds to the FWHM of the
        wavelengthwidth = 100e-9;

% this defines the run mode of the simulation,
        can be 'analyse' or
                %
                        'complete'.'analyse' means that sub results can be
                                saved using %
                        the statements in outputs_array.When complete is
                                specified,
        only %
                the final results will be saved using the outputs_array
                        statements.%
                runmode = 'analyse';
runmode = 'complete';


% this is the kind of source mode,
        can be 'steadystate' or sourcemode = 'pulsed';

% this determines whether or
        not to extract phasors in the volume of % the grid exphasorsvolume = 1;

% this determines whether or not to extract phasors around a %
                                     specified surface exphasorssurface = 0;

% this specifies a surface to extract the phasors at.These %
        quantities are in interior coordinate system;
% has the form[I0 I1 J0 J1 K0 K1] which defines the extremes of a %
        cuboid wihch defines the surface to extract phasors at %
        These should be set so that the interpolation scheme can work
                phasorsurface = [5 I - 5 5 J - 5 5 K - 5];

% could be '3' 'TE' or 'TM' dimension = '3';

% this defines the point about which the illumination is centred in %
        the so called 'interior' coordinate system %
        this means the the illumination is actually focussed on a point 25 %
        cells in front of the interface.illorigin =
        [floor(I / 2) floor(J / 2) interface.K0(1)];

% this is the width of the pulse in terms of the time step % gpulsewidth = 400;

% k_vec =
        linspace(2 * pi / (lambda + 100e-9), 2 * pi / (lambda - 100e-9), 1801);
% k_vec = k_vec(2 : end);

% f_ex_vec = 2.997924580105029e+08 * k_vec / 2 / pi;
% f_ex_vec = 2.997924580105029e+08 / lambda;
% intphasorssurface = 0;

% k_vec =
        linspace(2 * pi / (lambda + 100e-9), 2 * pi / (lambda - 100e-9), 1801);
k_vec = linspace(2 * pi / 1420e-9, 2 * pi / 1180e-9, 2048);
dk = diff(k_vec(1 : 2));
k_vec = [((-100 : -1) * dk + k_vec(1)) k_vec((1 : 100) * dk + k_vec(end))];

f_ex_vec = 2.997924580105029e+08 * k_vec / 2 / pi;

exdetintegral = 0;
k_det_obs = 5;
NA_det = (7e-3 / 2) / 36e-3;
beta_det = 25 / 36;
detmodevec = 1;
detsensefun = 'gaussian_d_telesto';
