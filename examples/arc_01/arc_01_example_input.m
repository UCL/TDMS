%% arc_01_example_input.m
% An example input file that can be passed to iteratefdtd_matrix to create
% an input file, which in turn may be passed to the tdms executable to
% perform the simulation.
% Section numbers (3.2.X) correspond to those provided in the PDF
% documentation.

%% Characteristic quantities
% These values are not read by iteratefdtd_matrix, but define
% characteristic scales for various quantities in the simulation.

% A characteristic wavelength from which we will derive other length-related input quantities,
% to pass to iteratefdtd_matrix and thus to tdms
lambda = 1300e-9;


%% Input flags to tdms executable

% Whether to use bandlimited interpolation over cubic interpolation
use_bli = 0;
% Whether to use the PSTD method over the FDTD method
use_pstd = 1;

%% 3.2.1 Grid

% Specifies the dimensionality of the simulation. IE whether we are
% simulating a full electromagnetic field ('3'), or just the 'TE' or 'TM'
% modes.
dimension = '3';

% Define the computational grid size, by giving the number of Yee cells along each coordinate axis.
% With J = 0, we define a 2D simulation rather than a 3D simulation.
I = 256;
J = 0;
K = 256;

% An array of length N_l, where N_l is the number of layers in the
% multilayer structure to be simulated.
% The i-th element is the z-index of the Yee cell at which an interface
% between the i-th and (i+1)-th layer occurs.
% If set to an empty vector [], the medium is assumed homogeneous.
multilayer = [];
% The remaining properties are either N_l-length vectors whose i-th element
% defines the material property of the i-th layer, or scalars in the event
% that we have a homogeneous medium.
epsr = 1.35^2;  % Relative permittivity of the layers
mur = 1;        % Relative permeability of the layers

% Size of the Yee cell in metres, in each of the coordinate directions.
delta = struct();
delta.x = lambda/4;
delta.y = lambda/4;
delta.z = lambda/4;

%% 3.2.4 FDTD specific
% This section appears first as we will be setting f_an in 3.2.2 based on
% the simulation timestep interval.

% Courant time step
dt = 2/sqrt(2)/pi*delta.x/(3e8/1.35)*.95;

% Number of time steps to perform
Nt = 500;

%% 3.2.2 Source

% Whether we are using a compact source condition,
% https://github.com/UCL/TDMS/issues/259
compactsource = 1;

% Frequency in Hertz of the incident EM field
f_an = asin( 2*pi/1300e-9*2.997924580105029e+08*dt/2)/(pi*dt);

% Define the planes where the incident waveforms are introduced.
% The variable interface has 6 members; I0, I1, J0, J1, K0, and K1, which are 1x2 vectors.
% The {I,J,K} indicates which coordinate direction {x,y,z} (respectively)
% the plane is perpendicular to.
% The {0,1} indicate which plane appears {earlier,later} in the computational grid.
% For each plane, the first entry is the constant {x,y,z} index of all Yee
% cells within the plane.
% The second entry is a boolean indicating whether or not an incident field is to be
% introduced at that particular plane.
interface.I0 = [5 0];   % Plane containing Yee cells w/ index of the form (5, j, k), at which NO incident field is applied.
interface.I1 = [I-5 0];
interface.J0 = [5 0];
interface.J1 = [J-5 0];
interface.K0 = [10 1];  % Plane containing Yee cells w/ index of the form (i, j, 10), at which an incident field is applied.
interface.K1 = [K-5 0];

% Function names (present in the MATLAB path) to use to generate the
% source-field on the incident planes.
% Whether one, both, or neither name should be provided depends on the
% manner in which you choose to define the source terms; through
% compactsource and a callable function, or a time-domain field in a pre-made .mat file.
% In this example, we call the efield_gauss_base function which creates a
% Gaussian E-field on the incident planes.
efname = 'efield_gauss_base';
hfname = '';

% The index of the Yee cell which fixes the origin of the Cartesian coordinate system of the grid.
% The incident field is defined relative to this coordinate system.
illorigin = [floor(I/2) floor(J/2) floor(K/2)];
% Sets the z-coordinate for the origin of the illorigin cell.
z_launch = 0;

% Spectral width of the modulating pulse of the incident field, in metres
wavelengthwidth = 120e-9;

%% 3.2.3 Simulation Type

% Specifies the type of source being used ('steadystate' or 'pulsed'). See section 3.2.3 of the
% documentation PDF.
sourcemode = 'pulsed';

% Defines the run mode of the simulation, if being run from MATLAB.
% If tdms is being run at the command-line, these values are ignored.
% 'analyse' : Sub-results can be saved using the statements in outputs_array
% 'complete': Only final results will be saved
runmode = 'complete';

%% 3.2.5 Output

% Not used as run mode is complete.
% If running TDMS from MATLAB, variables can be extracted from individual timesteps via these commands
outputs_array ={};

% Boolean indicating whether to extract phasors in the volume of the grid
exphasorsvolume = 1;

% Boolean indicating whether to extract phasors over a user-defined surface
exphasorssurface = 0;
% Specifies the user-defined surface to extract phasors over.
% phasorsurface has the form [I0 I1 J0 J1 K0 K1], which defines the
% extremes of a cuboid, whose surface will have phasors extracted over.
% These quantities are in the interior coordinate system.
phasorsurface = [5 I-5 1 1 20 K-5];

%% 3.2.6 Perfectly Matched Layer (PML)

% The order of the PML conductivity profile curve
n = 4;

% Maximum reflection at PML
R0 = 1e-7;

% Parameter mimicing the conductivity profile in the PML, applicable to dispersive materials only.
% See the documentation PDF, section 3.2.6.
kappa_max = 1;

% Number of PML cells in each direction.
Dxl = 10;
Dxu = 10;
Dyl = 0;
Dyu = 0;
Dzl = 10;
Dzu = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nlambda = 512;
lambda0 = 1300e-9;
dlambda = 170e-9;
b = 4*sqrt(log(2))*lambda0^2/(2*pi*3e8*dlambda);
omega0 = 2*pi*3e8/lambda0;

omega_min = omega0 - sqrt(4/b^2*log(10^3));
omega_max = omega0 + sqrt(4/b^2*log(10^3));

lambda_min = 3e8*2*pi/omega_max;
lambda_max = 3e8*2*pi/omega_min;


omega_vec = linspace(omega_min,omega_max,Nlambda);
k_vec = omega_vec/2.997924580105029e+08;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_ex_vec = asin( k_vec*2.997924580105029e+08*dt/2)/(pi*dt);

%ignore all below here
exdetintegral=0;
k_det_obs=10;
%k_obs = k_det_obs;
NA_det=(7e-3/2)/36e-3;
%NA = NA_det;
beta_det=25/36;
detmodevec=1:3;
detsensefun='gaussian_d_telesto_matlab';
air_interface = [];

det_trans_x = (-30:2:30)*1e-6;
det_trans_y = (-30:2:30)*1e-6;

illspecfun = '';
