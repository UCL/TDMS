%specify the interpolation method (1 or not present - cubic, 2 - bandlimited)
intmethod = 1;

%these are not involved in the formal input file spec
lambda = 1300e-9;

%size of Yee cell in metres
delta.x = lambda/4;
delta.y = lambda/4;
delta.z = lambda/4;


%define the grid size, a square of side 1.5 wavelengths
I = 256;
J = 0;
K = 256;
%K = 5500;

%order of the PML conductivity profile curve
n = 4;

%maximum reflection at PML
R0 = 1e-7;

%number of PML cells in each direction
Dxl = 10;
Dxu = 10;
Dyl = 0;
Dyu = 0;
Dzl = 10;
Dzu = 10;

%courant time step
dt = 2/sqrt(2)/pi*delta.x/(3e8/1.35)*.95;

%define the number of time steps
Nt = 500;
%Nt=12000;

%water
epsr = [1.35^2];
mur = [1];
kappa_max = [1];
multilayer = [];

%frequency in Hz
f_an = asin( 2*pi/1300e-9*2.997924580105029e+08*dt/2)/(pi*dt);

%This is where we define the planes where the incident waveforms
%are introduced. The variable interface has 6 members, I0, I1, J0,
%J1, K0, K1 which are 1x2 vectors. The first entry is the position
%of the plane, in local coordinates and the second entry os whether
%or not to apply the interface condition at that plane
interface.I0 = [5 0];
interface.I1 = [I-5 0];
interface.J0 = [5 0];
interface.J1 = [J-5 0];
interface.K0 = [10 1];
interface.K1 = [K-5 0];

%not used as run mode is complete
outputs_array ={};

%these are the function names used to generate the field
efname = 'efield_gauss_base';
hfname = 'hfield_focused_equiv';

%this is the z value at which the field is launched, in metres
z_launch = 0;


%this defines the point about which the illumination is centred in
%the so called 'interior' coordinate system
illorigin = [floor(I/2) floor(J/2) floor(K/2)];

%the wavelength width (in m). This corresponds to the FWHM of the
wavelengthwidth = 120e-9;

%this defines the run mode of the simulation, can be 'analyse' or
%'complete'. 'analyse' means that sub results can be saved using
%the statements in outputs_array. When complete is specified, only
%the final results will be saved using the outputs_array statements.
%runmode = 'analyse';
runmode = 'complete';


%this is the kind of source mode, can be 'steadystate' or
sourcemode = 'pulsed';

%this determines whether or not to extract phasors in the volume of
%the grid
exphasorsvolume = 1;

%this determines whether or not to extract phasors around a
%specified surface
exphasorssurface = 0;

%this specifies a surface to extract the phasors at. These
%quantities are in interior coordinate system;
%has the form [I0 I1 J0 J1 K0 K1] which defines the extremes of a
%cuboid wihch defines the surface to extract phasors at
%These should be set so that the interpolation scheme can work
phasorsurface = [5 I-5 1 1 20 K-5];

%could be '3' 'TE' or 'TM'
dimension = '3';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nlambda = 512;
lambda0 = 1300e-9;
dlambda = 170e-9;
b = 4*sqrt(log(2))*lambda0^2/(2*pi*3e8*dlambda);
omega0 = 2*pi*3e8/lambda0;
%omega1 = 2*pi*3e8/(lambda0-dlambda/2);

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
