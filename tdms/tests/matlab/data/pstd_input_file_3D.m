% Whether to use bandlimited interpolation over cubic interpolation
use_bli = 0;
% Whether to use the PSTD method over the FDTD method
use_pstd = 0;
% Whether we are using a compact source condition
compactsource = 1;

lambda = 1300e-9;

delta.x = lambda/4;
delta.y = lambda/4;
delta.z = lambda/4;

I = 128;
J = 128;
K = 64;

n = 4;

R0 = 1e-7;

Dxl = 10;
Dxu = 10;
Dyl = 0;
Dyu = 0;
Dzl = 10;
Dzu = 10;

dt = 2/sqrt(3)/pi*delta.x/(3e8/1.35)*.95;

Nt = 500;

epsr = 1.35^2;
mur = 1;
kappa_max = 1;
multilayer = [];

f_an = asin( 2*pi/1300e-9*2.997924580105029e+08*dt/2)/(pi*dt);

interface.I0 = [5 0];
interface.I1 = [I-5 0];
interface.J0 = [5 0];
interface.J1 = [J-5 0];
interface.K0 = [10 1];
interface.K1 = [K-5 0];

outputs_array ={};

efname = '';
hfname = '';

z_launch = 0;

illorigin = [floor(I/2) floor(J/2) floor(K/2)];

wavelengthwidth = 120e-9;

runmode = 'complete';

sourcemode = 'pulsed';

exphasorsvolume = 1;

exphasorssurface = 1;

phasorsurface = [5 I-5 5 J-5 5 K-5];

dimension = '3';

fieldsample.i = (I/2-4):(I/2+4);
fieldsample.j = (J/2-4):(J/2+4);
fieldsample.k = (K/2-4):(K/2+4);
fieldsample.n = [2 4];

[ii,jj,kk] = ndgrid((I/2-4):(I/2+4),(J/2-4):(J/2+4),(K/2-4):(K/2+4));
campssample.vertices = [ii(:) jj(:) kk(:)];
campssample.components = [1 2 3];

Nlambda = 16;
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

f_ex_vec = asin( k_vec*2.997924580105029e+08*dt/2)/(pi*dt);

exdetintegral=0;
k_det_obs=10;
NA_det=(7e-3/2)/36e-3;
beta_det=25/36;
detmodevec=1:3;
detsensefun='gaussian_d_telesto_matlab';
air_interface = [];
det_trans_x = (-30:2:30)*1e-6;
det_trans_y = (-30:2:30)*1e-6;
illspecfun = '';
