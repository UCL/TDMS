%This file is basically the same as input_file_rw4 except that a
%wider pulse width has been specified

%these are not involved in the formal input file spec
lambda = 3e8/1e9;

%size of Yee cell in metres
delta.x = lambda/20; 
delta.y = lambda/20;
delta.z = lambda/20;


%define the grid size
I = 20;
J = 20;
K = 20;

%initialise grid

%order of the PML conductivity profile curve
n = 4;

%maximu reflection at PML
R0 = 1e-5;

%number of PML cells in each direction
Dxl = 10;
Dxu = 10;
Dyl = 10;
Dyu = 10;
Dzl = 10;
Dzu = 10;

%time stemp - subject to restriction
%dt = 2.85e-11;
dt = 2.5e-11;

%define the number of time steps
Nt = 5;

%define the dielectric propertries of the interior
epsr = 1;
mur = 1;

%frequency in Hz
f_an = 3e8/lambda; 

%This is where we define the planes where the incident waveforms
%are introduced. The variable interface has 6 members, I0, I1, J0,
%J1, K0, K1 which are 1x2 vectors. The first entry is the position
%of the plane, in local coordinates and the second entry os whether
%or not to apply the interface condition at that plane
interface.I0 = [5 0];
interface.I1 = [I-5 0];
interface.J0 = [5 0];
interface.J1 = [J-5 0];
interface.K0 = [5 1];
interface.K1 = [K-5 0];

%the file to save fieldrecord to
%fieldrecord_save = 'fieldrecord_equal_pml_noupdate';

%file containing the grid material composition
material_file = 'material_file_tdisp';

%specifies which field values are to be output
%outputs_array{1} = {'/tmp/output1','interior',[floor(I/2) floor(J/2) floor(K/2)],[floor(I/2) floor(J/2) floor(K/2)],'xyz','accumulate'};
outputs_array{1} = {'/tmp/output','interior',[10 10 8],[14 14 12],'xyz','accumulate'};
%outputs_array{2} = {'/tmp/fdtd_output3/2output','global',[1 1 K1+1],[I+Dxu+Dxl+1 J+Dyl+Dyu+1 K1+1],'xyz','dump'};
%outputs_array{3} = {'/tmp/fdtd_output3/3output','global',[1 1 K1+2],[I+Dxu+Dxl+1 J+Dyl+Dyu+1 K1+2],'xyz','dump'};
%outputs_array{4} = {'/tmp/fdtd_output3/4output','global',[1 1 K1+3],[I+Dxu+Dxl+1 J+Dyl+Dyu+1 K1+3],'xyz','dump'};
%outputs_array{5} = {'/tmp/fdtd_output3/5output','global',[1 1 K1+4],[I+Dxu+Dxl+1 J+Dyl+Dyu+1 K1+4],'xyz','dump'};
%outputs_array{6} = {'/tmp/fdtd_output3/6output','global',[60 1 1],[60 J+Dyl+Dyu+1 K+Dzl+Dzu+1],'xyz','dump'};
%outputs_array{7} = {'/tmp/fdtd_output3/7output','global',[1 60 1],[I+Dxu+Dxl+1 60 K+Dzl+Dzu+1],'xyz','dump'};
%outputs_array{1} = {'/tmp/fdtd_output3/8output','global',[floor(I/2) floor(J/2) 1],[floor(I/2) floor(J/2) K+Dzl+Dzu+1],'xyz','dump'};
%outputs_array = {};
    
%these are the function names used to generate the field
efname = 'efield_plane_oblique_lin';
hfname = 'hfield_plane_oblique_lin';

%this is the z value at which the field is launched, in metres
z_launch = 0;


%this defines the point about which the illumination is centred in
%the so called 'interior' coordinate system
illorigin = [floor(I/2)+2 floor(J/2)+2 interface.K0(1)];

%this determines whether or not to apply the equations for
%consistency at the lower boundary
lower_boundary_update = 'false';

%the wavelength width (in m). This corresponds to the FWHM of the
%wavelength spectrum
wavelengthwidth = 0.03;

%this defines the run mode of the simulation, can be 'analyse' or
%'complete'. 'analyse' means that sub results can be saved using
%the statements in outputs_array. When complete is specified, only
%the final results will be saved using the outputs_array statements.
%runmode = 'analyse';
runmode = 'complete';


%this is the kind of source mode, can be 'steadystate' or
%'pulsed'. 
sourcemode = 'pulsed';
%ourcemode = 'steadystate';

%this determines whether or not to extract phasors in the volume of
%the grid
exphasorsvolume = 1;

%this determines whether or not to extract phasors around a
%specified surface
exphasorssurface = 1;

%this specifies a surface to extract the phasors at. These
%quantities are in interior coordinate system;
%has the form [I0 I1 J0 J1 K0 K1] which defines the extremes of a
%cuboid wihch defines the surface to extract phasors at
%These should be set so that the interpolation scheme can work
phasorsurface = [6 I-6 6 J-6 6 K-6];
