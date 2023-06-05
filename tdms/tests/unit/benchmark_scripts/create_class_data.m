%% This script creates the class_data.mat file for use in the hdf5 tests for reading in TDMS objects.
close all;
clear;

%% interface
% interface is a struct with fields {I, J, K} {0, 1}, corresponding to the planes at which a source field is introduced.
% Each field is a 1 - by - 2 array of doubles; the first element being the
% {I, J, K} index of the Yee cell in which the particular plane lies.
% The second element is cast to a bool and indicates whether any boundary
% conditions are to be applied on said plane.
interface = struct();
interface.I0 = [1 0]; % I0 in the I = 1 plane, no source condition
interface.I1 = [4 1]; % I1 in the I = 4 plane, source condition applied
interface.J0 = [2 0]; % J0 in the J = 2 plane, no source condition
interface.J1 = [5 0]; % J1 in the J = 5 plane, no source condition
interface.K0 = [3 1]; % K0 in the K = 3 plane, source condition applied
interface.K1 = [6 1]; %K1 in the K=6 plane, source condition applied

%% f_vec & f_vec_bad
% f_vec is a struct with two fields, fx_vec and fy_vec.
% These are each 1D vectors of doubles, usually they have the same length
% however for our testing purposes we will make them different lengths.
f_vec = struct();
f_vec.fx_vec = [0.25 0.5 0.75 1.];          % Row vector w/ 4 elements
f_vec.fy_vec = [-0.25; -0.5; -0.75; -1.];   % Col vector w/ 4 elements

% f_vec_bad will be used as the failure case in our tests.
f_vec_bad = struct();
f_vec_bad.fx_vec = zeros(2,2);      % Non-2D elements are not expected!
f_vec_bad.fy_vec = zeros(2,1);      % This is permitted, but the former field should error

%% phasorsurface
% This array is read into the Cuboid class. It is just an array of 6 integers (stored as doubles OFC) that correspond to Yee cell indices in the various axial directions
phasorsurface = [1 4 2 5 3 6];

%% dispersive_aux
% Structure array that populates a DispersiveMultiLayer object.
% Has 9 fields, which are all vectors (theoretically of the same length but for testing purposes nope):
% alpha, beta, gamma, kappa_x, kappa_y, kappa_z, sigma_x, sigma_y, sigma_z
dispersive_aux = struct();
dispersive_aux.alpha = 0:9;
dispersive_aux.beta = 0:9;
dispersive_aux.gamma = 0:9;
dispersive_aux.kappa_x = 0:9;
dispersive_aux.kappa_y = 0:9;
dispersive_aux.kappa_z = 0:9;
dispersive_aux.sigma_x = 0:9;
dispersive_aux.sigma_y = 0:9;
dispersive_aux.sigma_z = 0:9;

%% CMaterial, DMaterial, CCollection
% CMaterial has 9 elements, CCollection can have 6 or 9. DMaterial has 6.
% Fields are vectors, named {C,D}{a,b,c}{x,y,z}. c-fields are missing if only 6 fields are present.
Cmaterial = struct();
Cmaterial.Cax = 0:4;
Cmaterial.Cbx = 0:4;
Cmaterial.Ccx = 0:4;
Cmaterial.Cay = 0:4;
Cmaterial.Cby = 0:4;
Cmaterial.Ccy = 0:4;
Cmaterial.Caz = 0:4;
Cmaterial.Cbz = 0:4;
Cmaterial.Ccz = 0:4;

C = Cmaterial;
C = rmfield(C, "Ccx", "Ccy", "Ccz");

Dmaterial = struct();
Dmaterial.Dax = 5:9;
Dmaterial.Dbx = 5:9;
Dmaterial.Dcx = 5:9;
Dmaterial.Day = 5:9;
Dmaterial.Dby = 5:9;
Dmaterial.Dcy = 5:9;

%% save variables to the file we need
% Save the files to the expected filename for the unit tests to read the
% data back in.

save("unit_test_data/class_data.mat", "-v7.3");
