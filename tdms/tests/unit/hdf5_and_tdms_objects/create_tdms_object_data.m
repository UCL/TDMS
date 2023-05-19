%% This script creates the class_data.mat file for use in the hdf5 tests for reading in TDMS objects.
close all;
clear;

% % interface % interface is a struct with fields {
  I, J, K
} {0, 1}, corresponding to the planes at which a source field is introduced.%
                          Each field is a 1 -
                  by - 2 array of doubles;
the first element being the %
        {I, J, K} index of the Yee cell in which the particular plane lies.The %
        second element is cast to a bool and indicates whether any boundary %
        conditions are to be applied on said plane.interface = struct();
interface.I0 = [1 0];
% I0 in the I = 1 plane, no source condition interface.I1 = [4 1];
% I1 in the I = 4 plane, source condition applied interface.J0 = [2 0];
% J0 in the J = 2 plane, no source condition interface.J1 = [5 0];
% J1 in the J = 5 plane, no source condition interface.K0 = [3 1];
% K0 in the K = 3 plane, source condition applied interface.K1 = [6 1]; %K1 in the K=6 plane, source condition applied

%% save variables to the file we need
% Save the files to the expected filename for the unit tests to read the
% data back in.

save("class_data.mat", "-v7.3");
