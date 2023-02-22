function [] = run_bscan(test_directory, input_filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function generates the files used as input to the executeable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create directory into which to place the input files, if it doesn't exist already
dir_to_place_input_mats = test_directory;%strcat(test_directory,'/in')
if ~exist(dir_to_place_input_mats, 'dir')
    mkdir(dir_to_place_input_mats);
end

%% Generate the input file

%start by defining the coordinates of the computational grid
[x,y,z,lambda] = fdtd_bounds(input_filename);

%15 micron radius cylinder
rad = 15e-6;
%refractive index of cylinder
refind = 1.42;

%insert a cylinder at the origin
y = 0;
[X,Y,Z] = ndgrid(x,y,z);

%generate scattering matrix
I = zeros(size(X));
%set all Yee cells within the cylinder to have index of 1
I( (X.^2 + Z.^2) < rad^2 ) = 1;
I( (end-3):end,1,:) = 0;
I( :, 1, (end-3):end) = 0;
inds = find(I(:));
[ii,jj,kk] = ind2sub(size(I), inds);
composition_matrix = [ii jj kk ones(size(ii))];
material_matrix = [1 refind^2 1 0 0 0     0     0     0 0 0];

save('gridfile_cyl', 'composition_matrix', 'material_matrix');
%setup free space matrix and save
composition_matrix = [];
save('gridfile_fs', 'composition_matrix', 'material_matrix');

%generate tdms executable input files
iteratefdtd_matrix(input_filename,'filesetup',strcat(dir_to_place_input_mats,'/pstd_cyl_input'),'gridfile_cyl.mat','');
iteratefdtd_matrix(input_filename,'filesetup',strcat(dir_to_place_input_mats,'/pstd_fs_input'),'gridfile_fs.mat','');

end
