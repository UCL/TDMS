%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section generates the files used as input to the executeable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%start by defining the coordinates of the computational grid
addpath('../../tdms/tests/system/data/input_generation/matlab');
[x,y,z,lambda] = fdtd_bounds('pstd_input_file.m');

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
I( find( (X.^2 + Z.^2) < rad^2 )) = 1;
I( (end-3):end,1,:) = 0;
I( :, 1, (end-3):end) = 0;
inds = find(I(:));
[ii,jj,kk] = ind2sub(size(I), inds);
composition_matrix = [ii jj kk ones(size(ii))];
material_matrix = [1 refind^2 1 0 0 0     0     0     0 0 0];

save('gridfile_cyl', 'composition_matrix', 'material_matrix');
%setup free space matrix and save to directory gridfiles/
composition_matrix = [];
save('gridfile_fs', 'composition_matrix', 'material_matrix');

%generate tdms executable input files
[fdtdgrid, Exs,Eys,Ezs,Hxs,Hys,Hzs,grid_labels,camplitudes,vertices,facets,Id] = iteratefdtd_matrix('pstd_input_file.m','filesetup','in_pstd_cyl','gridfile_cyl.mat','');
[fdtdgrid, Exs,Eys,Ezs,Hxs,Hys,Hzs,grid_labels,camplitudes,vertices,facets,Id] = iteratefdtd_matrix('pstd_input_file.m','filesetup','in_pstd_fs','gridfile_fs.mat','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now run the executeables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval('!tdms in_pstd_fs.mat out_pstd_fs.mat');
eval('!tdms in_pstd_cyl.mat out_pstd_cyl.mat');

%plot the data
dat_cyl = load('out_pstd_cyl');
dat_fs = load('out_pstd_fs');

figure(1);clf;
subplot(2,1,1);
imagesc(dat_fs.z_i,dat_fs.x_i,abs(squeeze(dat_fs.Ex_i)));
axis equal;
title('Focussed beam in free space');

subplot(2,1,2);
imagesc(dat_fs.z_i,dat_fs.x_i,abs(squeeze(dat_cyl.Ex_i)));
title('Focussed beam with scattering cylinder');
axis equal;

%%
figure;
imagesc(dat_cyl.z_i, dat_cyl.x_i, abs(squeeze(dat_cyl.Hz_out)));
axis equal;
