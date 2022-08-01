%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section generates the files used as input to the executeable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%start by defining the coordinates of the computational grid
addpath('../../tdms/matlab');
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

%check whether the directories gridfiles/, in/, and out/ exist before we
%attempt to save to them.
%if they do not exist, create them
if ~exist('gridfiles', 'dir')
    mkdir gridfiles;
end %if (directory gridfiles/ exists)
if ~exist('in', 'dir')
    mkdir in;
end %if (directory in/ exists
if ~exist('out', 'dir')
    mkdir out;
end %if (directory out/ exists

save(sprintf('gridfiles/gridfile_cyl'), 'composition_matrix', 'material_matrix');
%setup free space matrix and save to directory gridfiles/
composition_matrix = [];
save(sprintf('gridfiles/gridfile_fs'), 'composition_matrix', 'material_matrix');

%generate C input files and save to directory in/
[fdtdgrid, Exs,Eys,Ezs,Hxs,Hys,Hzs,grid_labels,camplitudes,vertices,facets,Id] = iteratefdtd_matrix('pstd_input_file.m','filesetup','in/pstd_cyl','gridfiles/gridfile_cyl.mat','');
[fdtdgrid, Exs,Eys,Ezs,Hxs,Hys,Hzs,grid_labels,camplitudes,vertices,facets,Id] = iteratefdtd_matrix('pstd_input_file.m','filesetup','in/pstd_fs','gridfiles/gridfile_fs.mat','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now run the executeables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval('!tdms in/pstd_fs.mat out/pstd_fs.mat');
eval('!tdms in/pstd_cyl.mat out/pstd_cyl.mat');

%plot the data
dat_cyl = load('out/pstd_cyl');
dat_fs = load('out/pstd_fs');

figure(1);clf;
subplot(2,1,1);
imagesc(dat_fs.z_i,dat_fs.x_i,abs(squeeze(dat_fs.Ex_i)));
axis equal;
title('Focussed beam in free space');

subplot(2,1,2);
imagesc(dat_fs.z_i,dat_fs.x_i,abs(squeeze(dat_cyl.Ex_i)));
title('Focussed beam with scattering cylinder');
axis equal;
