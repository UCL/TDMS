function create_illumination_file(illfile_name, filename)
    %% Create a valid illumination file, saving it to the illfile_name, using the default gridfile settings and input values provided in filename.

    gridfile = 'gridfile_for_illsetup.mat';
    create_gridfile(gridfile);

    iteratefdtd_matrix(filename,'illsetup',illfile_name,'gridfile.mat','');
end
