function create_illumination_file(illfile_name, filename, gridfile)
    %% Create a valid illumination file, saving it to the illfile_name, using input values provided in filename.
    %% If no gridfile is provided, the default "testing" gridfile will be produced and used.
    if ~exist('gridfile', 'var')
        gridfile = 'gridfile_for_illsetup.mat';
    end
    create_gridfile(gridfile);

    iteratefdtd_matrix(filename,'illsetup',illfile_name,gridfile,'');
end
