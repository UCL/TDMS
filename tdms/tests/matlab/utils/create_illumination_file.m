function create_illumination_file(illfile_name, input_file, gridfile)
    %% Create a valid illumination file, saving it to illfile_name.
    %% Input values provided in the input_file will be used to setup the illumination.
    %% You may optionally provide the name of the intermediate gridfile that will need to be produced, should you need to reference it later.
    if ~exist('gridfile', 'var')
        gridfile = 'gridfile_for_illsetup.mat';
    end
    create_gridfile(gridfile);

    iteratefdtd_matrix(input_file,'illsetup',illfile_name,gridfile,'');
end
