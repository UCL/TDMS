function set_efname(input_file, nonempty)
    %% Overwrites the value of efname in the input_file provided.
    % If nonempty is 1, efname will be populated (and set to 'efield_gauss_base'),
    % otherwise efname will be set to the empty string.
    if nonempty
        replace_in_file(input_file, 'efname = ''''', 'efname = ''efield_gauss_base''');
    else
        replace_in_file(input_file, 'efname = ''efield_gauss_base''', 'efname = ''''');
    end
end
