function set_efname(filename, nonempty)
    %% Overwrites the value of efname in the file provided.
    % If nonempty is 1, efname will be populated (and set to 'efield_gauss_base'),
    % otherwise efname will be set to the empty string.
    if nonempty
        replace_in_file(filename, 'efname = ''''', 'efname = ''efield_gauss_base''');
    else
        replace_in_file(filename, 'efname = ''efield_gauss_base''', 'efname = ''''');
    end
end
