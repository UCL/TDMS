function set_compactsource(filename, value)
    %% Overwrites compactsource with value, in the input file provided
    if value == 0
        replace_in_file(filename, 'compactsource = 1', 'compactsource = 0');
    elseif value == 1
        replace_in_file(filename, 'compactsource = 0', 'compactsource = 1');
    end
end
