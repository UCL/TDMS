function set_compactsource(input_file, value)
    %% Overwrites compactsource with value, in the input_file provided
    if value == 0
        replace_in_file(input_file, 'compactsource = 1', 'compactsource = 0');
    elseif value == 1
        replace_in_file(input_file, 'compactsource = 0', 'compactsource = 1');
    end
end
