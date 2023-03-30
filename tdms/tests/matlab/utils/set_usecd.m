function set_usecd(filename, value)
    %% Overwrites usecd with value, in the input file provided
    if value == 0
        replace_in_file(filename, 'usecd = 1', 'usecd = 0');
    elseif value == 1
        replace_in_file(filename, 'usecd = 0', 'usecd = 1');
    end
end
