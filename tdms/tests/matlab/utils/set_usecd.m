function set_usecd(input_file, value)
    %% Overwrites usecd with value, in the input_file provided
    if value == 0
        replace_in_file(input_file, 'usecd = 1', 'usecd = 0');
    elseif value == 1
        replace_in_file(input_file, 'usecd = 0', 'usecd = 1');
    end
end
