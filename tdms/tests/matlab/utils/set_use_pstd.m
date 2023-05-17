function set_use_pstd(input_file, value)
    %% Overwrites use_pstd with value, in the input_file provided
    if value == 0
        replace_in_file(input_file, 'use_pstd = 1', 'use_pstd = 0');
    elseif value == 1
        replace_in_file(input_file, 'use_pstd = 0', 'use_pstd = 1');
    end
end
