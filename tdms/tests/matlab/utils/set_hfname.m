function set_hfname(input_file, nonempty)
    %% Overwrites the value of hfname in the input_file provided.
    % If nonempty is 1, hfname will be populated (and set to 'hfield_focused_equiv'),
    % otherwise hfname will be set to the empty string.
    if nonempty
        replace_in_file(input_file, 'hfname = ''''', 'hfname = ''hfield_focused_equiv''');
    else
        replace_in_file(input_file, 'hfname = ''hfield_focused_equiv''', 'hfname = ''''');
    end
end
