function set_hfname(filename, nonempty)
    %% Overwrites the value of hfname in the file provided.
    % If nonempty is 1, hfname will be populated (and set to 'hfield_focused_equiv'),
    % otherwise hfname will be set to the empty string.
    if nonempty
        replace_in_file(filename, 'hfname = ''''', 'hfname = ''hfield_focused_equiv''');
    else
        replace_in_file(filename, 'hfname = ''hfield_focused_equiv''', 'hfname = ''''');
    end
end
