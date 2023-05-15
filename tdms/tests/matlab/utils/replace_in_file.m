function replace_in_file(filename, old_string, new_string)
    %% Replaces all occurances of old_string in filename with new_string.

    % Open original file and read all the lines
    file  = fopen(filename, 'r');
    lines = fread(file,'*char')';
    fclose(file);

    % Open the file again in write mode and overwrite line-by-line, replacing the old_string with the new_string
    file  = fopen(filename,'w');
    lines = strrep(lines, old_string, new_string);
    fprintf(file,'%s',lines);
    fclose(file);
end
