function replace_in_file(filename, oldString, newString)
    %% Replace a particular string present in a file with another

    % Open original file and read all the lines
    file  = fopen(filename, 'r');
    lines = fread(file,'*char')';
    fclose(file);
    % Open the file again in write mode and overwrite line-by-line, replacing the oldString with the newString
    file  = fopen(filename,'w');
    lines = strrep(lines, oldString, newString);
    fprintf(file,'%s',lines);
    fclose(file);
end
