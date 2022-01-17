%function [res] = is_white_space(str)
%
%res is set to 1 if str is just a blank line (ie, composed of white space or just a new line).
%res is 0 otherwise.
function [res] = is_white_space(str)

res = 1;

counter = 1;
while counter < length(str)
    res = res & strncmp(str(counter),sprintf(' '),1);
    counter = counter + 1;
end
res = res & strncmp(str(counter),sprintf('\n'),1);
