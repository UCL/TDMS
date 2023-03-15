function [tf] = is_white_space(str)
    %% Returns 1 (true) if the str is entirely composed of whitespace and ends in a newline character.
    %% Usage is intended for when reading input files, and skipping over blank lines.
    % str   : String to check
    %
    % tf    : True/false result

% Logical array of indices where whitespace occurs
whitespace = isspace(str);
% If everything in the string except the whitespace is a single newline character at the end
if strcmp('\n', str(~whitespace)) && strcmp('\n', str(end-1:end))
    tf = 1;
else
    tf = 0;
end
end
