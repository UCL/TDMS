function [] = set_interpolation_method(target_file, save_as, switch_to)
% Copies the input file filename and saves a copy, but overrides the
%intmethod (interpolation method) flag in the copy.
% CLEARS WORKSPACE VARIABLES ON EXECUTION!
%   target_file : Input file to read from
%   save_as     : File name to save the new input file under
%   switch_to   : Either 1 (cubic) or 2 (band-limited). The interpolation
%   methods to specify in the new input file.

if ((switch_to ~= 1) && (switch_to ~= 2))
    error("Invalid interpolation method to switch to!");
end

load(target_file);
intmethod = switch_to;
save(save_as);
clear; close all;

end
