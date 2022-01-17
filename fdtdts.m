%function [dt_upper] = fdtdts(inputfile)
%
%Calculates the maximum allowable time step subject to the
%stability criterion
function [dt_upper] = fdtdts(inputfile)

[fid_input,message] = fopen(inputfile,'r');

%check if file was opened successfully
if fid_input== -1
    error(sprintf('File %s could not be opened for reading',input_file));
end

%proceed to_l read in config information
current_line = fgets(fid_input);

while current_line ~= -1
    eval(current_line);
    current_line = fgets(fid_input);
end

[epso muo c] = import_constants;

if exist('delta') ~= 1
    error('delta is not defined - cannot determine dt_upper');
end

dt_upper = 1/(c*sqrt(1/delta.x^2 + 1/delta.y^2 + 1/delta.z^2));
