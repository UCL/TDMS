%function [n] = fdtdminsteps(inputfile)
%
%estimates the minmum number of steps required in a fdtd simulation
%assuming the time step used is 0.95xmaximum time step and that the
%trailing edge of the guassian pulse travels at the speed of
%light. (dispersion free). Calculates the edge of the pulse travelling from
%interface to end wall and back.
function [n] = minsteps_fdtd(inputfile)

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

if exist('K') ~= 1
    error('K is not defined - cannot determine n');
end

if exist('delta') ~= 1
    error('delta is not defined - cannot determine n');
end

if exist('interface') ~= 1
    error('interface is not defined - cannot determine n');
end

if exist('f_an') ~= 1
    error('f_an is not defined - cannot determine n');
end


[to hwhm] = fdtdduration(inputfile);
[dt_upper] = fdtdts(inputfile);

%have to have pulse reaching close to 0 at the interface
t = 2*to;

%now this has to propagate from the interface to the other edge of
%the grid and then out again
t = t + (K-interface.K0(1))*delta.z/(c/sqrt(max(epsr))) + K*delta.z/(c/sqrt(max(epsr)));

n = ceil(t/(0.95*dt_upper));
