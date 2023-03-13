%function [x,y,z,lambda] = fdtd_bounds(input_file)
%
%input_file - file with input configuration information
%
%x, y and z are the grid labels of the interior space of the FDTD
%grid.
%
%lambda is the wavelength in dielectric material
function [x,y,z,lambda] = fdtd_bounds(input_file)

%input the configuration information
[fid_input,message] = fopen(input_file,'r');

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

%now need to check that all of the required variables have been set
variables = {'delta','I','J','K','z_launch','illorigin','lambda'};
must_abort = 0; %assumes all variables have been defined
for lvar = 1:length(variables)
    if exist(variables{lvar}) ~= 1
	if strncmp(variables(lvar),'z_launch',8)
	    fprintf(1,'Failed to define %s, setting it to 0\n',variables{lvar});
	    z_launch = 0;
	else
	    fprintf(1,'Failed to define %s\n',variables{lvar});
	    must_abort = 1;
	end
    end
end

if must_abort
    error('Not all variables were defined');
end


x = ((1:I) - illorigin(1))*delta.x;
y = ((1:J) - illorigin(2))*delta.y;
z = ((1:K) - illorigin(3))*delta.z + z_launch;
lambda = lambda;
