function [x,y,z,lambda] = fdtd_bounds(input_file)
	%% Compute vectors x, y, z defining the coordinates of the fdtd grid in the respective axial direction.
	%% That is, the Cartesian product x \otimes y \otimes z is the set of all gridpoints (xi, yj, zk) in the fdtd grid.
	% input_file	: The input file with configuration information.
	%
	% x, y, z		: Grid labels (spatial coodinates) of the interior space of the FDTD grid.
	% lambda		: Wavelength of light in the dielectric material.

%% Fetch the configuration information for this test

% Check that the input file can be found on the path
if isfile(input_file)
	% Run input file as a script to import variables into the workspace
	run(input_file);
else
	% Throw error - config file cannot be found
    error(sprintf('File %s could not be opened for reading',input_file));
end

% Check required variables have been set in the config file, and loaded
% Note that lambda will be returned once it is imported from the input file
required_variables = {'delta','I','J','K','illorigin','lambda'};
n_required_variables = length(required_variables);
i = 1;
% Search through all required variables and check they exist in the workspace
% Throw error (and break loop early) if they are not found
while i<=n_required_variables
	if ~exist(required_variables{i}, 'var')
		error(sprintf('Required variable %s not present in %s', required_variables{i}, input_file));
		break
	end
	i = i + 1;
end

% These variables have default values if they are not defined in the input_file
if ~exist('z_launch', 'var')
	% z_launch defines an additional "height" above the illumination origin
	z_launch = 0;
end

%% Now compute the fdtd grid coordinates along each axis. Values are (respective to each axial direction):
% I, J, K 		: Number of gridpoints
% delta,{x,y,z} : Spatial separation of the gridpoints
% illorigin		: Location of the origin of the source field
x = ((1:I) - illorigin(1))*delta.x;
y = ((1:J) - illorigin(2))*delta.y;
z = ((1:K) - illorigin(3))*delta.z + z_launch;

end
