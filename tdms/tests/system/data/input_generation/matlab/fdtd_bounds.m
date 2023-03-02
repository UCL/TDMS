function [x,y,z,lambda] = fdtd_bounds(input_file)
	%% Compute vectors x, y, z defining the coordinates of the fdtd grid in the respective axial direction.
	%% That is, the Cartesian product x \otimes y \otimes z is the set of all gridpoints (xi, yj, zk) in the fdtd grid.
	% input_file	: The input file with configuration information.
	%
	% x, y, z		: Grid labels (spatial coodinates) of the interior space of the FDTD grid.
	% lambda		: Wavelength of light in the dielectric material.

%% Fetch the configuration information for this test
[delta,I,J,K,illorigin,lambda,z_launch] = get_from_input_file(input_file, struct('z_launch', 0), ...
														'delta','I','J','K','illorigin','lambda','z_launch');

%% Now compute the fdtd grid coordinates along each axis. Values are (respective to each axial direction):
% I, J, K 		: Number of gridpoints
% delta,{x,y,z} : Spatial separation of the gridpoints
% illorigin		: Location of the origin of the source field
x = ((1:I) - illorigin(1))*delta.x;
y = ((1:J) - illorigin(2))*delta.y;
z = ((1:K) - illorigin(3))*delta.z + z_launch;
end
