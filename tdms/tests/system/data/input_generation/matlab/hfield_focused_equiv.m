function [H] = hfield_focused_equiv(X,Y,Z)
	%% Set the magnetic field source terms to be 0, over the computational grid X, Y, Z provided.
	%% This function is called by iteratefdtd_matrix to set the magnetic field source terms.
	% X, Y, Z	: Grid coordinates
	%
	% H			: Magnetic field over the grid

	H = cell(1, 3);
	for component = 1:3
		H{component} = zeros(size(X));
	end
end
