function [x, y, z] = yeeposition(i,j,k,delta,component)
	%% Calculate the position in cartesian space of the field component associated to the Yee cell with index (i,j,k).
	%% E-field components are offset by 0.5 * the Yee cell extent in their dimension.
	%% H-field components are offset by 0.5 * the Yee cell extent in the OTHER two dimensions.
	% i, j, k	: Indicies of the Yee cell.
	% delta		: Struct with members x, y, z, which are the Yee cell extent in the corresponding axial direction.
	% component	: String, the field component whose position we wish to know.
	%
	% x, y, z	: The cartesian coordinates of the field component.

%% Field component coordinate alligns with Yee cell centre for all other components
x = i * delta.x;
y = j * delta.y;
z = k * delta.z;
%% Add offsets based on which component we are examining
if strcmp(component,'Ex')
	x = x + 0.5*delta.x;
elseif strcmp(component,'Ey')
	y = y + 0.5*delta.y;
elseif strcmp(component,'Ez')
	z = z + 0.5*delta.z;
elseif strcmp(component,'Hx')
	y = y + 0.5*delta.y;
	z = z + 0.5*delta.z;
elseif strcmp(component,'Hy')
	x = x + 0.5*delta.x;
	z = z + 0.5*delta.z;
elseif strcmp(component,'Hz')
	x = x + 0.5*delta.x;
	y = y + 0.5*delta.y;
end
end
