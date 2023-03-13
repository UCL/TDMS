%function [x y z] = yeeposition(i,j,k,delta,component)
%
%Calculate the position in a cartesian grid of the field
%component specified by the yee cell index (i,j,k). delta
%is a struct with members:
%delta.x, delta.y and delta.z which are the yee cell
%dimension and component is the component of interest
function [x, y, z] = yeeposition(i,j,k,delta,component)

    if strcmp(component,'Ex')
	x = (i+0.5)*delta.x;
	y = j*delta.y;
	z = k*delta.z;
    elseif strcmp(component,'Ey')
	x = i*delta.x;
	y = (j+0.5)*delta.y;
	z = k*delta.z;
    elseif strcmp(component,'Ez')
	x = i*delta.x;
	y = j*delta.y;
	z = (k+0.5)*delta.z;
    elseif strcmp(component,'Hx')
	x = i*delta.x;
	y = (j+0.5)*delta.y;
	z = (k+0.5)*delta.z;
    elseif strcmp(component,'Hy')
	x = (i+0.5)*delta.x;
	y = j*delta.y;
	z = (k+0.5)*delta.z;
    elseif strcmp(component,'Hz')
	x = (i+0.5)*delta.x;
	y = (j+0.5)*delta.y;
	z = k*delta.z;
    end
