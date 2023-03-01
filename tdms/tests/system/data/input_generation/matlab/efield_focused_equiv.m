%function [E] = efield(X,Y,Z);
%
%This function is called by iteratefdtd_matrix to set the electric
%field source terms
%X, Y and Z should bu in microns
function [E] = efield_focused_equiv(X,Y,Z)
    E = {zeros(size(X)),zeros(size(X)),zeros(size(X))};
