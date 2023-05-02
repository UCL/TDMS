function [epso , muo , c] = import_constants()
    %% Returns the values of physical constants relevant to light propagation.
    %
    % eps0  : Permitivity of free space
    % mu0   : Permeability of free space
    % c     : Speed of light in vacuum

epso = 8.854187817e-12;
muo = 4.0 * pi * 1.0e-7;
c = 1.0 / sqrt(muo*epso);
end
