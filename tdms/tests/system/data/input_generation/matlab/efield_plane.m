function [E] = efield_plane(X,Y,Z)
    %% Sets up a plane-wave electric field propagating in the z-direction, over the spatial grid X, Y, Z.
    % X, Y, Z           : Vectors defining the coordinate grid, output of ndgrid(x, y, z)
    %
    % E                 : The gaussian electric field produced, defined over the spatial grid.

    %% Preallocate storage
    [m,n] = size(X);
    E = cell(1,2);
    E{1} = zeros(m,n);
    E{2} = zeros(m,n);

    %% Define constants
    lambda = 1300e-9;
    refind = 1.35;
    dz = lambda/10;
    k=2*pi/lambda;
    % dz/2 term due to the modified source condition
    vertices = [X(:) Y(:) Z(:)-dz/2];

    %% Compute field
    E{1} = 2*exp(sqrt(-1)*k*refind*(Z-dz/2));
    E{2} = 2*zeros(size(Z));
end
