function [I] = scattering_matrix(X, Y, Z, rad, obj_shape)
    %% Setup the scattering matrix for the given object shape and grid cooordinates.
    % X, Y, Z: Outputs of ndgrid, the coordinates of the Yee cells we will be using
    % rad: Radius of the sph object, or circular-face of a cyl object
    % obj_shape: The shape of the scattering object. Options are;
    %       sph : sphere
    %       cyl : cylinder

    I = zeros(size(X));

    if strcmp(obj_shape, 'sph')
        %set all Yee cells within the sphere to have index of 1
        I( X.^2 + Y.^2 + Z.^2 < rad^2 ) = 1;
        I((end-3):end, 1, :) = 0;
        I(:, 1, (end-3):end) = 0;
    elseif strcmp(obj_shape, 'cyl')
        %set all Yee cells within the cylinder to have index of 1
        I( X.^2 + Z.^2 < rad^2 ) = 1;
        I((end-3):end, 1, :) = 0;
        I(:, 1, (end-3):end) = 0;
    end
end
