function [composition_matrix] = composition_matrix_builder(file_with_coordinates, obj_shape, rad)
    %% Setup the composition matrix for the given object shape and grid cooordinates.
    % x, y, z: Vectors, the tensor/outer product of which forms the set of coordinates of the Yee cells we will be using
    % obj_shape: The shape of the scattering object. Options are;
    %       fs  : free-space, no obstacle
    %       sph : sphere
    %       cyl : cylinder
    %       sc  : point-source at the origin
    % rad: Radius of the sph object, or circular-face of a cyl object. Point sources and freespace ignore this input
    %
    % composition_matrix : The composition matrix for the simulation domain

    if strcmp(obj_shape, 'fs')
        composition_matrix = [];
    else
        % Obtain coordinates of the computational grid
        [x,y,z,lambda] = fdtd_bounds(file_with_coordinates);
        % Account for the possibility that we are defining a cylindrical object, so we need to "loose" a dimension
        if strcmp(obj_shape, 'cyl')
            y = 0;
        end
        % Generate Yee cell coords
        [X,Y,Z] = ndgrid(x,y,z);
        % Create scattering matrix of shape (x \otimes y \otimes z)
        % I_{jk} = 1 if this Yee cell is part of the scattering object, 0 otherwise
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
        elseif strcmp(obj_shape, 'sc')
            %set the Yee cell at the origin to have an index of 1
            I( (X==0) & (Y==0) & (Z==0) ) = 1;
        else
            error(sprintf('Error, unrecognised obstacle shape: %s', obj_shape));
        end
        inds = find(I(:));
        [ii,jj,kk] = ind2sub(size(I), inds);
        composition_matrix = [ii jj kk ones(size(ii))];
    end
end
