function [E] = efield_gauss_base(X, Y, Z, tight, gauss_pol_method, ntheta, nphi)
    %% Sets up a Gaussian electric field over the spatial grid X, Y, Z that is passed. Optional parameters can be passed to change the default behaviour.
    % X, Y, Z           : Vectors defining the coordinate grid, output of ndgrid(x, y, z)
    % tight             : Boolean, when set to true, a tighter Gaussian is produced but the number of lens angles is increased. Defaults to false.
    % gauss_pol_method  : Function handle wrapping gauss_pol_base - determines detailed shape of the Gaussian field.
    %                     gauss_pol_based is used as the default.
    % ntheta            : Override for the ntheta variable
    % nphi              : Override for the nphi variable
    % E                 : The gaussian electric field produced, defined over the spatial grid.

    %% Check for optional arguments
    if ~exist('tight', 'var')
        tight = false;
    end
    if ~exist('gauss_pol_method', 'var')
        gauss_pol_method = @(th, ph) gauss_pol_base(th, ph, false);
    end

    %% Prepare output field storage
    [m,n] = size(X);
    E = cell(1,2);
    E{1} = zeros(m,n);
    E{2} = zeros(m,n);

    %% Define constants
    lambda = 1300e-9;
    k=2*pi/lambda;
    refind = 1.35;
    % Adjust depending on tight-ness
    if tight
        dz = lambda/4;
        if ~exist('ntheta', 'var')
            ntheta = 100;
        end
        if ~exist('nphi', 'var')
            nphi = 100;
        end
    else
        dz = lambda/6;
        if ~exist('ntheta', 'var')
            ntheta = 200;
        end
        if ~exist('nphi', 'var')
            nphi = [];
        end
    end
    % dz/2 term due to the modified source condition
    vertices = [X(:) Y(:) Z(:)];
    nvec = refind;
    hvec = [];
    NA = 1;

    %% Compute the field
    % Compute the normalisation of the field
    [EpN,Em] = focstratfield([0 0 0],nvec,hvec,NA,lambda,ntheta,nphi,gauss_pol_method);
    % Then calculate the field at the interface
    [Ep,Em] = focstratfield(vertices,nvec,hvec,NA,lambda,ntheta,nphi,gauss_pol_method);
    % Multiply by due to the modified source condition
    E{1} = reshape(Ep(:,1)/EpN(1),size(X));
    E{2} = reshape(Ep(:,2)/EpN(1),size(X));
end
