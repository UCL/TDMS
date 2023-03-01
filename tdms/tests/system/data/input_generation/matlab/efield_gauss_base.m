function [E] = efield_gauss_base(X, Y, Z, tight, fstrat_method, gauss_pol_method)
    % Check for defaults MATLAB-style
    if ~exist('tight', 'var')
        tight = false;
    end
    if ~exist('fstrat_method', 'var')
        fstrat_method = @focstratfield_general_pol_2d;
    end
    if ~exist('gauss_pol_method', 'var')
        gauss_pol_method = @(th, ph) gauss_pol_base(th, ph, false);
    end

    [m,n] = size(X);
    E = cell(1,2);
    E{1} = zeros(m,n);
    E{2} = zeros(m,n);
    lambda = 1300e-9;
    k=2*pi/lambda;
    refind = 1.35;
    if tight
        dz = lambda/4;
        ntheta = 100;
        nphi = 100;
    else
        dz = lambda/6;
        ntheta = 200;
        nphi = [];
    end
    %dz/2 term due to the modified source condition
    vertices = [X(:) Y(:) Z(:)-dz/2];
    nvec = refind;
    hvec = [];
    NA = 1;

    %first calculate normalisation
    [EpN,Em] = fstrat_method([0 0 0],nvec,hvec,NA,lambda,ntheta,nphi,gauss_pol_method);

    %calculate the field at the interface
    [Ep,Em] = fstrat_method(vertices,nvec,hvec,NA,lambda,ntheta,nphi,gauss_pol_method);

    %factor of 2 due to the modified source condition
    E{1} = 2*reshape(Ep(:,1)/EpN(1),size(X));
    E{2} = 2*reshape(Ep(:,2)/EpN(1),size(X));
end
