%function [E] = efield(X,Y,Z);
%
%This function is called by iteratefdtd_matrix to set the electric
%field source terms
%X, Y and Z should be in microns
function [E] = efield_gauss_tight(X,Y,Z)

    [m,n] = size(X);
    E = cell(1,2);
    E{1} = zeros(m,n);
    E{2} = zeros(m,n);
    lambda = 1300e-9;
    refind = 1.35;
    dz = lambda/4;

    k=2*pi/lambda;

    %dz/2 term due to the modified source condition
    vertices = [X(:) Y(:) Z(:)-dz/2];
    nvec = refind;
    hvec = [];
    NA = 1;
    ntheta = 100;
    nphi = 100;

    %first calculate normalisation
    [EpN,Em] = focstratfield_general_pol([0 0 0],nvec,hvec,NA,lambda,ntheta,nphi,@gauss_pol_tight);

    %calculate the field at the interface
    [Ep,Em] = focstratfield_general_pol(vertices,nvec,hvec,NA,lambda,ntheta,nphi,@gauss_pol_tight);

    %factor of 2 due to the modified source condition
    E{1} = 2*reshape(Ep(:,1)/EpN(1),size(X));
    E{2} = 2*reshape(Ep(:,2)/EpN(1),size(X));
end

function [wx,wy] = gauss_pol_tight(th,ph);
    refind = 1.35;
    lambda = 1300e-9;
    OOES = 5e-6;

    k = 2*pi/lambda;
    W = 4/k/refind/OOES;

    wx = exp( -(sin(th)/W).^2 );
    wy = wx;
end
