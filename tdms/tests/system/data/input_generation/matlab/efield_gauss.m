%function [E] = efield(X,Y,Z);
%
%This function is called by iteratefdtd_matrix to set the electric
%field source terms
%X, Y and Z should bu in microns
function [E] = efield_gauss(X,Y,Z)


    [m,n] = size(X);
    E = cell(1,2);
    E{1} = zeros(m,n);
    E{2} = zeros(m,n);
    lambda = 1300e-9;
    refind = 1.35;
    dz = lambda/6;

    k=2*pi/lambda;

    %dz/2 term due to the modified source condition
    vertices = [X(:) Y(:) Z(:)-dz/2];
    nvec = refind;
    hvec = [];
    NA = 1;
    ntheta = 200;
    nphi = [];

    %first calculate normalisation
    [EpN,Em] = focstratfield_general_pol_2d([0 0 0],nvec,hvec,NA,lambda,ntheta,nphi,@gauss_pol);

    %calculate the field at the interface
    [Ep,Em] = focstratfield_general_pol_2d(vertices,nvec,hvec,NA,lambda,ntheta,nphi,@gauss_pol);

    %factor of 2 due to the modified source condition
    E{1} = 2*reshape(Ep(:,1)/EpN(1),size(X));
    E{2} = 2*reshape(Ep(:,2)/EpN(1),size(X));
