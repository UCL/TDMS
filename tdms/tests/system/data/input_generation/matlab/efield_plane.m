function [E] = efield_plane(X,Y,Z)

    [m,n] = size(X);
    E = cell(1,2);
    E{1} = zeros(m,n);
    E{2} = zeros(m,n);
    lambda = 1300e-9;
    refind = 1.35;
    dz = lambda/10;

    k=2*pi/lambda;

    %dz/2 term due to the modified source condition
    vertices = [X(:) Y(:) Z(:)-dz/2];

    E{1} = 2*exp(sqrt(-1)*k*refind*(Z-dz/2));
    E{2} = 2*zeros(size(Z));
end
