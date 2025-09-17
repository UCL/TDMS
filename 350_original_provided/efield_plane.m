% function[E] = efield(X, Y, Z);
% % This function is called by iteratefdtd_matrix to set the electric %
        field source terms % X,
        Y and Z should bu in microns function[E] = efield_plane(X, Y, Z)


                [m, n] = size(X);
E = cell(1, 2);
E{1} = zeros(m, n);
E{2} = zeros(m, n);
lambda = 1300e-9;

k = 2 * pi / lambda;

vertices = [X( :) Y( :) Z( :)];

E{1} = exp(sqrt(-1) * k * Z);
E{2} = zeros(size(Z));
