function [S] = multi_layer(zvec, nvec, lambda0, theta0, sp)
    %% Compute the characteristic matrix of a multilayer structure.
    % zvec      : [    z1, z2, ..., zN].
    % nvec      : [n0, n1, n2, ..., nN].
    % lambda0   : Wavelength of light in air.
    % theta0    : Angle of incidence of the plane wave, in radians.
    % sp        : Either 'TE' or 'TM'. 'TE' is taken to mean that the electric vector is perpendicular to the plane in which the wave propagates. 'TM' is when the electric field vector is in, or parallel to, the plane in which the wave propagates.
    %
    % S         : Characteristic matrix of the multilayer structure.

%% Check dimensionality of inputs
if numel(nvec)-numel(zvec) ~= 1
    error('nvec should have one more element that zvec');
end

%% Compute constant values
theta_vec=asin(nvec(1)*sin(theta0)./nvec);
st=sin(theta_vec);
ct=cos(theta_vec);

d=diff(zvec);
if numel(nvec)>2
    beta=ct(2:(end-1)).*nvec(2:(end-1)).*d*2*pi/lambda0;
end

if strcmp(sp,'TE')
    t_fresnel = 2*nvec(1:(end-1)).*ct(1:(end-1))./(nvec(1:(end-1)).*ct(1:(end-1))+nvec(2:end).*ct(2:end));
    r_fresnel = (nvec(1:(end-1)).*ct(1:(end-1))-nvec(2:end).*ct(2:end))./(nvec(1:(end-1)).*ct(1:(end-1))+nvec(2:end).*ct(2:end));
else
    t_fresnel = 2*nvec(1:(end-1)).*ct(1:(end-1))./(nvec(2:end).*ct(1:(end-1))+nvec(1:(end-1)).*ct(2:end));
    r_fresnel = (nvec(2:end).*ct(1:(end-1))-nvec(1:(end-1)).*ct(2:end))./(nvec(2:end).*ct(1:(end-1))+nvec(1:(end-1)).*ct(2:end));
end

%% Setup scattering matrix
S = 1/t_fresnel(1)*[1 r_fresnel(1);r_fresnel(1) 1];
for il=2:numel(zvec)
    P = [exp(-sqrt(-1)*beta(il-1)) 0;0 exp(sqrt(-1)*beta(il-1))];
    L = 1/t_fresnel(il)*[1 r_fresnel(il);r_fresnel(il) 1];
    S=S*P*L;
end
end
