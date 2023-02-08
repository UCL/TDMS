%[S] = multi_layer( zvec, nvec, lambda0, theta0, sp)
%
%Caculate the characteristic matrix of a multlayer structure
%specified as per the following:
%
% zvec =    [z1    z2   ... zN]
% nvec = [n0    n1    n2  ...  nN]
% lambda0: wavelength in air
% theta0: angle of incidence of plane wave in radians
% sp: 'TE' or 'TM'. 'TE' is taken to mean that the electric vector
%                        in perpendicular to the plane in which the
%                        wave propagates and 'TM' is the case when
%                        the electric field vector is in, or
%                        parallel to, the plane in which the wave
%                        propagates

function [S] = multi_layer( zvec, nvec, lambda0, theta0, sp)


%sp='TE';%TM %s in the TE case (perpendicular) and p the TM case (parallel), see Azzam and
        %Bashara page 271. Essentially, s is perpendicular to the
        %plane that the plane propagates in the p is parallel.

%function [S] = multi_layer( zvec, nvec, lambda0, theta0, sp)
    %perform some basic checks
	if numel(nvec)-numel(zvec) ~= 1
	    error('nvec should have one more element that zvec');
	end

theta_vec=asin(nvec(1)*sin(theta0)./nvec);

st=sin(theta_vec);
ct=cos(theta_vec);

d=diff(zvec);
if numel(nvec)>2
    beta=ct(2:(end-1)).*nvec(2:(end-1)).*d*2*pi/lambda0;
end

if strcmp(sp,'TE')%the 's' case
    t_fresnel = 2*nvec(1:(end-1)).*ct(1:(end-1))./(nvec(1:(end-1)).*ct(1:(end-1))+nvec(2:end).*ct(2:end));
    r_fresnel = (nvec(1:(end-1)).*ct(1:(end-1))-nvec(2:end).*ct(2:end))./(nvec(1:(end-1)).*ct(1:(end-1))+nvec(2:end).*ct(2:end));
else%the 'p' case
    t_fresnel = 2*nvec(1:(end-1)).*ct(1:(end-1))./(nvec(2:end).*ct(1:(end-1))+nvec(1:(end-1)).*ct(2:end));
    r_fresnel = (nvec(2:end).*ct(1:(end-1))-nvec(1:(end-1)).*ct(2:end))./(nvec(2:end).*ct(1:(end-1))+nvec(1:(end-1)).*ct(2:end));
end

%trans_s
%2.*n1*ct1/(n1*ct1 + n2*ct2);
%ref_s
%(n1*ct1 - n2*ct2)/(n1*ct1 + n2*ct2);

%trans_p
%2.*n1*ct1/(n2*ct1 + n1*ct2);
%ref_p
%(n2*ct1 - n1*ct2)/(n2*ct1 + n1*ct2);

S = 1/t_fresnel(1)*[1 r_fresnel(1);r_fresnel(1) 1];
for il=2:numel(zvec)
    %dP=exp(sqrt(-1)*beta(il-1));
    P = [exp(-sqrt(-1)*beta(il-1)) 0;0 exp(sqrt(-1)*beta(il-1))];
    L = 1/t_fresnel(il)*[1 r_fresnel(il);r_fresnel(il) 1];
    S=S*P*L;
end
%S1=S;
%[S] = multi_layer_03( zvec, nvec, lambda0);
%abs(S1-S)
%save ml_vars;
