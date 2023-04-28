function [saveas] = calc_field_tdfield(input_filename, saveas)

    if ~exist('saveas', 'var')
        saveas = 'eivars.mat'
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%All of these definitions need to be copied into the polfun (see below)
%general definitions
[epso , muo , c] = import_constants;
lambda0 = 1300e-9;
omega0 = 2*pi*c/lambda0;
%%%%%%%%%%%%%%%%%%%%%%%%%%

[ex_coords, ey_coords, tvec_E, fvec_E, f_an, hwhm, to_l] = getsourcecoords(input_filename);

lambdavec_E = c./fvec_E;

%the (ex_coords.z-1300e-9/6/2)/c term is necessary to compensate
%for the shift in time introduced by the phase ramp on the Exi and
%Eyi terms below
dz = 1300e-9/6/2;
H = -1*sqrt(-1)*exp(-pi*hwhm^2*(fvec_E - f_an).^2).*exp(sqrt(-1)*2*pi*fvec_E*(to_l-(ex_coords.z-dz)/(c/1.35)))*hwhm;
HN = H/max(abs(H));

Exi = zeros(numel(ex_coords.x),numel(ex_coords.y),numel(ex_coords.z),numel(fvec_E));
Eyi = zeros(numel(ex_coords.x),numel(ex_coords.y),numel(ex_coords.z),numel(fvec_E));

%Setup a plane wave, linearly polarised in the x-directin,
for il=1:numel(lambdavec_E)
    Exi(:,1,1,il) = 2*exp(sqrt(-1)*2*pi/lambdavec_E(il)*1.35*(ex_coords.z-dz));
end

exi = zeros(size(Exi));
eyi = zeros(size(Eyi));

for ix = 1:numel(ex_coords.x)
    et = real(fft(squeeze(Exi(ix,1,1,:)).*H.')*diff(fvec_E(1:2)));
    exi(ix,1,1,:) = [et(2:end).' et(1)];

    et = real(fft(squeeze(Eyi(ix,1,1,:)).*H.')*diff(fvec_E(1:2)));
    eyi(ix,1,1,:) = [et(2:end).' et(1)];
end

exi_store = exi;
exi = zeros( size(exi_store,1), size(exi_store,2), size(exi_store,4) );
exi(:,:,:) = exi_store(:,:,1,:);

eyi_store = eyi;
eyi = zeros( size(eyi_store,1), size(eyi_store,2), size(eyi_store,4) );
eyi(:,:,:) = eyi_store(:,:,1,:);


save saveas exi eyi;

end
