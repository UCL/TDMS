function [wx,wy] = gauss_pol_tight(th,ph);
%    global FWHM;
    refind = 1.35;
    lambda = 1300e-9;
    %one_on_e_squared = 50e-6;
    %FWHM = one_on_e_squared*sqrt(2*log(2))/2
    %FWHM = 0e-6;
%    FWHM = 40e-6;
    OOES = 5e-6;

    k = 2*pi/lambda;

    %W = FWHM/2/sqrt(2*log(2))*k*refind;
    W = 4/k/refind/OOES;

    wx = exp( -(sin(th)/W).^2 );
    wy = wx;
