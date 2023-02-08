function [wx,wy] = gauss_pol(th,ph);
%    global FWHM;
    refind = 1.35;
    lambda = 1300e-9;
    FWHM = 5e-6;

    k = 2*pi/lambda;
    W = FWHM/2/sqrt(2*log(2))*k*refind;
    
    wx = exp( -(W*sin(th)).^2 );
    wy = wx;
