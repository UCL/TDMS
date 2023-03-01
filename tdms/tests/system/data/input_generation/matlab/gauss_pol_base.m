function [wx, wy] = gauss_pol_base(th, ph, tight)
    % Roundabout MATLAB way of defining defaults for inputs
    if ~exist('tight', 'var')
        tight = false;
    end

    refind = 1.35;
    lambda = 1300e-9;
    k = 2*pi/lambda;

    if tight
        OOES = 5e-6;
        W = 4/k/refind/OOES;
        wx = exp(-(sin(th)/W).^2);
    else
        FWHM = 25e-6;
        W = FWHM/2/sqrt(2*log(2))*k*refind;
        wx = exp( -(W*sin(th)).^2 );
    end

    wy = wx;
end
