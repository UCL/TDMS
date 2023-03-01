function [wx, wy] = gauss_pol_base(th, ph, tight, var_val)
    % Roundabout MATLAB way of defining defaults for inputs
    if ~exist('tight', 'var')
        tight = false;
    end
    if ~exist('var_val', 'var')
        % Only one gets used, but it saves a comparison
        OOES = 5e-6
        FWHM = 25e-6
    else
        % Only one gets used, but it saves a comparison
        OOES = var_val;
        FWHM = var_val;
    end

    refind = 1.35;
    lambda = 1300e-9;
    k = 2*pi/lambda;

    if tight
        W = 4/k/refind/OOES;
        wx = exp(-(sin(th)/W).^2);
    else
        W = FWHM/2/sqrt(2*log(2))*k*refind;
        wx = exp( -(W*sin(th)).^2 );
    end

    wy = wx;
end
