function D = gaussian_detfun(X,Y,beta,mode)
%Note that this needs to be defined in the space of the sample, ie,
%    magnification = 4;%f2/f1
      MFD = 9.2e-6;
      %Wvec = [0.5e-6 MFD/2 MFD];
      %this is in the sample space
      xp = (0:30)*1e-6;
      D = exp( -( ((X-xp(mode))/(MFD/beta)).^2 + (Y/(MFD/beta)).^2) );
end
