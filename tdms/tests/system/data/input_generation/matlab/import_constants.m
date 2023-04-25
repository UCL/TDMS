%function [epso muo c] = import_constants
%
%This is where all physical constants are defined
function [epso , muo , c] = import_constants

epso = 8.854187817e-12;
muo = 4.0 * pi * 1.0e-7;
c = 1.0 / sqrt(muo*epso);
