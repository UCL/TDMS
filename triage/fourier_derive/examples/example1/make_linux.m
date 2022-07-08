%mex -v -f /opt/matlab/6.0/bin/cxxopts.sh  -I/opt/matlab/6.0/extern/include  -I/home/ptpc2/prmunro/code/ptws1/imaging/illumination_interface_general/generic_focussing GaussLaguerre.cpp /home/ptpc2/prmunro/code/ptws1/imaging/illumination_interface_general/generic_focussing/focussing_terms.cpp /home/ptpc2/prmunro/code/ptws1/matlablibrary/matlabio/matlabio.a
%mex -v -I./ -I/home/pmunro/UCL/imperial/code/ptws1/matlablibrary/matlabio3/  -I/  numerical_derivative.cpp -lfftw3 -lm
mex -v -I./ -I/home/pmunro/UCL/imperial/code/ptws1/matlablibrary/matlabio3/  -I/  test_mdfftw.cpp -lfftw3 -lm
