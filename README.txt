################################################################
#Introduction
################################################################

The majority of the code's functionality is contained within file:
iterater/iterater_OMP_PSTD.cpp

All source files used are listed in the make file:
iterater/MakefileOMP_PSTD_UCL2

Demonstration code is located in:
tests/arc_01/run_pstd_bscan.m
If this script is run in Matlab, and the executeable has been
correctly compiled, the demonstration simulation will run,
demonstrating all key steps in running a simulation.

################################################################
#Prerequisites for compiling the code
################################################################
Prior to compiling the code it is necessary to install the following:

-FFTW
-Matlab libraries as distributed with an installation of Matlab that
 are typically located in a directory such as /usr/local/MATLAB/R2019b/bin/glnxa64

################################################################
#Compilation instructions
################################################################
The following instructions correspond to a system running Ubuntu
linux. The Makefile and preliminary commands will need to be changed
according to the system's Matlab installation.

The code is compiled by moving into the "iterater" directory.

Then the following commands should be executed:

LD_LIBRARY_PATH=/usr/local/MATLAB/R2019b/bin/glnxa64:/usr/local/MATLAB/R2019b/sys/os/glnxa64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
make -f MakefileOMP_PSTD_UCL2

This will result in an executeable "prtmfdtdOMP_PSTD"

################################################################
#Running the demonstration code
################################################################
Once the executeable has been compiled, move into directory:
tests/arc_01

launch Matlab and run the Matlab script:
run_pstd_bscan.m

This script will generate the input to the executeable, run the
executeable and display sample output.
