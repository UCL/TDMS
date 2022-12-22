/**
 * @file input_output_names.h
 * @brief Constant values associated with the format of the input and output files and matrices
 */
#pragma once

#include <vector>
#include <string>

#define NMATRICES 49              //< number of input matrices
#define NOUTMATRICES_WRITE 23     //< number of output matrices to be written to output file
#define NOUTMATRICES_WRITE_ALL 25 //< number of output matrices to be written to output file
#define NOUTMATRICES_PASSED 31    //< number of output matrices passed by mexFunction

/*
These are the names of the arrays that we expect to recieve from the input file (and posibly gridfile) that is passed to tdms on the command line.
In the case where we are given an input file and a gridfile, the array names we expect to recieve from each is also recorded.
In the case where the -m (compressed output) flag is passed, we do not save all the possible output arrays and so have a complete list of outputs and a shortened list of those which are saved when using -m.
*/
namespace tdms_matrix_names {
    // all matrices that we could expect to get from an input file with a separate gridfile
    extern const std::vector<std::string> matrixnames_infile;
    // matrices we expect to get from an input file that also contains grid information
    extern const std::vector<std::string> matrixnames_input_with_grid;
    // matrices we expect to obtain from a separate gridfile
    extern const std::vector<std::string> matrixnames_gridfile;
    // oall output matrices we might want to write
    extern const std::vector<std::string> outputmatrices_all;
    // output matrices we want to write in a compressed (-m) output
    extern const std::vector<std::string> outputmatrices;
    // indices of the matrices to save when saving a complete system
    extern const int matricestosave_all[NOUTMATRICES_WRITE_ALL];
    // indices of the matrices to save when saving a compressed output
    extern const int matricestosave[NOUTMATRICES_WRITE];
}
