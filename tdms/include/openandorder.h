/**
 * @file openandorder.h
 * @brief Launch TDMS and file IO.
 *
 * Code for processing command line arguments, opening input files,  passing
 * matrices to the mexFunction and writing the output to the specified output
 * file.
 */
#pragma once

#include "argument_parser.h"
#include "input_matrices.h"
#include "mat_io.h"
#include "matrix_collection.h"

/**
 * @brief Save the resultant matrices into a file with name outputfilename
 * Throws a runtime error if the output file name cannot be opened.
 *
 * @param plhs the matrices themselves
 * @param _matricestosave the indexes of the matrices to save
 * @param _matrixnames names of the matrices to save
 * @param nmatrices the number of matrices to save
 * @param outputfilename name of the output file
 */
void saveoutput(mxArray **plhs, const int *_matricestosave, char *_matrixnames[], int nmatrices, const char *outputfilename);

/**
 * @brief Check that all input and output files can be accessed
 *
 * @param args the ArgumentNamespace
 */
void check_files_can_be_accessed(ArgumentNamespace &args);
