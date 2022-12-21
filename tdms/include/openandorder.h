/**
 * @file openandorder.h
 * @brief Launch TDMS and file IO.
 *
 * Code for processing command line arguments, opening input files,  passing
 * matrices to the mexFunction and writing the output to the specified output
 * file.
 */
#pragma once

#include "mat_io.h"
#include "argument_parser.h"
#include "matrix_collection.h"

/**
 * @brief Open the input mat file and check they are as expected
 *
 * @param mat_filename The MATLAB filename
 * @param matrix_names names of the matrices in the MATLAB file
 * @param matrix_ptrs pointers to arrays in C++ that will be populated by the MATLAB matrices
 * @param n_matrices the number of matrices in the MATLAB file
 */
void openandorder(const char *mat_filename,
                  const std::vector<std::string>& matrix_names,
                  const mxArray **matrix_ptrs, int n_matrices);

/**
 * @brief Save the resultant matrices into a file with name outputfilename
 * Throws a runtime error if the output file name cannot be opened.
 *
 * @param plhs the matrices themselves
 * @param matricestosave the indexes of the matrices to save
 * @param matrixnames names of the matrices to save
 * @param nmatrices the number of matrices to save
 * @param outputfilename name of the output file
 */
void saveoutput(mxArray **plhs, const int *matricestosave,
                const std::vector<std::string>& matrixnames, int nmatrices, const char *outputfilename);

/**
 * @brief Check that all input and output files can be accessed
 *
 * @param args the ArgumentNamespace
 */
void check_files_can_be_accessed(ArgumentNamespace &args);

/**
 * @brief Iterate through the matrix names and assign the pointer to each matrix
 * into the appropriate entry of pointers
 *
 * @param expected The names of the matrices in the MATLAB file.
 * @param actual The actual matrices in the MATLAB file.
 * @param pointers Pointers to the matrices.
 */
void assign_matrix_pointers(MatrixCollection &expected, MatFileMatrixCollection &actual, const mxArray **pointers);
