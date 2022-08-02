#include "mat_io.h"
#include "argument_parser.h"
#include "matrix_collection.h"


void openandorder(const char *mat_filename, char **matrix_names, const mxArray **matrix_ptrs, int n_matrices);

void saveoutput(mxArray **plhs, const int *matricestosave, char *matrixnames[], int nmatrices, const char *outputfilename);

void check_files_can_be_accessed(ArgumentNamespace &args);

void assign_matrix_pointers(MatrixCollection &expected, MatFileMatrixCollection &actual, const mxArray **pointers);
