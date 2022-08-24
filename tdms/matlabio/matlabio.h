#include <complex>
#include <string>
#include "mat_io.h"

double ****castMatlab4DArray(double *array, int nrows, int ncols, int nlayers, int nblocks);
void freeCastMatlab4DArray(double ****castArray, int nlayers, int nblocks);
double ***castMatlab3DArray(double *array, int nrows, int ncols, int nlayers);
void freeCastMatlab3DArray(double ***castArray, int nlayers);
double **castMatlab2DArray(double *array, int nrows, int ncols);
void freeCastMatlab2DArray(double **castArray);
unsigned char ***castMatlab3DArrayUint8(unsigned char *array, int nrows, int ncols, int nlayers);
void freeCastMatlab3DArrayUint8(unsigned char ***castArray, int nlayers);
int ***castMatlab3DArrayInt(int *array, int nrows, int ncols, int nlayers);
void freeCastMatlab3DArrayInt(int ***castArray, int nlayers);
int **castMatlab2DArrayInt(int *array, int nrows, int ncols);
void freeCastMatlab2DArrayInt(int **castArray);

void assert_is_struct(const mxArray* ptr, const std::string &name);

void assert_num_fields_equals(int num, const mxArray* ptr, const std::string &name);

/**
 * Get a pointer to a matrix within a struct with a given name. Throws a runtime error if the
 * resulting tensor is not two dimensional.
 * @param ptr Pointer to the struct
 * @param name Name of the attribute
 * @param struct_name Name of the struct
 * @return Pointer to the matrix
 */
mxArray* ptr_to_matrix_in(const mxArray* ptr, const std::string &name, const std::string &struct_name);

/**
 * Get a pointer to an array within a struct with a given name. Throws a runtime error if the
 * resulting tensor is not a vector with dimensions 1xN.
 * @param ptr Pointer to the struct
 * @param name Name of the attribute
 * @param struct_name Name of the struct
 * @return Pointer to the vector
 */
mxArray* ptr_to_vector_in(const mxArray* ptr, const std::string &name, const std::string &struct_name);

mxArray* ptr_to_vector_or_empty_in(const mxArray* ptr, const std::string &name, const std::string &struct_name);
