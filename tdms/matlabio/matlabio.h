#include "complex"

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

