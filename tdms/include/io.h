#include "matio.h"


int openandorder(char *matfilename, char **matrixnames, const mxArray **matrixptrs, int nmatrices);

int saveoutput(mxArray **plhs, int *matricestosave, char *matrixnames[], int nmatrics, char *outputfilename);

void freememory( int nmatrices, mxArray **matrixptrs);
