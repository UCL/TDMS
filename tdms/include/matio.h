#ifndef MIOFLAG
#include "../matlabio/Mat_io.h"
#include "tdms_iterator.h"
#else
#include "mat.h"
#include "mex.h"
#include "matrix.h"
#endif

#ifdef CLUSTER
typedef int mwSize; 
#endif
