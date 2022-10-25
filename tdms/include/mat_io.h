/**
 * @file mat_io.h
 * @brief Includes MATLAB headers for I/O.
 */
#pragma once

#include <mat.h>
#include <mex.h>
#include <matrix.h>
// https://stackoverflow.com/questions/27069676/how-to-dynamically-create-an-array-of-mxarray-in-a-matlab-mex-file

#ifdef CLUSTER
typedef int mwSize;
#endif
