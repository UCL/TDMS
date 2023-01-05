#pragma once

#include "matrix.h"

/* MATLAB DANGER ZONE: Pointers that are created via malloc-like commands, but are linked to would-be outputs. These are destroyed when their parent arrays are destroyed via mxFree, however this can't be done by this class since we need to pass back the output data before we delete it. Need to check that mxFree is called on plhs after writing outputs! */
class Iterator_IntermediateMATLABVariables {
    public:
      mxArray *mx_Idx, *mx_Idy;//< Hold the arrays in the Idx and Idy fields of plhs[26] (Id output)
      double **Idx_re, **Idy_re, **Idx_im, **Idy_im;//< Point to the real (re) and imaginary (im) parts of the data in the Id output
};
