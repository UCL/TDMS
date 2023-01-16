#include "id_variables.h"

#include "matlabio.h"

using namespace std;

void IDVariables::link_to_pointer(mxArray *&id_pointer, int _n_frequencies, int _n_det_modes) {
  // save array sizes for deallocation
  n_frequencies = _n_frequencies;
  n_det_modes = _n_det_modes;
  // flag memory assignment
  memory_assigned = true;

  // now construct the fields for the structure array, and assign pointers to the data
  int ndims = 2;
  int dims[2] = {n_det_modes, n_frequencies};

  mx_Idx = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
  Idx_re = cast_matlab_2D_array(mxGetPr(mx_Idx), dims[0], dims[1]);
  Idx_im = cast_matlab_2D_array(mxGetPi(mx_Idx), dims[0], dims[1]);

  mx_Idy = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
  Idy_re = cast_matlab_2D_array(mxGetPr(mx_Idy), dims[0], dims[1]);
  Idy_im = cast_matlab_2D_array(mxGetPi(mx_Idy), dims[0], dims[1]);

  Idx = (complex<double> **) malloc(sizeof(complex<double> *) * n_frequencies);
  Idy = (complex<double> **) malloc(sizeof(complex<double> *) * n_frequencies);

  // zero the arrays we've created
  for (int ifx = 0; ifx < n_frequencies; ifx++) {
    Idx[ifx] = (complex<double> *) malloc(sizeof(complex<double>) * dims[0]);
    Idy[ifx] = (complex<double> *) malloc(sizeof(complex<double>) * dims[0]);
    for (int im = 0; im < n_det_modes; im++) {
      Idx[ifx][im] = 0.;
      Idy[ifx][im] = 0.;
      Idx_re[ifx][im] = 0.;
      Idx_im[ifx][im] = 0.;
      Idy_re[ifx][im] = 0.;
      Idy_im[ifx][im] = 0.;
    }
  }

  // set the Idx, Idy arrays as the field values for the output struct
  mxSetField(id_pointer, 0, "Idx", mx_Idx);
  mxSetField(id_pointer, 0, "Idy", mx_Idy);
}

IDVariables::~IDVariables() {
  if (memory_assigned) {
    free_cast_matlab_2D_array(Idx_re);
    free_cast_matlab_2D_array(Idx_im);
    free_cast_matlab_2D_array(Idy_re);
    free_cast_matlab_2D_array(Idy_im);
    for (int ifx = 0; ifx < n_frequencies; ifx++) {
      free(Idx[ifx]);
      free(Idy[ifx]);
    }
    free(Idx);
    free(Idy);
  }
}
