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

  x_ptr = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
  x_real = cast_matlab_2D_array(mxGetPr(x_ptr), dims[0], dims[1]);
  x_imag = cast_matlab_2D_array(mxGetPi(x_ptr), dims[0], dims[1]);

  y_ptr = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
  y_real = cast_matlab_2D_array(mxGetPr(y_ptr), dims[0], dims[1]);
  y_imag = cast_matlab_2D_array(mxGetPi(y_ptr), dims[0], dims[1]);

  x = (complex<double> **) malloc(sizeof(complex<double> *) * n_frequencies);
  y = (complex<double> **) malloc(sizeof(complex<double> *) * n_frequencies);

  // zero the arrays we've created
  for (int ifx = 0; ifx < n_frequencies; ifx++) {
    x[ifx] = (complex<double> *) malloc(sizeof(complex<double>) * dims[0]);
    y[ifx] = (complex<double> *) malloc(sizeof(complex<double>) * dims[0]);
    for (int im = 0; im < n_det_modes; im++) {
      x[ifx][im] = 0.;
      y[ifx][im] = 0.;
      x_real[ifx][im] = 0.;
      x_imag[ifx][im] = 0.;
      y_real[ifx][im] = 0.;
      y_imag[ifx][im] = 0.;
    }
  }

  // set the Idx, Idy arrays as the field values for the output struct
  mxSetField(id_pointer, 0, "Idx", x_ptr);
  mxSetField(id_pointer, 0, "Idy", y_ptr);
}

IDVariables::~IDVariables() {
  if (memory_assigned) {
    free_cast_matlab_2D_array(x_real);
    free_cast_matlab_2D_array(x_imag);
    free_cast_matlab_2D_array(y_real);
    free_cast_matlab_2D_array(y_imag);
    for (int ifx = 0; ifx < n_frequencies; ifx++) {
      free(x[ifx]);
      free(y[ifx]);
    }
    free(x);
    free(y);
  }
}
