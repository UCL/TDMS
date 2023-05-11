#include "arrays/dtilde.h"

#include "globals.h"
#include "matlabio.h"
#include "matrix.h"

using namespace std;
using tdms_math_constants::IMAGINARY_UNIT;

void DTilde::set_component(Tensor3D<complex<double>> &tensor,
                           const mxArray *ptr, const string &name, int n_rows,
                           int n_cols) {

  auto element = ptr_to_nd_array_in(ptr, 3, name, "D_tilde");

  auto dims = (int *) mxGetDimensions(element);
  int n_det_modes = dims[0];

  if (dims[1] != n_rows || dims[2] != n_cols) {
    throw runtime_error("D_tilde.{x, y} has final dimensions " +
                        to_string(dims[1]) + "x" + to_string(dims[2]) +
                        " but it needed to be " + to_string(n_rows) + "x" +
                        to_string(n_cols));
  }

  complex<double> ***p =
          (complex<double> ***) malloc(sizeof(complex<double> **) * n_cols);
  for (int j = 0; j < n_cols; j++) {
    p[j] = (complex<double> **) malloc(sizeof(complex<double> *) * n_rows);
    for (int i = 0; i < n_rows; i++) {
      p[j][i] =
              (complex<double> *) malloc(sizeof(complex<double>) * n_det_modes);
    }
  }

  auto temp_re =
          cast_matlab_3D_array(mxGetPr(element), dims[0], dims[1], dims[2]);
  auto temp_im =
          cast_matlab_3D_array(mxGetPi(element), dims[0], dims[1], dims[2]);

  for (int k = 0; k < n_det_modes; k++)
    for (int j = 0; j < n_cols; j++)
      for (int i = 0; i < n_rows; i++) {
        p[j][i][k] = temp_re[j][i][k] + IMAGINARY_UNIT * temp_im[j][i][k];
      }

  free_cast_matlab_3D_array(temp_re, n_cols);
  free_cast_matlab_3D_array(temp_im, n_cols);
  tensor.initialise(p, n_cols, n_rows, n_det_modes, true);
}

void DTilde::initialise(const mxArray *ptr, int n_rows, int n_cols) {

  if (mxIsEmpty(ptr)) { return; }

  assert_is_struct_with_n_fields(ptr, 2, "D_tilde");
  set_component(x, ptr, "Dx_tilde", n_rows, n_cols);
  set_component(y, ptr, "Dy_tilde", n_rows, n_cols);
  n_det_modes =
          mxGetDimensions(ptr_to_nd_array_in(ptr, 3, "Dx_tilde", "D_tilde"))[0];
}
