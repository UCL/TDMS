#include "arrays/tdms_matrix.h"

#include <stdexcept>
#include <string>

using std::runtime_error;
using std::to_string;

GratingStructure::GratingStructure(const mxArray *ptr, int I_tot) {

  if (mxIsEmpty(ptr)) { return; }

  auto dims = mxGetDimensions(ptr);
  if (mxGetNumberOfDimensions(ptr) != 2 || dims[0] != 2 ||
      dims[1] != (I_tot + 1)) {
    throw runtime_error("structure should have dimension 2 x (I_tot+1) ");
  }

  int **matlab_buffer =
          cast_matlab_2D_array((int *) mxGetPr(ptr), 2, I_tot + 1);
  initialise(matlab_buffer, 2, I_tot + 1, true);
  free_cast_matlab_2D_array(matlab_buffer);
}

void Pupil::initialise_from_matlab(const mxArray *ptr, int n_rows, int n_cols) {

  if (mxIsEmpty(ptr)) { return; }

  auto dims = (int *) mxGetDimensions(ptr);

  if (mxGetNumberOfDimensions(ptr) != 2 || dims[0] != n_rows ||
      dims[1] != n_cols) {
    throw runtime_error("Pupil has dimension " + to_string(dims[0]) + "x" +
                        to_string(dims[1]) + " but it needed to be " +
                        to_string(n_rows) + "x" + to_string(n_cols));
  }

  double **matlab_buffer = cast_matlab_2D_array(mxGetPr(ptr), n_rows, n_cols);
  initialise(matlab_buffer, n_rows, n_cols, true);
  free_cast_matlab_2D_array(matlab_buffer);
}
