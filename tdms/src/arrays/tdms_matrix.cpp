#include "arrays/tdms_matrix.h"

#include <stdexcept>

using std::runtime_error;

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
