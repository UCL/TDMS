#include "arrays/tdms_matrix.h"

#include <stdexcept>
#include <string>

#include <spdlog/spdlog.h>

using std::runtime_error;
using std::to_string;

GratingStructure::GratingStructure(const mxArray *ptr, int I_tot) {

  if (mxIsEmpty(ptr)) { return; }

  auto dims = mxGetDimensions(ptr);
  if (mxGetNumberOfDimensions(ptr) != 2 || dims[0] != 2 ||
      dims[1] != (I_tot + 1)) {
    throw runtime_error("structure should have dimension 2 x " +
                        to_string(I_tot + 1));
  }

  int **matlab_buffer =
          cast_matlab_2D_array((int *) mxGetPr(ptr), 2, I_tot + 1);
  initialise(matlab_buffer, 2, I_tot + 1, true);
  free_cast_matlab_2D_array(matlab_buffer);
}

void Vertices::initialise_from_matlab(const mxArray *ptr) {

  auto element = ptr_to_matrix_in(ptr, "vertices", "campssample");
  if (mxIsEmpty(element)) { return; }

  auto dims = mxGetDimensions(element);

  if (dims[1] != 3) {
    throw runtime_error("Second dimension in campssample.vertices must be 3");
  } else {
    spdlog::info("Found vertices ({0:d} x 3)", dims[0]);
  }

  int **matlab_buffer =
          cast_matlab_2D_array((int *) mxGetPr(element), dims[0], dims[1]);
  initialise(matlab_buffer, dims[0], 3, true);
  free_cast_matlab_2D_array(matlab_buffer);

  // Decrement stored index values, because of MATLAB->C indexing
  std::for_each(data_.begin(), data_.end(), [](int &n) { n--; });
}

ijk Vertices::index_in_row(int row) const {
  return {operator()(row, 0), operator()(row, 1), operator()(row, 2)};
}
