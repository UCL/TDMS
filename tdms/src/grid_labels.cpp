#include "grid_labels.h"


GridLabels::GridLabels(const mxArray *ptr){

  x = mxGetPr(ptr_to_vector_in(ptr, "x_grid_labels", "grid_labels"));
  y = mxGetPr(ptr_to_vector_in(ptr, "y_grid_labels", "grid_labels"));
  z = mxGetPr(ptr_to_vector_in(ptr, "z_grid_labels", "grid_labels"));
}
