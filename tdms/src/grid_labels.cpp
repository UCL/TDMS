#include "grid_labels.h"

#include "matlabio.h"

GridLabels::GridLabels(const mxArray *ptr){

  x = mxGetPr(ptr_to_vector_in(ptr, "x_grid_labels", "grid_labels"));
  y = mxGetPr(ptr_to_vector_in(ptr, "y_grid_labels", "grid_labels"));
  z = mxGetPr(ptr_to_vector_in(ptr, "z_grid_labels", "grid_labels"));
}

void GridLabels::initialise_from(GridLabels &labels_to_copy_from, int i_l, int i_u, int j_l,
                                 int j_u, int k_l, int k_u) {
  for (int i = i_l; i <= i_u; i++) { x[i - i_l] = labels_to_copy_from.x[i]; }
  for (int j = j_l; j <= j_u; j++) { y[j - j_l] = labels_to_copy_from.y[j]; }
  for (int k = k_l; k <= k_u; k++) { z[k - k_l] = labels_to_copy_from.z[k]; }
}
