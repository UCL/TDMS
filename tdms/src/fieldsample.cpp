#include "fieldsample.h"

void FieldSample::set_from(const mxArray *ptr) {

  if (mxIsEmpty(ptr)) { return; }

  assert_is_struct_with_n_fields(ptr, 4, "fieldsample");
  i = Vector<int>(ptr_to_vector_or_empty_in(ptr, "i", "fieldsample"));
  j = Vector<int>(ptr_to_vector_or_empty_in(ptr, "j", "fieldsample"));
  k = Vector<int>(ptr_to_vector_or_empty_in(ptr, "k", "fieldsample"));
  n = Vector<double>(ptr_to_vector_or_empty_in(ptr, "n", "fieldsample"));

  int n_dims = 4;
  if (all_vectors_are_non_empty()) {
    int dims[4] = {i.size(), j.size(), k.size(), n.size()};
    mx = mxCreateNumericArray(n_dims, (const mwSize *) dims, mxDOUBLE_CLASS, mxREAL);
    tensor = cast_matlab_4D_array(mxGetPr(mx), i.size(), j.size(), k.size(), n.size());
  } else {
    int dims[4] = {0, 0, 0, 0};
    mx = mxCreateNumericArray(n_dims, (const mwSize *) dims, mxDOUBLE_CLASS, mxREAL);
  }
}

void FieldSample::extract(ElectricSplitField &E_split, PerfectlyMatchedLayer &pml,
                          int n_simulation_timesteps) {
  double Ex_temp = 0., Ey_temp = 0., Ez_temp = 0.;

/* Extract the (electric) field at each of the vertices.
Since the split-field has already been computed, we can do this in parallel by reading the values from the split field and interpolating to the vertices independently of each other.
*/
#pragma omp parallel default(shared) private(Ex_temp, Ey_temp, Ez_temp)
  {
#pragma omp for
    for (int kt = 0; kt < k.size(); kt++) {
      for (int jt = 0; jt < j.size(); jt++) {
        for (int it = 0; it < i.size(); it++) {
          CellCoordinate current_cell(i[it] + pml.Dxl - 1, j[jt] + pml.Dyl - 1,
                                      k[kt] + pml.Dzl - 1);
          Ex_temp = E_split.interpolate_to_centre_of(AxialDirection::X, current_cell);
          if (current_cell.j() != 0) {
            Ey_temp = E_split.interpolate_to_centre_of(AxialDirection::Y, current_cell);
          } else {
            Ey_temp = E_split.yx[current_cell] + E_split.yz[current_cell];
          }
          Ez_temp = E_split.interpolate_to_centre_of(AxialDirection::Z, current_cell);
          for (int nt = 0; nt < n.size(); nt++)
            tensor[nt][kt][jt][it] +=
                    pow(Ex_temp * Ex_temp + Ey_temp * Ey_temp + Ez_temp * Ez_temp, n[nt] / 2.) /
                    n_simulation_timesteps;
        }
      }
    }
  }
}

FieldSample::~FieldSample() {
  if (all_vectors_are_non_empty()) { free_cast_matlab_4D_array(tensor, k.size(), n.size()); }
}
