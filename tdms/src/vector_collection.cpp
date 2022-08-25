#include <utility>
#include "vector_collection.h"
#include "utils.h"


using namespace std;

void XYZVectors::set_ptr(const char c, double* ptr){
  switch (c) {
    case 'x': {x = ptr; break;}
    case 'y': {y = ptr; break;}
    case 'z': {z = ptr; break;}
    default: throw std::runtime_error("Have no element " + to_string(c));
  }
}

CMaterial::CMaterial(const mxArray *ptr) {

  assert_num_fields_equals(9, ptr, "Cmaterials");

  init_xyz_vectors(ptr, a, "Ca");
  init_xyz_vectors(ptr, b, "Cb");
  init_xyz_vectors(ptr, c, "Cc");
}

void MaterialCollection::init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays, const string &prefix) {

  for (char component : {'x', 'y', 'z'}) {
    auto element = ptr_to_vector_in(ptr, prefix + component, "material");
    arrays.set_ptr(component, mxGetPr(element));
  }
}

void CCollection::init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays, const string &prefix) {

  for (char component : {'x', 'y', 'z'}) {

    auto element = ptr_to_matrix_in(ptr, prefix + component, "C");
    is_multilayer = mxGetDimensions(element)[0] != 1;
    arrays.set_ptr(component, mxGetPr(element));
  }
}

CCollection::CCollection(const mxArray *ptr) {

  auto num_fields = mxGetNumberOfFields(ptr);
  if (num_fields != 6 && num_fields != 9) {
    throw runtime_error("C should have 6 or 9 members, it has " + to_string(num_fields));
  }

  init_xyz_vectors(ptr, a, "Ca");
  init_xyz_vectors(ptr, b, "Cb");

  if (num_fields == 9) {
    is_disp_ml = true;
    init_xyz_vectors(ptr, c, "Cc");
  }
}

DMaterial::DMaterial(const mxArray *ptr) {

  assert_num_fields_equals(6, ptr, "Dmaterials");

  init_xyz_vectors(ptr, a, "Da");
  init_xyz_vectors(ptr, b, "Db");
}

DCollection::DCollection(const mxArray *ptr) {

  assert_num_fields_equals(6, ptr, "D");

  init_xyz_vectors(ptr, a, "Da");
  init_xyz_vectors(ptr, b, "Db");
}

void DCollection::init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays, const string &prefix) {

  for (char component : {'x', 'y', 'z'}) {
    auto element = ptr_to_matrix_in(ptr, prefix + component, "D");
    arrays.set_ptr(component, mxGetPr(element));
  }
}

DispersiveMultiLayer::DispersiveMultiLayer(const mxArray *ptr) {

  if (mxIsEmpty(ptr)) {
    return;
  }
  assert_is_struct_with_n_fields(ptr, 9, "dispersive_aux");

  alpha = mxGetPr(ptr_to_vector_in(ptr, "alpha", "dispersive_aux"));
  beta = mxGetPr(ptr_to_vector_in(ptr, "beta", "dispersive_aux"));
  gamma = mxGetPr(ptr_to_vector_in(ptr, "gamma", "dispersive_aux"));
  kappa.x = mxGetPr(ptr_to_matrix_in(ptr, "kappa_x", "dispersive_aux"));
  kappa.y = mxGetPr(ptr_to_matrix_in(ptr, "kappa_y", "dispersive_aux"));
  kappa.z = mxGetPr(ptr_to_matrix_in(ptr, "kappa_z", "dispersive_aux"));
  sigma.x = mxGetPr(ptr_to_matrix_in(ptr, "sigma_x", "dispersive_aux"));
  sigma.y = mxGetPr(ptr_to_matrix_in(ptr, "sigma_y", "dispersive_aux"));
  sigma.z = mxGetPr(ptr_to_matrix_in(ptr, "sigma_z", "dispersive_aux"));
}
