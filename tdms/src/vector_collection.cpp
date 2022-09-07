#include <utility>
#include "vector_collection.h"


using namespace std;

void XYZVectors::set_ptr(const char c, double* ptr){
  switch (c) {
    case 'x': {x = ptr; break;}
    case 'y': {y = ptr; break;}
    case 'z': {z = ptr; break;}
    default: throw std::runtime_error("Have no element " + std::string(1, c));
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

    auto element = ptr_to_2d_array_in(ptr, prefix + component);
    if (mxGetDimensions(element)[0] != 1) {
      throw runtime_error("Incorrect dimension on material: " + prefix + component);
    }
    arrays.set_ptr(component, mxGetPr(element));
  }
}

void CCollection::init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays, const string &prefix) {

  for (char component : {'x', 'y', 'z'}) {

    auto element = ptr_to_2d_array_in(ptr, prefix + component);
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

    auto element = ptr_to_2d_array_in(ptr, prefix + component);
    arrays.set_ptr(component, mxGetPr(element));
  }
}
