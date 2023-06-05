#include "arrays.h"

#include <iostream>
#include <spdlog/spdlog.h>
#include <utility>

#include "globals.h"
#include "utils.h"

using namespace std;
using namespace tdms_math_constants;

void XYZVectors::set_ptr(const char c, double *ptr) {
  switch (c) {
    case 'x': {
      x = ptr;
      break;
    }
    case 'y': {
      y = ptr;
      break;
    }
    case 'z': {
      z = ptr;
      break;
    }
    default:
      throw std::runtime_error("Have no element " + std::string(1, c));
  }
}
void XYZVectors::set_ptr(AxialDirection d, double *ptr) {
  switch (d) {
    case AxialDirection::X: {
      x = ptr;
      break;
    }
    case AxialDirection::Y: {
      y = ptr;
      break;
    }
    case AxialDirection::Z: {
      z = ptr;
      break;
    }
    default:
      throw std::runtime_error("Have no element " + to_string(d));
  }
}

bool XYZVectors::all_elements_less_than(double comparison_value,
                                        int vector_length,
                                        AxialDirection component,
                                        int buffer_start) const {
  double *component_pointer;
  switch (component) {
    case AxialDirection::X:
      component_pointer = x;
      break;
    case AxialDirection::Y:
      component_pointer = y;
      break;
    case AxialDirection::Z:
      component_pointer = z;
      break;
    default:
      throw runtime_error("Error - component not recognised");
      break;
  }
  for (int index = buffer_start; index < vector_length; index++) {
    if (component_pointer[buffer_start + index] > comparison_value) {
      return false;
    }
  }
  return true;
}
bool XYZVectors::all_elements_less_than(double comparison_value, int nx, int ny,
                                        int nz) const {
  if (!all_elements_less_than(comparison_value, nx, AxialDirection::X)) {
    return false;
  }
  if (!all_elements_less_than(comparison_value, ny, AxialDirection::Y)) {
    return false;
  }
  if (!all_elements_less_than(comparison_value, nz, AxialDirection::Z)) {
    return false;
  }
  return true;
}

CMaterial::CMaterial(const mxArray *ptr) {

  assert_num_fields_equals(9, ptr, "Cmaterials");

  init_xyz_vectors(ptr, a, "Ca");
  init_xyz_vectors(ptr, b, "Cb");
  init_xyz_vectors(ptr, c, "Cc");
}

void MaterialCollection::init_xyz_vectors(const mxArray *ptr,
                                          XYZVectors &arrays,
                                          const string &prefix) {

  for (char component : {'x', 'y', 'z'}) {
    auto element = ptr_to_vector_in(ptr, prefix + component, "material");
    arrays.set_ptr(component, mxGetPr(element));
  }
}

void CCollection::init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays,
                                   const string &prefix) {

  for (char component : {'x', 'y', 'z'}) {

    auto element = ptr_to_matrix_in(ptr, prefix + component, "C");
    is_multilayer = mxGetDimensions(element)[0] !=
                    1;// this only matters when we check the 'z' component
                      // right? No point re-setting it each time when it's not
                      // used? TODO: check this.
    arrays.set_ptr(component, mxGetPr(element));
  }
}

CCollection::CCollection(const mxArray *ptr) {

  auto num_fields = mxGetNumberOfFields(ptr);
  if (num_fields != 6 && num_fields != 9) {
    throw runtime_error("C should have 6 or 9 members, it has " +
                        to_string(num_fields));
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

void DCollection::init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays,
                                   const string &prefix) {

  for (char component : {'x', 'y', 'z'}) {
    auto element = ptr_to_matrix_in(ptr, prefix + component, "D");
    arrays.set_ptr(component, mxGetPr(element));
  }
}

bool DispersiveMultiLayer::is_dispersive(double near_zero_tolerance) const {
  for (double gamma_val : gamma) {
    if (fabs(gamma_val) > near_zero_tolerance) {
      // non-zero attenuation constant of a Yee cell implies media is dispersive
      return true;
    }
  }
  return false;
}

FrequencyExtractVector::FrequencyExtractVector(const mxArray *ptr,
                                               double omega_an) {

  if (mxIsEmpty(ptr)) {
    n = 1;
    vector = (double *) malloc(sizeof(double));
    vector[0] = omega_an / 2. / DCPI;

  } else {
    auto dims = mxGetDimensions(ptr);
    auto n_dims = mxGetNumberOfDimensions(ptr);

    if (n_dims != 2 || !(dims[0] == 1 || dims[1] == 1)) {
      throw runtime_error("f_ex_vec should be a vector with N>0 elements");
    }
    // compute the number of elements prior to displaying
    n = std::max(dims[0], dims[1]);
    spdlog::info("f_ex_vec has ndims={} N={}", n_dims, n);
    vector = (double *) mxGetPr(ptr);
  }
}

double FrequencyExtractVector::max() {
  double tmp = -DBL_MAX;
  for (int i = 0; i < n; i++) { tmp = std::max(tmp, vector[i]); }
  return tmp;
}

void FieldComponentsVector::initialise(const mxArray *ptr) {

  auto element = ptr_to_matrix_in(ptr, "components", "campssample");
  if (mxIsEmpty(element)) { return; }

  auto dims = mxGetDimensions(element);
  vector = (int *) mxGetPr((mxArray *) element);
  n = max(dims[0], dims[1]);
}

int FieldComponentsVector::index(int value) {

  for (int i = 0; i < n; i++) {
    if (vector[i] == value) return i;
  }

  return -1;
}

void DetectorSensitivityArrays::initialise(int n_rows, int n_cols) {

  v = (fftw_complex *) fftw_malloc(n_rows * n_cols * sizeof(fftw_complex));
  plan = fftw_plan_dft_2d(n_cols, n_rows, v, v, FFTW_FORWARD, FFTW_MEASURE);

  cm = (complex<double> **) malloc(sizeof(complex<double> *) * n_rows);
  for (int j = 0; j < n_rows; j++) {
    cm[j] = (complex<double> *) malloc(sizeof(complex<double>) * n_cols);
  }
  cerr << "Ex_t_cm has size " << n_rows << "x" << n_cols << endl;
}

DetectorSensitivityArrays::~DetectorSensitivityArrays() {
  fftw_free(v);
  fftw_destroy_plan(plan);
}
