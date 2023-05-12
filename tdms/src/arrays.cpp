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
      throw std::runtime_error("Have no element " + to_string(c));
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

DispersiveMultiLayer::DispersiveMultiLayer(const mxArray *ptr) {

  if (mxIsEmpty(ptr)) { return; }
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

bool DispersiveMultiLayer::is_dispersive(int K_tot,
                                         double near_zero_tolerance) {
  for (int i = 0; i < K_tot; i++) {
    if (fabs(gamma[i]) > near_zero_tolerance) {
      // non-zero attenuation constant of a Yee cell implies media is dispersive
      return true;
    }
  }
  return false;
}

GratingStructure::GratingStructure(const mxArray *ptr, int I_tot) {

  if (mxIsEmpty(ptr)) { return; }

  auto dims = mxGetDimensions(ptr);
  if (mxGetNumberOfDimensions(ptr) != 2 || dims[0] != 2 ||
      dims[1] != (I_tot + 1)) {
    throw runtime_error("structure should have dimension 2 x (I_tot+1) ");
  }

  matrix = cast_matlab_2D_array((int *) mxGetPr(ptr), 2, I_tot + 1);
}

GratingStructure::~GratingStructure() {
  free_cast_matlab_2D_array(matrix);
  // prevent double free when calling ~Matrix, superclass destructor
  matrix = nullptr;
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

void Pupil::initialise(const mxArray *ptr, int n_rows, int n_cols) {

  if (mxIsEmpty(ptr)) { return; }

  auto dims = (int *) mxGetDimensions(ptr);

  if (mxGetNumberOfDimensions(ptr) != 2 || dims[0] != n_rows ||
      dims[1] != n_cols) {
    throw runtime_error("Pupil has dimension " + to_string(dims[0]) + "x" +
                        to_string(dims[1]) + " but it needed to be " +
                        to_string(n_rows) + "x" + to_string(n_cols));
  }

  matrix = cast_matlab_2D_array(mxGetPr(ptr), n_rows, n_cols);
  this->n_cols = n_cols;
  this->n_rows = n_rows;
}

Pupil::~Pupil() {
  free_cast_matlab_2D_array(matrix);
  matrix = nullptr;
}

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

  auto p = (complex<double> ***) malloc(sizeof(complex<double> **) * n_cols);
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
  tensor.initialise(p, n_cols, n_rows, n_det_modes);
}

void DTilde::initialise(const mxArray *ptr, int n_rows, int n_cols) {

  if (mxIsEmpty(ptr)) { return; }

  assert_is_struct_with_n_fields(ptr, 2, "D_tilde");
  set_component(x, ptr, "Dx_tilde", n_rows, n_cols);
  set_component(y, ptr, "Dy_tilde", n_rows, n_cols);
  n_det_modes =
          mxGetDimensions(ptr_to_nd_array_in(ptr, 3, "Dx_tilde", "D_tilde"))[0];
}

void IncidentField::set_component(Tensor3D<double> &component,
                                  const mxArray *ptr, const std::string &name) {

  if (mxIsEmpty(mxGetField(ptr, 0, name.c_str()))) {
    spdlog::info("{} not present", name);
    return;
  }

  auto element = ptr_to_nd_array_in(ptr, 3, name, "tdfield");
  auto dims = mxGetDimensions(element);
  int N = dims[0], M = dims[1], O = dims[2];
  component.initialise(cast_matlab_3D_array(mxGetPr(element), N, M, O), O, M,
                       N);
  component.is_matlab_initialised = true;

  cerr << "Got tdfield, dims=(" + to_string(N) + "," + to_string(M) + "," +
                  to_string(O) + ")"
       << endl;
}

IncidentField::IncidentField(const mxArray *ptr) {

  assert_is_struct_with_n_fields(ptr, 2, "tdfield");
  set_component(x, ptr, "exi");
  set_component(y, ptr, "eyi");
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

void Vertices::initialise(const mxArray *ptr) {

  auto element = ptr_to_matrix_in(ptr, "vertices", "campssample");
  if (mxIsEmpty(element)) { return; }

  auto dims = mxGetDimensions(element);
  int n_vertices = n_rows = dims[0];
  n_cols = dims[1];

  if (n_cols != 3) {
    throw runtime_error("Second dimension in campssample.vertices must be 3");
  }

  spdlog::info("Found vertices ({0:d} x 3)", n_vertices);
  matrix = cast_matlab_2D_array((int *) mxGetPr(element), n_vertices, n_cols);

  for (int j = 0; j < n_vertices; j++)// decrement index for MATLAB->C indexing
    for (int k = 0; k < n_cols; k++) { matrix[k][j] -= 1; }
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

EHVec::~EHVec() {
  if (has_elements()) {
    for (int i = 0; i < n_rows; i++) fftw_free(matrix[i]);
    free(matrix);
  }
  matrix = nullptr;
}
