#include <iostream>
#include <utility>
#include "arrays.h"
#include "globals.h"
#include "numeric.h"
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

GratingStructure::GratingStructure(const mxArray *ptr, int I_tot) {

  if (mxIsEmpty(ptr)) {
    return;
  }

  auto dims = mxGetDimensions(ptr);
  if (mxGetNumberOfDimensions(ptr) != 2 || dims[0] != 2 || dims[1] != (I_tot + 1)){
    throw runtime_error("structure should have dimension 2 x (I_tot+1) ");
  }

  matrix = cast_matlab_2D_array((int *) mxGetPr(ptr), 2, I_tot + 1);
}

GratingStructure::~GratingStructure() {
  free_cast_matlab_2D_array(matrix);
}

template <typename T>
Vector<T>::Vector(const mxArray *ptr) {
  n = (int)mxGetNumberOfElements(ptr);
  vector = (T*) malloc((unsigned) (n * sizeof(T)));

  auto matlab_ptr = mxGetPr(ptr);
  for (int i = 0; i < n; i++){
    vector[i] = (T) matlab_ptr[i];
  }
}

FrequencyExtractVector::FrequencyExtractVector(const mxArray *ptr, double omega_an) {

  if (mxIsEmpty(ptr)) {
    n = 1;
    vector = (double *) malloc(sizeof(double));
    vector[0] = omega_an / 2. / dcpi;

  } else {
    auto dims = mxGetDimensions(ptr);
    auto n_dims = mxGetNumberOfDimensions(ptr);

    if (n_dims != 2 || !(dims[0] == 1 || dims[1] == 1)){
      throw runtime_error("f_ex_vec should be a vector with N>0 elements");
    }
    cerr << "f_ex_vec has ndims=" << n_dims << " N=" << dims[0] << endl;

    n = max(dims[0], dims[1]);
    vector = (double *) mxGetPr(ptr);
  }
}

void FrequencyVectors::initialise(const mxArray *ptr){

  if (mxIsEmpty(ptr)) {
    return;
  }

  assert_is_struct_with_n_fields(ptr, 2, "f_vec");
  x = Vector<double>(ptr_to_vector_in(ptr, "fx_vec", "f_vec"));
  y = Vector<double>(ptr_to_vector_in(ptr, "fy_vec", "f_vec"));
}

void Pupil::initialise(const mxArray *ptr, int n_rows, int n_cols) {

  if (mxIsEmpty(ptr)){
    return;
  }

  auto dims = (int *)mxGetDimensions(ptr);
  
  if (mxGetNumberOfDimensions(ptr) != 2 || dims[0] != n_rows || dims[1] != n_cols){
    throw runtime_error("Pupil has dimension "+ to_string(dims[0]) + "x"
                        + to_string(dims[1]) + " but it needed to be " +
                        to_string(n_rows) + "x" + to_string(n_cols));
  }

  matrix = cast_matlab_2D_array(mxGetPr(ptr), n_rows, n_cols);
  this->n_cols = n_cols;
  this->n_rows = n_rows;
}

Pupil::~Pupil() {
  free_cast_matlab_2D_array(matrix);
}

template<typename T>
Tensor3D<T>::Tensor3D(T*** tensor, int n_layers, int n_cols, int n_rows){
  this->tensor = tensor;
  this->n_layers = n_layers;
  this->n_cols = n_cols;
  this->n_rows = n_rows;
}

Tensor3D<complex<double>> DTilde::component_in(const mxArray *ptr, const string &name,
                                               int n_rows, int n_cols){

  auto element = ptr_to_nd_array_in(ptr, 3, name, "D_tilde");

  auto dims = (int *)mxGetDimensions(element);
  int n_det_modes = dims[0];

  if (dims[1] != n_rows || dims[2] != n_cols){
    throw runtime_error("D_tilde.{x, y} has final dimensions "+ to_string(dims[1]) + "x"
                        + to_string(dims[2]) + " but it needed to be " +
                        to_string(n_rows) + "x" + to_string(n_cols));
  }

  auto p = (complex<double> ***) malloc(sizeof(complex<double> **) * n_cols);
  for (int j = 0; j < n_cols; j++) {
    p[j] = (complex<double> **) malloc(sizeof(complex<double> *) * n_rows);
    for (int i = 0; i < n_rows; i++) {
      p[j][i] = (complex<double> *) malloc(sizeof(complex<double>) * n_det_modes);
    }
  }

  auto temp_re = cast_matlab_3D_array(mxGetPr(element), dims[0], dims[1], dims[2]);
  auto temp_im = cast_matlab_3D_array(mxGetPi(element), dims[0], dims[1], dims[2]);

  for (int k = 0; k < n_det_modes; k++)
    for (int j = 0; j < n_cols; j++)
      for (int i = 0; i < n_rows; i++) {
        p[j][i][k] = temp_re[j][i][k] + I * temp_im[j][i][k];
      }

  free_cast_matlab_3D_array(temp_re, n_cols);
  free_cast_matlab_3D_array(temp_im, n_cols);

  return {p, n_cols, n_rows, n_det_modes};
}

void DTilde::initialise(const mxArray *ptr, int n_rows, int n_cols) {

  if (mxIsEmpty(ptr)){
    return;
  }

  assert_is_struct_with_n_fields(ptr, 2, "D_tilde");
  x = component_in(ptr, "Dx_tilde", n_rows, n_cols);
  y = component_in(ptr, "Dy_tilde", n_rows, n_cols);
  n_det_modes = mxGetDimensions(ptr_to_nd_array_in(ptr, 3, "Dx_tilde", "D_tilde"))[0];
}

Tensor3D<double> IncidentField::component_in(const mxArray *ptr, const std::string &name){

  if (mxIsEmpty(mxGetField(ptr, 0, name.c_str()))) {
    cerr << name+" not present" << endl;
    return {};
  }

  auto element = ptr_to_nd_array_in(ptr, 3, name, "tdfield");
  auto dims = mxGetDimensions(element);
  int N = dims[0], M = dims[1], O = dims[2];
  auto field = Tensor3D<double>(cast_matlab_3D_array(mxGetPr(element), N, M, O), O, M, N);
  field.is_matlab_initialised = true;

  cerr << "Got tdfield, dims=("+to_string(N)+","+to_string(M)+","+to_string(O)+")" << endl;

  return field;
}

IncidentField::IncidentField(const mxArray *ptr){

  assert_is_struct_with_n_fields(ptr, 2, "tdfield");
  x = component_in(ptr, "exi");
  y = component_in(ptr, "eyi");
}

FieldSample::FieldSample(const mxArray *ptr){

  if (mxIsEmpty(ptr)){
    return;
  }

  assert_is_struct_with_n_fields(ptr, 4, "fieldsample");
  i = Vector<int>(ptr_to_vector_or_empty_in(ptr, "i", "fieldsample"));
  j = Vector<int>(ptr_to_vector_or_empty_in(ptr, "j", "fieldsample"));
  k = Vector<int>(ptr_to_vector_or_empty_in(ptr, "k", "fieldsample"));
  n = Vector<double>(ptr_to_vector_or_empty_in(ptr, "n", "fieldsample"));

  int n_dims = 4;
  if (all_vectors_are_non_empty()){
    int dims[4] = {i.size(), j.size(), k.size(), n.size()};
    mx = mxCreateNumericArray(n_dims, (const mwSize *) dims, mxDOUBLE_CLASS, mxREAL);
    tensor = cast_matlab_4D_array(mxGetPr(mx), i.size(), j.size(), k.size(), n.size());
  } else {
    int dims[4] = {0, 0, 0, 0};
    mx = mxCreateNumericArray(n_dims, (const mwSize *) dims, mxDOUBLE_CLASS, mxREAL);
  }
}

FieldSample::~FieldSample() {
  if (all_vectors_are_non_empty()) {
    free_cast_matlab_4D_array(tensor, k.size(), n.size());
  }
}
