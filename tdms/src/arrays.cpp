#include "arrays.h"

#include <iostream>
#include <spdlog/spdlog.h>
#include <utility>

#include "globals.h"
#include "utils.h"

using namespace std;
using namespace tdms_math_constants;

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
