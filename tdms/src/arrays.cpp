#include "arrays.h"

#include <iostream>
#include <spdlog/spdlog.h>
#include <utility>

#include "globals.h"
#include "utils.h"

using namespace std;
using namespace tdms_math_constants;

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
