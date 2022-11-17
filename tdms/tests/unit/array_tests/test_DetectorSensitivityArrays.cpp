/**
 * @file test_DetectorSensitivityArrays.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the DetectorSensitivityArrays class and its subclasses
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>
#include <fftw3.h>

#include <cmath>

#include "arrays.h"
#include "globals.h"

using namespace std;
using namespace tdms_math_constants;

TEST_CASE("DetectorSensitivityArrays") {

  mxArray *matlab_input;
  DetectorSensitivityArrays dsa;
  const int n_rows = 8, n_cols = 4;

  // default constructor should set everything to nullptrs
  // destructor uses fftw destroy, which handles nullptrs itself
  bool all_are_nullptrs = (dsa.cm==nullptr) && (dsa.plan==nullptr) && (dsa.v==nullptr);
  REQUIRE(all_are_nullptrs);

  // now let's construct a non-trival object
  dsa.initialise(n_rows, n_cols);
  // we should be able to assign complex doubles to the elements of cm now
  // to test the plan executation, we'd need an analytic DFT to hand
  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < n_cols; j++) {
      dsa.cm[i][j] = (double) i + IMAGINARY_UNIT * (double) j;
      dsa.v[j * n_rows + i][0] = dsa.cm[i][j].real();
      dsa.v[j * n_rows + i][1] = dsa.cm[i][j].imag();
    }
  }
  // we can call the fftw_plan execution, which should place the 2D FFT into dsa.v
  // simply checking executation is sufficient, as fftw should cover whether the FFT is actually meaningful in what it puts out
  REQUIRE_NOTHROW(fftw_execute(dsa.plan));
}
