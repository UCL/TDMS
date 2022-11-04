/**
 * @file test_DetectorSensitivityArrays.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the DetectorSensitivityArrays class and its subclasses
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"
#include "globals.h"

using namespace std;
using namespace tdms_math_constants;

const double tol = 1e-16;

TEST_CASE("DetectorSensitivityArrays: allocation and deallocation") {
  SPDLOG_INFO("== Testing DetectorSensitivityArrays class");
  // default constructor should set everything to nullptrs
  // destructor uses fftw destroy, which handles nullptrs itself
  DetectorSensitivityArrays empty_dsa;
  bool all_are_nullptrs = (empty_dsa.cm==nullptr) && (empty_dsa.plan==nullptr) && (empty_dsa.v==nullptr);
  REQUIRE(all_are_nullptrs);

  // now let's construct a non-trival object
  DetectorSensitivityArrays dsa;
  const int n_rows = 8, n_cols = 4;
  dsa.initialise(n_rows, n_cols);
  // we should be able to assign complex doubles to the elements of cm now
  for (int i = 0; i < n_rows; i++) {
    for (int j =0; j < n_cols; j++) {
        CHECK_NOTHROW(dsa.cm[i][j] = (double)i + IMAGINARY_UNIT * (double)j);
    }
  }
  // we can call the fftw_plan execution, which should place the 2D FFT into dsa.v
  // simply checking executation is sufficient, as fftw should cover whether the FFT is actually meaningful in what it puts out
  REQUIRE_NOTHROW(fftw_execute(dsa.plan));
  // https://en.wikipedia.org/wiki/Multidimensional_transform#Multidimensional_Fourier_transform <- choose entries accordingly

  // destructor will then clear the memory
}
