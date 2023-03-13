/**
 * @file test_numerical_derivative.cpp
 * @brief Tests of the numerical differentiation/FFT functions.
 */
#include "numerical_derivative.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <fftw3.h>
#include <spdlog/spdlog.h>

#include "unit_test_utils.h"

using Catch::Approx;
using tdms_tests::near_zero;

// fftw_complex is typdef to a double[2] - first element is Re, second Im.
const int REAL = 0, IMAG = 1;

/**
 * @brief Test complex number multilication (internally used for numerical
 * derivative).
 */
TEST_CASE("Element-by-element multiplication of array of complex numbers") {

  // setup
  fftw_complex a[3] = {
          {0, 1},
          {1, 1},
          {0, 1},
  };
  fftw_complex b[3] = {
          {0, 1},
          {1, 1},
          {1, -1},
  };
  fftw_complex expected[3] = {
          {-1, 0},// i² = -1
          {0, 2}, // (1+i)² = 2i
          {1, 1}  // i(1-i) = 1+i
  };

  // call
  fftw_complex output[3];
  complex_mult_vec(a, b, output, 3);

  // check the output matches the expected values
  for (int i = 0; i < 3; i++) {
    REQUIRE(output[i][REAL] == expected[i][REAL]);
    REQUIRE(output[i][IMAG] == expected[i][IMAG]);
  }
}

/**
 * @brief Test the numerical derivative function itself.
 *
 * The cosine should have negative sine as a derivative. This test involves a
 * bit of trickery because the normalisation is arbitrary. So the ratio of the
 * expected function (-sinθ) and the output (k dcosθ / dx) should be a constant
 * (k) modulo the places where either function is a zero. So calculate the mean
 * ratio of the expected to the output, and check the ratio ~= it's mean value.
 */
TEST_CASE("Numerical derivative") {

  // array size
  const int NSAMPLES = 64;

  // setup buffers and fft plans
  fftw_complex sampled_cosine[NSAMPLES], output[NSAMPLES], dk[NSAMPLES] = {0.};
  double minus_sine[NSAMPLES];
  fftw_plan pf = fftw_plan_dft_1d(NSAMPLES, sampled_cosine, output,
                                  FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan pb = fftw_plan_dft_1d(NSAMPLES, sampled_cosine, output,
                                  FFTW_BACKWARD, FFTW_ESTIMATE);

  // sample cosθ and its derivative (-sinθ)
  for (int i = 0; i < NSAMPLES; i++) {
    double theta = 100. * i / (double) NSAMPLES * M_PI;
    sampled_cosine[i][REAL] = std::cos(theta);
    minus_sine[i] = -1.0 * std::sin(theta);
  }

  // call
  first_derivative(sampled_cosine, output, dk, NSAMPLES, pf, pb);

  // output has arbitrary normalisation so take the ratio at each sample
  double ratio[NSAMPLES], total = 0;
  int zeroes = 0;
  for (int i = 0; i < NSAMPLES; i++) {

    // skip the zeros of either function - send to NaN
    if (near_zero(output[i][REAL]) || near_zero(minus_sine[i])) {
      ratio[i] = std::numeric_limits<double>::quiet_NaN();
      zeroes++;
    } else {
      ratio[i] = minus_sine[i] / output[i][REAL];
      total += ratio[i];
    }
    spdlog::trace("expected: {} \t calculate: {} \t ratio: {}", minus_sine[i],
                  output[i][REAL], ratio[i]);
  }
  double mean = total / (NSAMPLES - zeroes);

  // every element of the ratio should be close to the mean
  for (int i = 0; i < NSAMPLES; i++)
    if (not std::isnan(ratio[i])) REQUIRE(ratio[i] == Approx(mean));
}
