/**
 * @file test_interpolation_functions.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests the performance of the interpolation functions, using 1D data
 * mimicing a coordinate axes
 */
#include "interpolation_methods.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <sstream>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "globals.h"
#include "unit_test_utils.h"

using namespace tdms_math_constants;
using Catch::Approx;
using tdms_tests::euclidean;
using tdms_tests::is_close_or_better;
using namespace std;

/**
 * @brief In the case when cubic interpolation is to be used, check that all
 * polynomial fields up to cubic order are interpolated exactly (to within
 * machine error)
 *
 * Checks are run on both the old interp{1,2,3} functions and newer const
 * Interp_scheme instances. Old cubic methods will be redundant upon integration
 * of BLi into the codebase.
 */
TEST_CASE("Cubic interpolation: exact for polynomials of degree <= 3") {
  // equidistant points
  double x[] = {0., 1., 2., 3.};
  // test acceptence tolerance. Allow for FLOP imprecision and rounding errors
  // error should be \approx 4 max(c_i) x^2 __DBL__EPSILON, so order 40 *
  // __DBL__EPSILON__
  double tol = 4e1 * __DBL_EPSILON__;

  // constant field
  double c0 = M_PI;
  // array that will be passed to Interp_scheme::interpolate
  double interp_data[4] = {c0, c0, c0, c0};

  double v1 = c0, v2 = c0, v3 = c0, v4 = c0;
  double v12 = c0, v23 = c0, v34 = c0;

  // check Interp_scheme class method
  CHECK(v12 == Approx(CBFst.interpolate(interp_data)).epsilon(tol));
  CHECK(v23 == Approx(CBMid.interpolate(interp_data)).epsilon(tol));
  CHECK(v34 == Approx(CBLst.interpolate(interp_data)).epsilon(tol));

  // linear
  double c1 = -2.7182818;
  v1 += c1 * x[0];
  interp_data[0] += c1 * x[0];
  v2 += c1 * x[1];
  interp_data[1] += c1 * x[1];
  v3 += c1 * x[2];
  interp_data[2] += c1 * x[2];
  v4 += c1 * x[3];
  interp_data[3] += c1 * x[3];
  v12 += c1 * (x[1] + x[0]) / 2.;
  v23 += c1 * (x[2] + x[1]) / 2.;
  v34 += c1 * (x[3] + x[2]) / 2.;

  // check Interp_scheme class method
  CHECK(v12 == Approx(CBFst.interpolate(interp_data)).epsilon(tol));
  CHECK(v23 == Approx(CBMid.interpolate(interp_data)).epsilon(tol));
  CHECK(v34 == Approx(CBLst.interpolate(interp_data)).epsilon(tol));

  // quadratic
  double c2 = 9.81;
  v1 += c2 * x[0] * x[0];
  interp_data[0] += c2 * x[0] * x[0];
  v2 += c2 * x[1] * x[1];
  interp_data[1] += c2 * x[1] * x[1];
  v3 += c2 * x[2] * x[2];
  interp_data[2] += c2 * x[2] * x[2];
  v4 += c2 * x[3] * x[3];
  interp_data[3] += c2 * x[3] * x[3];
  v12 += c2 * (x[1] + x[0]) * (x[1] + x[0]) / 4.;
  v23 += c2 * (x[2] + x[1]) * (x[2] + x[1]) / 4.;
  v34 += c2 * (x[3] + x[2]) * (x[3] + x[2]) / 4.;

  // check Interp_scheme class method
  CHECK(v12 == Approx(CBFst.interpolate(interp_data)).epsilon(tol));
  CHECK(v23 == Approx(CBMid.interpolate(interp_data)).epsilon(tol));
  CHECK(v34 == Approx(CBLst.interpolate(interp_data)).epsilon(tol));

  // cubic
  double c3 = 4.2;
  v1 += c3 * x[0] * x[0] * x[0];
  interp_data[0] += c3 * x[0] * x[0] * x[0];
  v2 += c3 * x[1] * x[1] * x[1];
  interp_data[1] += c3 * x[1] * x[1] * x[1];
  v3 += c3 * x[2] * x[2] * x[2];
  interp_data[2] += c3 * x[2] * x[2] * x[2];
  v4 += c3 * x[3] * x[3] * x[3];
  interp_data[3] += c3 * x[3] * x[3] * x[3];
  v12 += c3 * (x[1] + x[0]) * (x[1] + x[0]) * (x[1] + x[0]) / 8.;
  v23 += c3 * (x[2] + x[1]) * (x[2] + x[1]) * (x[2] + x[1]) / 8.;
  v34 += c3 * (x[3] + x[2]) * (x[3] + x[2]) * (x[3] + x[2]) / 8.;

  // check Interp_scheme class method
  CHECK(v12 == Approx(CBFst.interpolate(interp_data)).epsilon(tol));
  CHECK(v23 == Approx(CBMid.interpolate(interp_data)).epsilon(tol));
  CHECK(v34 == Approx(CBLst.interpolate(interp_data)).epsilon(tol));
}

/**
 * @brief The hard-coded numerical values for the interpolation constant should
 * all sum to the same value
 *
 * Note - the coefficients are not required to sum to unity!
 */
TEST_CASE("BLi: interpolation-coefficient sums match") {
  /* Tolerance to accept imprecision to
    max. 16 FLOPs implies max discrepency of 8*_DBL_EPSILON__ error
    (8 additions, error in each number max __DBL_EPSILON__).
    tol should be \approx 10 * __DBL_EPSILON__ ,
    after accounting for funny-business multiplying everything by 1
    */
  double tol = 1e1 * __DBL_EPSILON__;
  double coeff_sums[8];
  double a[8] = {1., 1., 1., 1., 1., 1., 1., 1.};

  // fast way to sum the coefficients is to "interpolate" the constant function
  // 1
  coeff_sums[0] = BL0.interpolate(a);
  coeff_sums[1] = BL1.interpolate(a);
  coeff_sums[2] = BL2.interpolate(a);
  coeff_sums[3] = BL3.interpolate(a);
  coeff_sums[4] = BL4.interpolate(a);
  coeff_sums[5] = BL5.interpolate(a);
  coeff_sums[6] = BL6.interpolate(a);
  coeff_sums[7] = BL7.interpolate(a);
  // now check that the entries of coeff_sums are the same
  int n = 8;
  while (--n > 0 && abs(coeff_sums[n] - coeff_sums[0]) < tol)
    ;
  // we only reach the end of the while loop, IE get to n==0, when all elements
  // in the array are the same
  REQUIRE(n == 0);
}

/**
 * @brief We will check that BLi interpolation data gives comparible error to
 * the equivalent functions in MATLAB
 *
 *  * For 100 sample points, we will use BLi to interpolate the following
 * complex-valued function with 100 sample points: real part of sin(2\pi x) imag
 * part of pulse function.
 *
 * Interoplation will then be tested against over the range [0,1], the max
 * element-wise error (by absolute value) will be determined. We will then check
 * that this is of the same order of magnitude as the error produced by
 * MATLAB, 5.35317432e-04.
 *
 * For 100 sample points, we will use BLi to interpolate the following functions
 * with 100 sample points:
 * - The constant function 1    : range 0,1 : max. element-wise error
 * (MATLAB) 2.82944733e-04
 * - sin(2\pi x)                : range 0,1 : max. element-wise error
 * (MATLAB) 2.63468327e-04
 * - pulse function             : range 0,1 : max. element-wise error
 * (MATLAB) 4.87599933e-04
 * - complex function           : range 0,1 : max. element-wise error
 * (MATLAB) 5.35317432e-04
 *
 * MATLAB norm-errors:
 * constant function norm error: 2.81526454e-03
 * sin function norm error:      1.85845330e-03
 * pulse function norm error:    9.65609916e-04
 * complex fn norm error:        2.08892374e-03
 *
 * Error values produced from benchmark_test_interpolation_functions.m
 * The complex function has real part sin(2\pi x) and its imaginary part is the
 * pulse function.
 */

// Hardcode MATLAB errors from benchmark_test_interpolation_functions.m
inline const double ML_const_fn_max_pointwise_error = 2.82944733e-04,
                    ML_sin_max_pointwise_error = 2.63468327e-04,
                    ML_pulse_max_pointwise_error = 4.87599933e-04,
                    ML_complex_fn_max_pointwise_error = 5.35317432e-04;
inline const double ML_cont_fn_norm_error = 2.81526454e-03,
                    ML_sin_norm_error = 1.85845330e-03,
                    ML_pulse_norm_error = 9.65609916e-04,
                    ML_complex_fn_norm_error = 2.08892374e-03;

// Evalutes f(x) = 1, for consistency with test sectioning
inline double constant_1(double x) { return 1.; }

// Evaluates f(x) = sin(2 \pi x)
inline double s2pi(double x) { return sin(2. * DCPI * x); }

/**
 * @brief Evaluates the smooth pulse/ mollifier Kernel function, supported
 * between 0 and 1
 *
 * The smooth mollifier \phi(x) is the function
 * \phi(x) =
 * \begin{cases}
 *  0                   &   when |x| >= 1,
 *  e^{-1/(1-|x|^2)}    &   when |x| < 1
 * \end{cases}
 * This function is compactly supported over the interval [-1,1].
 * By evaluating pulse(x) = \phi(3[2x-1]), we obtain a smooth function with
 * support in [1/3,2/3] \subset [0,1]
 *
 * @param x Point of evaluation
 * @return double Evaluted value
 */
inline double pulse(double x) {
  double absxhat = abs(3. * (2. * x - 1.));
  if (absxhat >= 1) {
    return 0.;
  } else {
    return exp(-1. / (1 - absxhat * absxhat));
  }
}

// Evalutes the complex function we will be using as a benchmark, f(x) =
// sin(2\pi x) + i * pulse(x)
inline complex<double> complex_fn(double x) {
  return s2pi(x) + IMAGINARY_UNIT * pulse(x);
}

TEST_CASE("BLi: MATLAB benchmarking") {
  // setup test logging information
  stringstream logging_string;
  logging_string << scientific << setprecision(8);

  int nSamples = 100;//< number of datapooints in this dimension
  double spacing = 1. / (double) (nSamples - 1);//< spacing between datapoints
  double xi[nSamples];                          //< coordinates of the samples
  double xi5[nSamples - 1];//< coordinates of the interpolation points
  // setup cell centres (xi5) and data sample positions (xi)
  for (int i = 0; i < nSamples - 1; i++) {
    xi[i] = ((double) i) * spacing;
    xi5[i] = xi[i] + spacing / 2.;
  }
  xi[nSamples - 1] = 1.;

  double f_errors[nSamples - 1];//< pointwise interpolation errors
  double max_error;             //< maximum pointwise interpolation error
  double norm_error;            //< norm-error of interpolation
  double MATLAB_max_error;      //< MATLAB max-error benchmark
  double MATLAB_norm_error;     //< MATLAB norm-error benchmark

  SECTION("Real-valued functions") {
    double f_data[nSamples];            //< function data at xi
    double f_exact[nSamples - 1];       //< function exact values at xi5
    double f_interp[nSamples - 1];      //< interpolated values at xi5
    double (*analytic_function)(double);//< the (test) function to interpolate

    SECTION("Constant function") {
      logging_string << "Constant function | ";
      MATLAB_max_error = ML_const_fn_max_pointwise_error;
      MATLAB_norm_error = ML_cont_fn_norm_error;
      analytic_function = &constant_1;
    }
    SECTION("sin(2 pi x)") {
      logging_string << "sin(2 pi x)       | ";
      MATLAB_max_error = ML_sin_max_pointwise_error;
      MATLAB_norm_error = ML_sin_norm_error;
      analytic_function = &s2pi;
    }
    SECTION("Pulse function") {
      logging_string << "pulse function    | ";
      MATLAB_max_error = ML_pulse_max_pointwise_error;
      MATLAB_norm_error = ML_pulse_norm_error;
      analytic_function = &pulse;
    }

    // Setup sample points (xi), Yee cell centres (xi5), and function values at
    // these points (sampled data & exact values)
    for (int i = 0; i < nSamples - 1; i++) {
      f_data[i] = analytic_function(xi[i]);
      f_exact[i] = analytic_function(xi5[i]);
    }
    f_data[nSamples - 1] = analytic_function(xi[nSamples - 1]);
    // perform interpolation
    for (int i = 0; i < nSamples - 1; i++) {
      InterpolationScheme scheme = best_scheme(nSamples, i + 1);
      f_interp[i] = scheme.interpolate(
              f_data, i + 1 - scheme.number_of_datapoints_to_left);
      // Compare interpolated values to the true values
      f_errors[i] = abs(f_exact[i] - f_interp[i]);
    }
  }
  SECTION("Complex-valued functions") {
    complex<double> f_data[nSamples];      // function data at xi
    complex<double> f_exact[nSamples - 1]; // function exact values at xi5
    complex<double> f_interp[nSamples - 1];// interpolated values at xi5

    logging_string << "complex function   | ";
    MATLAB_max_error = ML_complex_fn_max_pointwise_error;
    MATLAB_norm_error = ML_complex_fn_norm_error;

    // Setup sample points (xi), Yee cell centres (xi5), and function values at
    // these points (sampled data & exact values)
    for (int i = 0; i < nSamples - 1; i++) {
      f_data[i] = complex_fn(xi[i]);
      f_exact[i] = complex_fn(xi5[i]);
    }
    f_data[nSamples - 1] = complex_fn(xi[nSamples - 1]);
    // perform interpolation
    for (int i = 0; i < nSamples - 1; i++) {
      InterpolationScheme scheme = best_scheme(nSamples, i + 1);
      f_interp[i] = scheme.interpolate(
              f_data, i + 1 - scheme.number_of_datapoints_to_left);
      // Compare interpolated values to the true values
      f_errors[i] = abs(f_exact[i] - f_interp[i]);
    }
  }

  // Compare maximum pointwise error
  max_error = *max_element(f_errors, f_errors + nSamples - 1);
  REQUIRE(is_close_or_better(max_error, MATLAB_max_error));
  // Compare norm-errors
  norm_error = euclidean(f_errors, nSamples - 2);
  REQUIRE(is_close_or_better(norm_error, MATLAB_norm_error));

  // report test results
  logging_string << "Max ptwise error: " << max_error << " ("
                 << MATLAB_max_error << ") | ";
  logging_string << "Norm error: " << norm_error << " (" << MATLAB_norm_error
                 << ")";
  SPDLOG_INFO(logging_string.str());
}
