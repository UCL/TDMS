/**
 * @file test_BLi_vs_cubic_interpolation.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests the performance of bandlimited interpolation against cubic
 * interpolation
 */
#include "interpolation_methods.h"

#include <algorithm>
#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "input_flags.h"
#include "unit_test_utils.h"

using namespace std;
using tdms_flags::InterpolationMethod;
using tdms_tests::euclidean;
using tdms_tests::order_of_magnitude;

// function to interpolate
inline double f_BLi_vs_Cubic(double x) { return 1. / ((10. * x * x) + 1.); }

/**
 * @brief We will check that BLi gives a better approximation than cubic
 * interpolation.
 *
 * The metric we will use is the norm of the error-vector.
 * As the cell-size becomes larger (across the same domain), BLi should begin to
 * perform better than cubic interpolation.
 *
 * For the following cell sizes, over the range -2 to 2, and for the function
 * f(x) = 1/(10x^2+1), MATLAB predicts the following errors: Cell size    | BLi
 * err           | Cubic err 0.25 | 1.25953429e-02    | 3.28105674e-02 0.1
 * | 7.10587711e-04    | 5.49601748e-03 0.05 | 7.80325615e-04 | 5.53430587e-04
 * 0.01        | 1.97562511e-03    | 2.06917814e-06 Values computed via the
 * benchmark_test_BLi_vs_cubic.m script.
 *
 * We will test the following:
 * - The norm-error of BLi is lesser/greater than cubic in each case
 * - The order of the magnitude (of BOTH) errors is the same as that obtained
 * from MATLAB (for validation reasons)
 */
TEST_CASE("Benchmark: BLi is better than cubic interpolation") {
  SPDLOG_INFO("===== Benchmarking BLi against cubic interpolation =====");
  // number of cell sizes to test at
  const int n_trials = 4;
  // Yee cell sizes to test across
  double cell_sizes[n_trials] = {0.25, 0.1, 0.05, 0.01};
  // Flags for whether or not we expect BLi to perform better
  bool BLi_is_better[n_trials] = {true, true, false, false};
  // MATLAB errors for BLi
  double MATLAB_BLi_errs[n_trials] = {1.25953429e-02, 7.10587711e-04,
                                      7.80325615e-04, 1.97562511e-03};
  // MATLAB errors for cubic interp
  double MATLAB_cub_errs[n_trials] = {3.28105674e-02, 5.49601748e-03,
                                      5.53430587e-04, 2.06917814e-06};
  // coordinate of the LHS of the first Yee cell
  double x_lower = -2.;
  // spatial extent of the domain will be [x_lower, x_lower+extent]
  double extent = 4.;
  // norm errors recorded across runs
  double BLi_norm_errs[n_trials], cub_norm_errs[n_trials];

  for (int trial = 0; trial < n_trials; trial++) {
    // Yee cell "size" to use
    double cellSize = cell_sizes[trial];
    // the number of datapoints that we have, the number of "cells" we have is
    // one less than this
    int n_datapts = round(extent / cellSize);
    CHECK(n_datapts > 8);// BLi cannot run otherwise

    // the centres of the Yee cells
    double cell_centres[n_datapts - 1];
    // the exact field values at the Yee cell centres
    double exact_field_values[n_datapts - 1];
    // the spatial coordinates of the samples of our field
    double field_positions[n_datapts];
    // the field values at the sample positions
    double field_samples[n_datapts];
    for (int i = 0; i < n_datapts; i++) {
      if (i != n_datapts - 1) {
        cell_centres[i] = x_lower + (((double) i) + 0.5) * cellSize;
      }
      field_positions[i] = x_lower + ((double) i) * cellSize;

      exact_field_values[i] = f_BLi_vs_Cubic(cell_centres[i]);
      field_samples[i] = f_BLi_vs_Cubic(field_positions[i]);
    }

    // the interpolated field values via BLi
    double BLi_interp[n_datapts - 1];
    // the interpolated field values via cubic interp
    double cub_interp[n_datapts - 1];

    // pointwise-interpolation error in BLi
    double BLi_err[n_datapts - 1];
    // pointwise-interpolation error in cubic
    double cub_err[n_datapts - 1];

    // perform interpolation - this is manual since best_interp_scheme will
    // always want to do BLi cubic interpolation only changes at the first and
    // last cells
    cub_interp[0] = CBFst.interpolate(field_samples);
    for (int i = 1; i < n_datapts - 2; i++) {
      cub_interp[i] = CBMid.interpolate(
              field_samples, i + 1 - CBMid.number_of_datapoints_to_left);
    }
    cub_interp[n_datapts - 2] = CBLst.interpolate(
            field_samples, n_datapts - 1 - CBLst.number_of_datapoints_to_left);

    // BLi interpolation
    for (int i = 0; i < n_datapts - 1; i++) {
      // we checked earlier that BLi should always be available, so this is
      // always a BLi scheme
      InterpolationScheme use_scheme =
              best_scheme(n_datapts - 1, i, InterpolationMethod::BandLimited);
      // to be on the safe side, check that we have a BLi scheme. BLi schemes
      // are always strictly better than CUBIC_INTERP_MIDDLE, and not equal or
      // worse.
      CHECK(use_scheme.is_better_than(CUBIC_INTERP_MIDDLE));
      // interpolate using the cell data we have provided
      // i + 1 - use_scheme.number_of_datapoints_to_left is the index of
      // field_samples to read as the point v[0] in the schemes
      BLi_interp[i] = use_scheme.interpolate(
              field_samples, i + 1 - use_scheme.number_of_datapoints_to_left);
    }

    // compare to exact values - again ignore cell 0 for the time being
    for (int i = 0; i < n_datapts - 1; i++) {
      BLi_err[i] = abs(BLi_interp[i] - exact_field_values[i]);
      cub_err[i] = abs(cub_interp[i] - exact_field_values[i]);
    }
    // compute square-norm error
    BLi_norm_errs[trial] = euclidean(BLi_err, n_datapts - 2);
    cub_norm_errs[trial] = euclidean(cub_err, n_datapts - 2);

    // value-testing commences: the better method should be superior, and the
    // worse method should be worse. assert this is the case for each cellSize
    // we used.
    CHECK((BLi_norm_errs[trial] <= cub_norm_errs[trial]) ==
          BLi_is_better[trial]);
    SPDLOG_INFO("BLi err: {0:.8e} | Cubic err : {1:.8e} | Expected BLi to be "
                "better: {2:d}",
                BLi_norm_errs[trial], cub_norm_errs[trial],
                BLi_is_better[trial]);

    // in all cases, assert that the order of magnitude of error is what we
    // expect, if we had used MATLAB
    CHECK(order_of_magnitude(BLi_norm_errs[trial]) ==
          order_of_magnitude(MATLAB_BLi_errs[trial]));
    CHECK(order_of_magnitude(cub_norm_errs[trial]) ==
          order_of_magnitude(MATLAB_cub_errs[trial]));
  }
}
