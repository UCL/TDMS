/**
 * @file test_BLi_vs_cubic_interpolation.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests the performance of band-limited interpolation against cubic interpolation
 */
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "globals.h"
#include "interpolation_methods.h"
#include "unit_test_utils.h"

using namespace std;
using tdms_math_constants::DCPI;
using namespace tdms_tests;

// function to interpolate
inline double f_BLi_vs_Cubic(double x) { return sin(2. * DCPI * x) / ((10. * x * x) + 1.); }

/**
 * @brief We will check that BLi gives a better approximation than cubic interpolation.
 *
 * The metric we will use is the norm of the error-vector.
 * As the cell-size becomes larger (across the same domain), BLi should begin to perform better than cubic interpolation.
 *
 * For the following cell sizes, over the range -2 to 2, and for the function f(x) = sin(2\pi x)/(10x^2+1), MATLAB predicts the following errors:
 * Cellsize | (MAX error) BLi :     Cubic       |   (RMSD) BLi    :      Cubic     | %pts BLi better
 * 1.0e-01  |  1.48931677e-03 : 9.96898136e-03  |  2.66761900e-06 : 1.71214373e-04 | 82.05%
 * 5.0e-02  |  8.02870428e-05 : 9.35367215e-04  |  2.31881524e-08 : 8.61216185e-07 | 77.22%
 * 1.0e-02  |  1.75965585e-04 : 1.61196447e-06  |  6.79644463e-08 : 2.41007234e-12 | 0.00%
 * 5.0e-03  |  1.88927787e-04 : 1.00885869e-07  |  7.68604981e-08 : 9.44110925e-15 | 0.25%
 * 1.0e-03  |  1.93280583e-04 : 1.61583136e-10  |  7.99275135e-08 : 2.41912139e-20 | 0.00%
 * 5.0e-04  |  1.93418193e-04 : 1.00993103e-11  |  8.00251513e-08 : 9.44996450e-23 | 0.00%
 * Values computed via the benchmark_test_BLi_vs_cubic.m script.
 *
 * We will test the following:
 * - The RMSD error is_close_or_better than that obtained from MATLAB for both BLi and cubic
 */
TEST_CASE("Benchmark: BLi against cubic interpolation") {
  // number of cell sizes to test at
  const int n_trials = 6;
  // Yee cell sizes to test across
  double cell_sizes[n_trials] = {1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4};
  // Flags for whether or not we expect BLi have superior RMSD
  bool BLi_is_better[n_trials] = {true, true, false, false, false, false};
  // MATLAB errors for BLi
  double MATLAB_BLi_RMSD[n_trials] = {2.66761900e-06, 2.31881524e-08, 6.79644463e-08,
                                      7.68604981e-08, 7.99275135e-08, 8.00251513e-08};
  // MATLAB errors for cubic interp
  double MATLAB_cub_RMSD[n_trials] = {1.71214373e-04, 8.61216185e-07, 2.41007234e-12,
                                      9.44110925e-15, 2.41912139e-20, 9.44996450e-23};
  // coordinate of the LHS of the first Yee cell
  double x_lower = -2.;
  // spatial extent of the domain will be [x_lower, x_lower+extent]
  double extent = 4.;
  // norm errors recorded across runs
  double BLi_RMSD[n_trials], cub_RMSD[n_trials];

  // to record log so it prints out a bit nicer to the screen
  stringstream logging_string;
  logging_string << scientific << setprecision(8);
  logging_string << "\nBenchmarking BLi against Cubic interpolation\n";
  logging_string << "   BLi RMSD      (benchmark)    |   Cubic RMSD     (benchmark)    | Correct "
                    "scheme better? \n";

  for (int trial = 0; trial < n_trials; trial++) {
    // Yee cell "size" to use
    double cellSize = cell_sizes[trial];
    // the number of datapoints that we have, the number of "cells" we have is one less than this
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
      if (i != n_datapts - 1) { cell_centres[i] = x_lower + (((double) i) + 0.5) * cellSize; }
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

    // perform interpolation - this is manual since best_interp_scheme will always want to do BLi
    // cubic interpolation only changes at the first and last cells
    cub_interp[0] = CBFst.interpolate(field_samples);
    for (int i = 1; i < n_datapts - 2; i++) {
      cub_interp[i] = CBMid.interpolate(field_samples, i + 1 - CBMid.number_of_datapoints_to_left);
    }
    cub_interp[n_datapts - 2] =
            CBLst.interpolate(field_samples, n_datapts - 1 - CBLst.number_of_datapoints_to_left);

    // BLi interpolation
    bool check_using_BLi = true;
    for (int i = 0; i < n_datapts - 1; i++) {
      // we checked earlier that BLi should always be available, so this is always a BLi scheme
      InterpolationScheme use_scheme = best_scheme(n_datapts - 1, i);
      // to be on the safe side, check that we have a BLi scheme. BLi schemes are always strictly better than CUBIC_INTERP_MIDDLE, and not equal or worse.
      check_using_BLi = check_using_BLi && use_scheme.is_better_than(CUBIC_INTERP_MIDDLE);
      // interpolate using the cell data we have provided
      // i + 1 - use_scheme.number_of_datapoints_to_left is the index of field_samples to read as the point v[0] in the schemes
      BLi_interp[i] = use_scheme.interpolate(field_samples,
                                             i + 1 - use_scheme.number_of_datapoints_to_left);
    }
    // require that we used BLi for all of these interpolations
    REQUIRE(check_using_BLi);

    // compute RMSD
    BLi_RMSD[trial] =
            relative_mean_square_difference(BLi_interp, exact_field_values, n_datapts - 1);
    cub_RMSD[trial] =
            relative_mean_square_difference(cub_interp, exact_field_values, n_datapts - 1);

    // value-testing commences: the better method should be superior, and the worse method should be worse.
    // assert this is the case for each cellSize we used.
    bool expected_scheme_is_better = ((BLi_RMSD[trial] <= cub_RMSD[trial]) == BLi_is_better[trial]);
    CHECK(expected_scheme_is_better);

    // update the log
    //SPDLOG_INFO("BLi RMSD: {0:.8e} | Cubic RMSD : {1:.8e} | Expected BLi to be better: {2:d}", BLi_RMSD[trial], cub_RMSD[trial], BLi_is_better[trial]);
    logging_string << BLi_RMSD[trial] << " (" << MATLAB_BLi_RMSD[trial] << ") | ";
    logging_string << cub_RMSD[trial] << " (" << MATLAB_cub_RMSD[trial] << ") | ";
    if (expected_scheme_is_better) {
      logging_string << "Yes";
    } else {
      logging_string << "No";
    }
    logging_string << "\n";

    // assert error is better or close
    CHECK(is_close_or_better(BLi_RMSD[trial], MATLAB_BLi_RMSD[trial], 1e-8));
    CHECK(is_close_or_better(cub_RMSD[trial], MATLAB_cub_RMSD[trial], 1e-6));
  }

  SPDLOG_INFO(logging_string.str());
}
