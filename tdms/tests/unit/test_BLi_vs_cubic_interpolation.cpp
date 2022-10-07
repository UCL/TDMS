/**
 * @file test_BLi_vs_cubic_interpolation.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests the performance of band-limited interpolation against cubic interpolation
 */
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <iostream>

#include "interpolation_methods.h"

using namespace std;

// Computes the 2-norm of the vector v from buffer start to buffer end
inline double norm(double *v, int end, int start = 0) {
    double norm_val = 0.;
    for (int i = start; i <= end; i++) {
        norm_val += v[i] * v[i];
    }
    return sqrt(norm_val);
}

// function to interpolate
inline double f_BLi_vs_Cubic(double x) {
    return 1. / ((10. * x * x) + 1.);
}

// computes order of magnitude of a double
int orderOfMagnitude(double x) {
    return floor(log10(x));
}

/**
 * @brief We will check that BLi gives a better approximation than cubic interpolation.
 *
 * The metric we will use is the norm of the error-vector.
 * As the cell-size becomes larger (across the same domain), BLi should begin to perform better than cubic interpolation.
 *
 * For the following cell sizes, over the range -2 to 2, and for the function f(x) = 1/(10x^2+1), MATLAB predicts the following errors:
 * Cell size    | BLi err           | Cubic err
 *  0.25        | 1.25953429e-02    | 3.28105674e-02
 *  0.1         | 7.10587711e-04    | 5.49601748e-03
 *  0.05        | 7.80325615e-04    | 5.53430587e-04
 *  0.01        | 1.97562511e-03    | 2.06917814e-06
 * Values computed via the benchmark_test_BLi_vs_cubic.m script.
 *
 * We will test the following:
 * - The norm-error of BLi is lesser/greater than cubic in each case
 * - The order of the magnitude (of BOTH) errors is the same as that obtained from MATLAB (for validation reasons)
 */
TEST_CASE("Benchmark: BLi is better than cubic interpolation") {
    cout << "===== Benchmarking BLi against cubic interpolation =====" << endl;
    // number of cell sizes to test at
    const int n_trials = 4;
    // Yee cell sizes to test across
    double cell_sizes[n_trials] = {0.25, 0.1, 0.05, 0.01};
    // Flags for whether or not we expect BLi to perform better
    bool BLi_is_better[n_trials] = {true, true, false, false};
    // MATLAB errors for BLi
    double MATLAB_BLi_errs[n_trials] = {1.25953429e-02,
                                        7.10587711e-04,
                                        7.80325615e-04,
                                        1.97562511e-03};
    // MATLAB errors for cubic interp
    double MATLAB_cub_errs[n_trials] = {3.28105674e-02,
                                        5.49601748e-03,
                                        5.53430587e-04,
                                        2.06917814e-06};
    // coordinate of the LHS of the first Yee cell
    double x_lower = -2.;
    // spatial extent of the domain will be [x_lower, x_lower+extent]
    double extent = 4.;
    // norm errors recorded across runs
    double BLi_norm_errs[n_trials], cub_norm_errs[n_trials];

    for (int trial = 0; trial < n_trials; trial++) {
        // Yee cell "size" to use
        double cellSize = cell_sizes[trial];
        // the number of Yee cells that we want
        int n_YCs = ceil(extent / cellSize);
        CHECK(n_YCs > 8); // BLi cannot run otherwise

        // the centres of the Yee cells
        double cell_centres[n_YCs];
        // the exact field values at the Yee cell centres
        double exact_field_values[n_YCs];
        // the spatial coordinates of the samples of our field
        double field_positions[n_YCs];
        // the field values at the sample positions
        double field_samples[n_YCs];
        for (int i = 0; i < n_YCs; i++) {
            cell_centres[i] = x_lower + (((double)i) + 0.5) * cellSize;
            field_positions[i] = x_lower + (i + 1) * cellSize;

            exact_field_values[i] = f_BLi_vs_Cubic(cell_centres[i]);
            field_samples[i] = f_BLi_vs_Cubic(field_positions[i]);
        }

        // the interpolated field values via BLi
        double BLi_interp[n_YCs];
        // the interpolated field values via cubic interp
        double cub_interp[n_YCs];

        // pointwise-interpolation error in BLi
        double BLi_err[n_YCs];
        // pointwise-interpolation error in cubic
        double cub_err[n_YCs];

        // perform interpolation - currently DO NOT interpolate to cell 0's centre

        // have to manually do cubic interpolation since best_interp_scheme will always want to do BLi
        // cubic interpolation only changes at the first and last cells
        cub_interp[1] = CBFst.interpolate(field_samples, 1 - CBFst.number_of_datapoints_to_left);
        for (int i = 2; i <= n_YCs - 2; i++) {
            cub_interp[i] = CBMid.interpolate(field_samples, i - CBMid.number_of_datapoints_to_left);
        }
        cub_interp[n_YCs - 1] = CBLst.interpolate(field_samples, n_YCs - 1 - CBLst.number_of_datapoints_to_left);

        // BLi interpolation
        for (int i = 1; i < n_YCs; i++) {
            // we checked earlier that BLi should always be available (CHECK(n_YCs>8)), so this is always a BLi scheme
            interpScheme use_scheme = best_interp_scheme(n_YCs, i);
            // to be on the safe side, check that we have a BLi scheme. BLi schemes are always strictly better than CUBIC_INTERP_MIDDLE, and not equal or worse.
            CHECK(use_scheme.is_better_than(CUBIC_INTERP_MIDDLE));
            // interpolate using the cell data we have provided
            // i - (use_scheme.index + 1) is the index of field_samples to read as the point v[0] in the scheme
            BLi_interp[i] = use_scheme.interpolate(field_samples, i - use_scheme.number_of_datapoints_to_left);
        }

        // compare to exact values - again ignore cell 0 for the time being
        for (int i = 1; i < n_YCs; i++) {
            BLi_err[i] = abs(BLi_interp[i] - exact_field_values[i]);
            cub_err[i] = abs(cub_interp[i] - exact_field_values[i]);
        }
        // compute square-norm error
        BLi_norm_errs[trial] = norm(BLi_err, n_YCs - 1, 1);
        cub_norm_errs[trial] = norm(cub_err, n_YCs - 1, 1);

        // value-testing commences: the better method should be superior, and the worse method should be worse.
        // assert this is the case for each cellSize we used.
        CHECK((BLi_norm_errs[trial] <= cub_norm_errs[trial]) == BLi_is_better[trial]);

        // in all cases, assert that the order of magnitude of error is what we expect, if we had used MATLAB
        CHECK(floor(log10(BLi_norm_errs[trial])) == floor(log10(MATLAB_BLi_errs[trial])));
        CHECK(orderOfMagnitude(cub_norm_errs[trial]) == orderOfMagnitude(MATLAB_cub_errs[trial]));
    }
}
