/**
 * @file test_interpolation_functions.cpp
 * @author William Graham (william.graham@ucl.ac.uk)
 * @brief Tests the result of the interpolation functions.
 */
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>

#include "globals.h"
#include "interpolation_methods.h"

using namespace std;
using namespace TDMS_MATH_CONSTANTS;

/**
 * @brief In the case when cubic interpolation is to be used, check that all polynomial fields up to cubic order are interpolated exactly (to within machine error).
 *
 * Checks are run on both the old interp{1,2,3} functions and newer const Interp_scheme instances. Old cubic methods will be redundant upon integration of BLi into the codebase.
 */
TEST_CASE("test_interpolation_functions: testing that cubic interpolation is exact") {
    // equidistant points
    double x[] = {0., 1., 2., 3.};
    // test acceptence tolerance. Allow for FLOP imprecision and rounding errors
    // error should be \approx 4 max(c_i) x^2 __DBL__EPSILON, so order 40 * __DBL__EPSILON__
    double tol = 4e1 * __DBL_EPSILON__;
    // array that will be passed to Interp_scheme::interpolate
    double interp_data[4];

    // constant field
    double c0 = 3.1415;
    interp_data[0] = c0;
    interp_data[1] = c0;
    interp_data[2] = c0;
    interp_data[3] = c0;

    double v1 = c0, v2 = c0, v3 = c0, v4 = c0;
    double v12 = c0, v23 = c0, v34 = c0;

    // check old interp methods
    CHECK(abs(v12 - interp2(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v23 - interp1(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v34 - interp3(v1, v2, v3, v4)) <= tol);
    // check Interp_scheme class method
    CHECK(abs(v12 - CBFst.interpolate(interp_data)) <= tol);
    CHECK(abs(v23 - CBMid.interpolate(interp_data)) <= tol);
    CHECK(abs(v34 - CBLst.interpolate(interp_data)) <= tol);

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

    // check old interp methods
    CHECK(abs(v12 - interp2(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v23 - interp1(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v34 - interp3(v1, v2, v3, v4)) <= tol);
    // check Interp_scheme class method
    CHECK(abs(v12 - CBFst.interpolate(interp_data)) <= tol);
    CHECK(abs(v23 - CBMid.interpolate(interp_data)) <= tol);
    CHECK(abs(v34 - CBLst.interpolate(interp_data)) <= tol);

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

    // check old interp methods
    CHECK(abs(v12 - interp2(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v23 - interp1(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v34 - interp3(v1, v2, v3, v4)) <= tol);
    // check Interp_scheme class method
    CHECK(abs(v12 - CBFst.interpolate(interp_data)) <= tol);
    CHECK(abs(v23 - CBMid.interpolate(interp_data)) <= tol);
    CHECK(abs(v34 - CBLst.interpolate(interp_data)) <= tol);

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

    // check old interp methods
    CHECK(abs(v12 - interp2(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v23 - interp1(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v34 - interp3(v1, v2, v3, v4)) <= tol);
    // check Interp_scheme class method
    CHECK(abs(v12 - CBFst.interpolate(interp_data)) <= tol);
    CHECK(abs(v23 - CBMid.interpolate(interp_data)) <= tol);
    CHECK(abs(v34 - CBLst.interpolate(interp_data)) <= tol);
}

/**
 * @brief The hard-coded numerical values for the interpolation constant should all sum to the same value
 *
 * Note - the coefficients are not required to sum to unity!
 */
TEST_CASE("bandlimited_interpolation: check that the interpolation constant values all sum to the same value") {
    /* Tolerance to accept imprecision to
    max. 16 FLOPs implies max discrepency of 8*_DBL_EPSILON__ error
    (8 additions, error in each number max __DBL_EPSILON__).
    tol should be \approx 10 * __DBL_EPSILON__ ,
    after accounting for funny-business multiplying everything by 1
    */
    double tol = 1e1 * __DBL_EPSILON__;
    double coeff_sums[8];
    double a[8] = {1., 1., 1., 1., 1., 1., 1., 1.};

    // fast way to sum the coefficients is to "interpolate" the constant function 1
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
    // we only reach the end of the while loop, IE get to n==0, when all elements in the array are the same
    REQUIRE(n == 0);
}

/**
 * @brief We will check that BLi interpolation over real-valued data gives comparible error to the equivalent functions in MATLAB
 *
 * For 100 sample points, we will use BLi to interpolate the following functions with 100 sample points:
 * - The constant function 1    : range 0,1         : max. element-wise error (MATLAB) 2.82944733e-04
 * - sin(2\pi x)                : range 0,1         : max. element-wise error (MATLAB) 2.63468327e-04
 * - pulse function             : range 0,1         : max. element-wise error (MATLAB) 4.87599933e-04
 *
 * We will then compare the maximum of the pointwise error between the interpolated values and exact values to the same quantity computed via MATLAB's interp function, and determine whether the order of magnitude of the errors is the same.
 *
 * For readability, we break this out into separate test cases for each function.
 *
 */

/**
 * @brief Test BLi performance on the constant function.
 * - The constant function 1    : range 0,1         : max. element-wise error (MATLAB) 2.82944733e-04
 */
TEST_CASE("(real-valued) Band limited interpolation: constant function") {
    int nSamples = 100;                         // number of "Yee cells" in this dimension
    double const_fn_data[nSamples];             // function data (only need 8 data points, but simulate 100 cells)
    fill_n(const_fn_data, nSamples, 1);         // initalise to f(x)=1
    double const_fn_interp[nSamples];           // interpolated values
    double const_fn_errors[nSamples];           // error in interpolated values
    double const_fn_max_error;                  // maximum error

    double const_fn_MATLAB_error = 2.82944733e-04;  // the MATLAB error
    
    // Yee cell 0 has no value "to the left" - this will change with BL_TO_CELL_0 being included.
    // Also recall that best_interp_scheme(nSamples, i) returns the scheme that interpolates to the centre of cell i

    // constant function interpolation
    const_fn_interp[1] = best_interp_scheme(nSamples, 1).interpolate(const_fn_data);
    const_fn_interp[2] = best_interp_scheme(nSamples, 2).interpolate(const_fn_data);
    const_fn_interp[3] = best_interp_scheme(nSamples, 3).interpolate(const_fn_data);
    for (int i=4; i<nSamples-4; i++) {
        // need to offset now so that the correct sample points are provided
        const_fn_interp[i] = best_interp_scheme(nSamples, i).interpolate(const_fn_data, i-4);
    }
    const_fn_interp[nSamples-4] = best_interp_scheme(nSamples, nSamples-4).interpolate(const_fn_data, nSamples-8);
    const_fn_interp[nSamples-3] = best_interp_scheme(nSamples, nSamples-3).interpolate(const_fn_data, nSamples-8);
    const_fn_interp[nSamples-2] = best_interp_scheme(nSamples, nSamples-2).interpolate(const_fn_data, nSamples-8);
    const_fn_interp[nSamples-1] = best_interp_scheme(nSamples, nSamples-1).interpolate(const_fn_data, nSamples-8);

    // Compare interpolated values to the true values, which are just f(x)=1.
    // NOTE: cont_fn_interp[0] is uninitialised, since we don't presently interpolate to cell 0's centre, and is thus skipped
    for (int i=0; i<nSamples-1; i++) {
        const_fn_errors[i] = abs( 1. - const_fn_interp[i+1] );
    }
    // get maximum error
    const_fn_max_error = *max_element( const_fn_errors, const_fn_errors + nSamples - 1 );

    // compare O.o.Mag of error - fail if we are orders of magnitude out
    REQUIRE(floor(log10(const_fn_max_error)) <= floor(log10(const_fn_MATLAB_error)));
    // compare absolute error - flag (and fail, but less harshly) if we are doing worse than we expect (but are close)
    CHECK(const_fn_max_error <= const_fn_MATLAB_error);
    // report test results
    cout << fixed << scientific << setprecision(8);
    cout << "BLi: constant function interpolation" << endl;
    cout << "Error: " << const_fn_max_error << " | Benchmark: " << const_fn_MATLAB_error << endl;
}

inline double s2pi(double x) {
    return sin(2. * M_PI * x);
}

/**
 * @brief Test BLi performance on the sine function.
 * - sin(2\pi x)                : range 0,1         : max. element-wise error (MATLAB) 2.63468327e-04
 */
TEST_CASE("(real-valued) Band limited interpolation: sin(2 pi x)") {
    int nSamples = 100;                             // number of "Yee cells" in this dimension
    double spacing = 1. / (double)(nSamples - 1);   // spacing between Yee cell centres
    double xi[nSamples];                            // positions of the "field components"
    double xi5[nSamples];                           // positions of the "Yee cell" centres, xi5[i] = centre of cell i
    double f_data[nSamples];                        // function data at xi
    double f_exact[nSamples];                       // function exact values at xi5
    double f_interp[nSamples];                      // interpolated values at xi5
    double f_errors[nSamples];                      // error at xi5
    double max_error;                               // maximum error across xi5 points

    double sin_MATLAB_error = 2.63468327e-04;       // the MATLAB error

    // setup the sample points, Yee cell centres, and function values (for sampling and exactness)
    for (int i = 0; i < nSamples; i++)
    {
        xi[i] = ((double)i) / ((double)(nSamples - 1));
        xi5[i] = xi[i] - spacing / 2.;

        f_data[i] = s2pi(xi[i]);
        f_exact[i] = s2pi(xi5[i]);
    }

    // Yee cell 0 has no value "to the left" - this will change with BL_TO_CELL_0 being included.
    // also recall that best_interp_scheme(nSamples, i) returns the scheme that interpolates to the centre of cell i, IE, to position xi5[i].

    // sin function interpolation
    f_interp[1] = best_interp_scheme(nSamples, 1).interpolate(f_data);
    f_interp[2] = best_interp_scheme(nSamples, 2).interpolate(f_data);
    f_interp[3] = best_interp_scheme(nSamples, 3).interpolate(f_data);
    for (int i = 4; i < nSamples - 4; i++)
    {
        // need to offset now so that the correct sample points are provided
        f_interp[i] = best_interp_scheme(nSamples, i).interpolate(f_data, i - 4);
    }
    f_interp[nSamples - 4] = best_interp_scheme(nSamples, nSamples - 4).interpolate(f_data, nSamples - 8);
    f_interp[nSamples - 3] = best_interp_scheme(nSamples, nSamples - 3).interpolate(f_data, nSamples - 8);
    f_interp[nSamples - 2] = best_interp_scheme(nSamples, nSamples - 2).interpolate(f_data, nSamples - 8);
    f_interp[nSamples - 1] = best_interp_scheme(nSamples, nSamples - 1).interpolate(f_data, nSamples - 8);

    // Compare interpolated values to the true values. 
    // NOTE: index 0 is invalid as we currently don't interpolate to here
    for (int i = 0; i < nSamples - 1; i++)
    {
        f_errors[i] = abs(f_exact[i + 1] - f_interp[i + 1]);
    }

    // get maximum error
    max_error = *max_element(f_errors, f_errors + nSamples - 1);

    // compare O.o.Mag of error - fail if we are orders of magnitude out
    REQUIRE(floor(log10(max_error)) <= floor(log10(sin_MATLAB_error)));
    // compare absolute error - flag (and fail, but less harshly) if we are doing worse than we expect (but are close)
    CHECK(max_error < sin_MATLAB_error);
    // report test results
    cout << fixed << scientific << setprecision(8);
    cout << "BLi: sine function interpolation" << endl;
    cout << "Error: " << max_error << " | Benchmark: " << sin_MATLAB_error << endl;
}

/**
 * @brief Evaluates the smooth pulse/ mollifier Kernel function, supported between 0 and 1
 *
 * The smooth mollifier \phi(x) is the function
 * \phi(x) =
 * \begin{cases}
 *  0                   &   when |x| >= 1,
 *  e^{-1/(1-|x|^2)}    &   when |x| < 1
 * \end{cases}
 * This function is compactly supported over the interval [-1,1].
 * By evaluating pulse(x) = \phi(3[2x-1]), we obtain a smooth function with support in [1/3,2/3] \subset [0,1]
 *
 * @param x Point of evaluation
 * @return double Evaluted value
 */
inline double pulse(double x) {
    double absxhat = abs(3. * (2.*x - 1.));
    if (absxhat >= 1) {
        return 0.;
    }
    else {
        return exp( -1. / (1 - absxhat*absxhat) );
    }
}

/**
 * @brief Test BLi performance on the compact pulse.
 * - pulse function             : range 0,1         : max. element-wise error (MATLAB) 4.87599933e-04
 */
TEST_CASE("(real-valued) Band limited interpolation: compact pulse") {
    int nSamples = 100;                           // number of "Yee cells" in this dimension
    double spacing = 1. / (double)(nSamples - 1); // spacing between Yee cell centres
    double xi[nSamples];                          // positions of the "field components"
    double xi5[nSamples];                         // positions of the "Yee cell" centres, xi5[i] = centre of cell i
    double f_data[nSamples];                      // constant function data at xi
    double f_exact[nSamples];                     // constant function exact values at xi5
    double f_interp[nSamples];                    // interpolated values at xi5
    double f_errors[nSamples];                    // error at xi5
    double max_error;                             // maximum error across xi5 points

    double pulse_MATLAB_error = 4.87599933e-04;   // the MATLAB error

    // setup the sample points, Yee cell centres, and function values (for sampling and exactness)
    for (int i = 0; i < nSamples; i++)
    {
        xi[i] = ((double)i) / ((double)(nSamples - 1));
        xi5[i] = xi[i] - spacing / 2.;

        f_data[i] = pulse(xi[i]);
        f_exact[i] = pulse(xi5[i]);
    }

    // Yee cell 0 has no value "to the left" - this will change with BL_TO_CELL_0 being included.
    // also recall that best_interp_scheme(nSamples, i) returns the scheme that interpolates to the centre of cell i, IE, to position xi5[i].

    // pulse function interpolation
    f_interp[1] = best_interp_scheme(nSamples, 1).interpolate(f_data);
    f_interp[2] = best_interp_scheme(nSamples, 2).interpolate(f_data);
    f_interp[3] = best_interp_scheme(nSamples, 3).interpolate(f_data);
    for (int i = 4; i < nSamples - 4; i++)
    {
        // need to offset now so that the correct sample points are provided
        f_interp[i] = best_interp_scheme(nSamples, i).interpolate(f_data, i - 4);
    }
    f_interp[nSamples - 4] = best_interp_scheme(nSamples, nSamples - 4).interpolate(f_data, nSamples - 8);
    f_interp[nSamples - 3] = best_interp_scheme(nSamples, nSamples - 3).interpolate(f_data, nSamples - 8);
    f_interp[nSamples - 2] = best_interp_scheme(nSamples, nSamples - 2).interpolate(f_data, nSamples - 8);
    f_interp[nSamples - 1] = best_interp_scheme(nSamples, nSamples - 1).interpolate(f_data, nSamples - 8);

    // compare interpolated values to the true values. NOTE: index 0 is invalid as we currently don't interpolate to here
    for (int i = 0; i < nSamples - 1; i++)
    {
        f_errors[i] = abs(f_exact[i + 1] - f_interp[i + 1]);
    }

    // get maximum error
    max_error = *max_element(f_errors, f_errors + nSamples - 1);
    // compare O.o.Mag of error - fail if we are orders of magnitude out
    REQUIRE(floor(log10(max_error)) <= floor(log10(pulse_MATLAB_error)));
    // compare absolute error - flag (and fail, but less harshly) if we are doing worse than we expect (but are close)
    CHECK(max_error < pulse_MATLAB_error);
    // report test results
    cout << fixed << scientific << setprecision(8);
    cout << "BLi: compact pulse function interpolation" << endl;
    cout << "Error: " << max_error << " | Benchmark: " << pulse_MATLAB_error << endl;
}

/**
 * @brief We will check that BLi interpolation over complex-valued data gives comparible error to the equivalent functions in MATLAB
 *
 * For 100 sample points, we will use BLi to interpolate the following complex-valued function with 100 sample points:
 * real part of sin(2\pi x)
 * imag part of pulse function.
 *
 * Interoplation will then be tested against over the range [0,1], the max element-wise error (by absolute value) will be determined. We will then check that this is of the same order of magnitude as the error produced by MATLAB, 5.35317432e-04.
 */
TEST_CASE("(complex-valued) Band limited interpolation") {
    int nSamples = 100;                            // number of "Yee cells" in this dimension
    double spacing = 1. / (double)(nSamples - 1);  // spacing between Yee cell centres
    double xi[nSamples];                           // positions of the "field components"
    double xi5[nSamples];                          // positions of the "Yee cell" centres, xi5[i] = centre of cell i
    
    complex<double> f_data[nSamples];              // constant function data at xi
    complex<double> f_exact[nSamples];             // constant function exact values at xi5
    complex<double> f_interp[nSamples];            // interpolated values at xi5
    double f_abs_errors[nSamples];                 // error at xi5

    double max_error;                              // maximum error across xi5 points
    double MATLAB_error = 5.35317432e-04;          // the MATLAB error

    // setup the sample points, Yee cell centres, and function values (for sampling and exactness)
    for (int i = 0; i < nSamples; i++) {
        xi[i] = ((double)i) / ((double)(nSamples - 1));
        xi5[i] = xi[i] - spacing / 2.;

        f_data[i] = s2pi(xi[i]) + (pulse(xi[i]) * IMAGINARY_UNIT);
        f_exact[i] = s2pi(xi5[i]) + (pulse(xi5[i]) * IMAGINARY_UNIT);
    }

    // Yee cell 0 has no value "to the left" - this will change with BL_TO_CELL_0 being included.
    // also recall that best_interp_scheme(nSamples, i) returns the scheme that interpolates to the centre of cell i, IE, to position xi5[i].

    // pulse function interpolation
    f_interp[1] = best_interp_scheme(nSamples, 1).interpolate(f_data);
    f_interp[2] = best_interp_scheme(nSamples, 2).interpolate(f_data);
    f_interp[3] = best_interp_scheme(nSamples, 3).interpolate(f_data);
    for (int i = 4; i < nSamples - 4; i++) {
        // need to offset now so that the correct sample points are provided
        f_interp[i] = best_interp_scheme(nSamples, i).interpolate(f_data, i - 4);
    }
    f_interp[nSamples - 4] = best_interp_scheme(nSamples, nSamples - 4).interpolate(f_data, nSamples - 8);
    f_interp[nSamples - 3] = best_interp_scheme(nSamples, nSamples - 3).interpolate(f_data, nSamples - 8);
    f_interp[nSamples - 2] = best_interp_scheme(nSamples, nSamples - 2).interpolate(f_data, nSamples - 8);
    f_interp[nSamples - 1] = best_interp_scheme(nSamples, nSamples - 1).interpolate(f_data, nSamples - 8);

    // compare interpolated values to the true values. NOTE: index 0 is invalid as we currently don't interpolate to here
    for (int i = 0; i < nSamples - 1; i++) {
        f_abs_errors[i] = abs(f_exact[i + 1] - f_interp[i + 1]);
    }

    // get maximum error
    max_error = *max_element(f_abs_errors, f_abs_errors + nSamples - 1);
    // compare O.o.Mag of error - fail if we are orders of magnitude out
    REQUIRE(floor(log10(max_error)) <= floor(log10(MATLAB_error)));
    // compare absolute error - flag (and fail, but less harshly) if we are doing worse than we expect (but are close)
    CHECK(max_error < MATLAB_error);
    // report test results
    cout << fixed << scientific << setprecision(8);
    cout << "BLi: complex function interpolation" << endl;
    cout << "Error: " << max_error << " | Benchmark: " << MATLAB_error << endl;
}
