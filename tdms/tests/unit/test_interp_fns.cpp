# include "catch2/catch_test_macros.hpp"
# include "interpolation_methods.h"

# include <algorithm>
# include <cmath>
# include <iostream>

using namespace std;

/**
 * @brief In the case when cubic interpolation is to be used, check that all polynomial fields up to cubic order are interpolated exactly (to within machine error)
 *
 */
TEST_CASE("cubic_interpolation: cubic interpolation is exact")
{

    // equidistant points
    double x[] = {0., 1., 2., 3.};
    // test acceptence tolerance. Allow for FLOP imprecision and rounding errors
    // error should be \approx 4 max(c_i) x^2 __DBL__EPSILON, so order 40 * __DBL__EPSILON__
    double tol = 4e1 * __DBL_EPSILON__;

    // constant field
    double c0 = 3.1415;
    double v1 = c0, v2 = c0, v3 = c0, v4 = c0;
    double v12 = c0, v23 = c0, v34 = c0;

    CHECK(abs(v12 - interp2(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v23 - interp1(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v34 - interp3(v1, v2, v3, v4)) <= tol);

    // linear
    double c1 = -2.7182818;
    v1 += c1 * x[0];
    v2 += c1 * x[1];
    v3 += c1 * x[2];
    v4 += c1 * x[3];
    v12 += c1 * (x[1] + x[0]) / 2.;
    v23 += c1 * (x[2] + x[1]) / 2.;
    v34 += c1 * (x[3] + x[2]) / 2.;

    CHECK(abs(v12 - interp2(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v23 - interp1(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v34 - interp3(v1, v2, v3, v4)) <= tol);

    // quadratic
    double c2 = 9.81;
    v1 += c2 * x[0] * x[0];
    v2 += c2 * x[1] * x[1];
    v3 += c2 * x[2] * x[2];
    v4 += c2 * x[3] * x[3];
    v12 += c2 * (x[1] + x[0]) * (x[1] + x[0]) / 4.;
    v23 += c2 * (x[2] + x[1]) * (x[2] + x[1]) / 4.;
    v34 += c2 * (x[3] + x[2]) * (x[3] + x[2]) / 4.;

    CHECK(abs(v12 - interp2(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v23 - interp1(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v34 - interp3(v1, v2, v3, v4)) <= tol);

    // cubic
    double c3 = 4.2;
    v1 += c3 * x[0] * x[0] * x[0];
    v2 += c3 * x[1] * x[1] * x[1];
    v3 += c3 * x[2] * x[2] * x[2];
    v4 += c3 * x[3] * x[3] * x[3];
    v12 += c3 * (x[1] + x[0]) * (x[1] + x[0]) * (x[1] + x[0]) / 8.;
    v23 += c3 * (x[2] + x[1]) * (x[2] + x[1]) * (x[2] + x[1]) / 8.;
    v34 += c3 * (x[3] + x[2]) * (x[3] + x[2]) * (x[3] + x[2]) / 8.;

    CHECK(abs(v12 - interp2(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v23 - interp1(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v34 - interp3(v1, v2, v3, v4)) <= tol);
}

/**
 * @brief The hard-coded numerical values for the interpolation constant should all sum to the same value
 *
 * Note - the coefficients are not required to sum to unity!
 */
TEST_CASE("bandlimited_interpolation: coefficient sum")
{

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
    while (--n > 0 && abs(coeff_sums[n] - coeff_sums[0]) < tol);
    // we only reach the end of the while loop, IE get to n==0, when all elements in the array are the same
    REQUIRE(n == 0);
}

/**
 * @brief We will check that BLi interpolation gives comparible error to the equivalent functions in MATLAB
 *
 * For 100 sample points, we will use BLi to interpolate the following functions with 100 sample points:
 * - The constant function 1    : range 0,1         : max. element-wise error (MATLAB) 2.82944733e-04
 * - sin(2\pi x)                : range 0,1         : max. element-wise error (MATLAB) 2.63468327e-04
 * - e^(-1/(2|x-0.5|-1))^2)     : range 0,1         : max. element-wise error (MATLAB) 1.57403846e-03
 *
 * We will then compare the maximum of the pointwise error between the interpolated values and exact values to the same quantity computed via MATLAB's interp function, and determine whether the order of magnitude of the errors is the same.
 *
 * For readability, we break this out into separate test cases for each function.
 * 
 */

// BLi: constant function interpolation
TEST_CASE("bandlimited_interpolation: order of error, constant function")
{

    int nSamples = 100;
    double field_positions[nSamples], Yee_cell_centres[nSamples];

    // constant function variables
    double const_fn_vals[nSamples], const_fn_interp[nSamples], const_fn_exact[nSamples], const_fn_errors[nSamples-1], const_fn_max_error, const_fn_MATLAB_error;

    // input the MATLAB error
    const_fn_MATLAB_error = 2.82944733e-04;

    // setup the sample points, Yee cell centres, and function values (for sampling and exactness)
    for (int i = 0; i < nSamples; i++)
    {
        field_positions[i] = ((double)i / (double)nSamples);
        Yee_cell_centres[i] = field_positions[i] - 1./(2.*(double)nSamples);

        const_fn_vals[i] = 1.; const_fn_exact[i] = 1.;
    }

    // Yee cell 0 has no value "to the left" - this will change with BL_TO_CELL_0 being included.
    // also recall that best_interp_scheme(nSamples, i) returns the scheme that interpolates to the centre of cell i, IE, to position Yee_cell_centres[i].

    // constant function interpolation
    const_fn_interp[1] = best_interp_scheme(nSamples, 1).interpolate(const_fn_vals);
    const_fn_interp[2] = best_interp_scheme(nSamples, 2).interpolate(const_fn_vals);
    const_fn_interp[3] = best_interp_scheme(nSamples, 3).interpolate(const_fn_vals);
    for (int i=4; i<nSamples-4; i++) {
        // need to offset now so that the correct sample points are provided
        const_fn_interp[i] = best_interp_scheme(nSamples, i).interpolate(const_fn_vals, i-4);
    }
    const_fn_interp[nSamples-4] = best_interp_scheme(nSamples, nSamples-4).interpolate(const_fn_vals, nSamples-8);
    const_fn_interp[nSamples-3] = best_interp_scheme(nSamples, nSamples-3).interpolate(const_fn_vals, nSamples-8);
    const_fn_interp[nSamples-2] = best_interp_scheme(nSamples, nSamples-2).interpolate(const_fn_vals, nSamples-8);
    const_fn_interp[nSamples-1] = best_interp_scheme(nSamples, nSamples-1).interpolate(const_fn_vals, nSamples-8);

    // compare interpolated values to the true values. NOTE: index 0 is invalid as we currently don't interpolate to here
    for (int i=0; i<nSamples-1; i++) {
        const_fn_errors[i] = abs( const_fn_exact[i+1] - const_fn_interp[i+1] );
    }

    // get maximum error
    const_fn_max_error = *max_element( const_fn_errors, const_fn_errors + nSamples - 1 );

    // compare O.o.Mag of error
    CHECK(floor(log10(const_fn_max_error)) <= floor(log10(const_fn_MATLAB_error)));
    cout << "const_fn \t | " + to_string(const_fn_max_error) + "\t\t | " + to_string(const_fn_MATLAB_error) + "\n";
}

// BLi: sine function interpolation
TEST_CASE("bandlimited_interpolation: order of error, sine function")
{

    int nSamples = 100;
    double field_positions[nSamples], Yee_cell_centres[nSamples];

    // sin(2\pi x) function variables
    double sin_vals[nSamples], sin_interp[nSamples], sin_exact[nSamples], sin_errors[nSamples - 1], sin_max_error, sin_MATLAB_error;

    // input the MATLAB error
    sin_MATLAB_error = 2.63468327e-04;

    // setup the sample points, Yee cell centres, and function values (for sampling and exactness)
    for (int i = 0; i < nSamples; i++)
    {
        field_positions[i] = ((double)i / (double)nSamples);
        Yee_cell_centres[i] = field_positions[i] - 1. / (2. * (double)nSamples);

        sin_vals[i] = sin(field_positions[i]);
        sin_exact[i] = sin(Yee_cell_centres[i]);
    }

    // Yee cell 0 has no value "to the left" - this will change with BL_TO_CELL_0 being included.
    // also recall that best_interp_scheme(nSamples, i) returns the scheme that interpolates to the centre of cell i, IE, to position Yee_cell_centres[i].

    // sin function interpolation
    sin_interp[1] = best_interp_scheme(nSamples, 1).interpolate(sin_vals);
    sin_interp[2] = best_interp_scheme(nSamples, 2).interpolate(sin_vals);
    sin_interp[3] = best_interp_scheme(nSamples, 3).interpolate(sin_vals);
    for (int i = 4; i < nSamples - 4; i++)
    {
        // need to offset now so that the correct sample points are provided
        sin_interp[i] = best_interp_scheme(nSamples, i).interpolate(sin_vals, i - 4);
    }
    sin_interp[nSamples - 4] = best_interp_scheme(nSamples, nSamples - 4).interpolate(sin_vals, nSamples - 8);
    sin_interp[nSamples - 3] = best_interp_scheme(nSamples, nSamples - 3).interpolate(sin_vals, nSamples - 8);
    sin_interp[nSamples - 2] = best_interp_scheme(nSamples, nSamples - 2).interpolate(sin_vals, nSamples - 8);
    sin_interp[nSamples - 1] = best_interp_scheme(nSamples, nSamples - 1).interpolate(sin_vals, nSamples - 8);

    // compare interpolated values to the true values. NOTE: index 0 is invalid as we currently don't interpolate to here
    for (int i = 0; i < nSamples - 1; i++)
    {
        sin_errors[i] = abs(sin_exact[i + 1] - sin_interp[i + 1]);
    }

    // get maximum error
    sin_max_error = *max_element(sin_errors, sin_errors + nSamples - 1);

    // compare O.o.Mag of error
    CHECK(floor(log10(sin_max_error)) <= floor(log10(sin_MATLAB_error)));
    cout << "sin(2 pi x) \t | " + to_string(sin_max_error) + "\t\t | " + to_string(sin_MATLAB_error) + "\n";
}

// BLi: pulse function interpolation
TEST_CASE("bandlimited_interpolation: order of error, compact pulse")
{

    int nSamples = 100;
    double field_positions[nSamples], Yee_cell_centres[nSamples];

    // pulse function e^(-1/(2|x-0.5|-1))^2) variables
    double pulse_vals[nSamples], pulse_interp[nSamples], pulse_exact[nSamples], pulse_errors[nSamples], pulse_max_error, pulse_MATLAB_error;

    // input the MATLAB error
    pulse_MATLAB_error = 1.57403846e-03;

    // setup the sample points, Yee cell centres, and function values (for sampling and exactness)
    for (int i = 0; i < nSamples; i++)
    {
        field_positions[i] = ((double)i / (double)nSamples);
        Yee_cell_centres[i] = field_positions[i] - 1. / (2. * (double)nSamples);

        if (2. * abs(field_positions[i] - 0.5) >= 1)
        {
            pulse_vals[i] = 0.;
        }
        else
        {
            double x = 2. * abs(field_positions[i] - 0.5);
            pulse_vals[i] = exp(-1. / ((x - 1.) * (x - 1.)));
        }
        if (2. * abs(Yee_cell_centres[i] - 0.5) >= 1)
        {
            pulse_exact[i] = 0.;
        }
        else
        {
            double x = 2. * abs(Yee_cell_centres[i] - 0.5);
            pulse_exact[i] = exp(-1. / ((x - 1.) * (x - 1.)));
        }
    }

    // Yee cell 0 has no value "to the left" - this will change with BL_TO_CELL_0 being included.
    // also recall that best_interp_scheme(nSamples, i) returns the scheme that interpolates to the centre of cell i, IE, to position Yee_cell_centres[i].

    // pulse function interpolation
    pulse_interp[1] = best_interp_scheme(nSamples, 1).interpolate(pulse_vals);
    pulse_interp[2] = best_interp_scheme(nSamples, 2).interpolate(pulse_vals);
    pulse_interp[3] = best_interp_scheme(nSamples, 3).interpolate(pulse_vals);
    for (int i = 4; i < nSamples - 4; i++)
    {
        // need to offset now so that the correct sample points are provided
        pulse_interp[i] = best_interp_scheme(nSamples, i).interpolate(pulse_vals, i - 4);
    }
    pulse_interp[nSamples - 4] = best_interp_scheme(nSamples, nSamples - 4).interpolate(pulse_vals, nSamples - 8);
    pulse_interp[nSamples - 3] = best_interp_scheme(nSamples, nSamples - 3).interpolate(pulse_vals, nSamples - 8);
    pulse_interp[nSamples - 2] = best_interp_scheme(nSamples, nSamples - 2).interpolate(pulse_vals, nSamples - 8);
    pulse_interp[nSamples - 1] = best_interp_scheme(nSamples, nSamples - 1).interpolate(pulse_vals, nSamples - 8);

    // compare interpolated values to the true values. NOTE: index 0 is invalid as we currently don't interpolate to here
    for (int i = 0; i < nSamples - 1; i++)
    {
        pulse_errors[i] = abs(pulse_exact[i + 1] - pulse_interp[i + 1]);
    }

    // get maximum error
    pulse_max_error = *max_element(pulse_errors, pulse_errors + nSamples - 1);

    // compare O.o.Mag of error
    CHECK(floor(log10(pulse_max_error)) <= floor(log10(pulse_MATLAB_error)));
    cout << "pulse \t | " + to_string(pulse_max_error) + "\t\t | " + to_string(pulse_MATLAB_error) + "\n";
}

/**
 * @brief We will check that BLi gives a better approximation than cubic interpolation.
 * 
 * This is not a clear-cut scenario: there are certain function classes (specifically, low-order polynomials, or slowly-varying non-compactly supported functions) in which we expect cubic interpolation to be superior.
 * However, there are other function classes on which we expect BLi to be superior, which in particular encompasses the functional forms of the fields we expect to encounter in the simulations.
 * 
 * As such, we will perform both BLi and cubic interpolation on the function
 * f:[0,1] \mapsto \mathbb{R}, f(x) = [placeholder sin(2\pi x)]
 * We will check that, in each case, the interpolated value given by the BLi scheme is a superior approximation to its cubic counterpart.
 *  
 */
TEST_CASE("Benchmark: BLi is better than cubic interpolation") {
    // number of samples to use
    int nSamples = 100;
    if (nSamples < 8) {
        // BLi cannot run, throw error
        throw runtime_error("nSamples < 8 - cannot interpolate using BLi!");
    }
    // point-spacing
    double spacing = 1./(double) (nSamples-1);

    // coordinates of the "field components" we have data for
    double xi[nSamples]; 
    // coordinates of the "Yee cell" centres - recall that we currently do not interpolate to the centre of Yee cell 0
    double xi5[nSamples-1];

    // field values
    double f_data[nSamples];
    // BLi interpolated values, BLi_interp[i] approximates f(xi5[i]) via BLi
    double BLi_interp[nSamples-1];
    // cubic interpolated values, cubic_interp[i] approximates f(xi5[i]) via cubic interp
    double cubic_interp[nSamples-1];
    // exact function values
    double f_exact[nSamples-1];

    // errors
    double BLi_err[nSamples-1], cubic_err[nSamples-1];

    // populate arrays
    for(int i=0; i<nSamples-1; i++) {
        xi[i] = 0 + i*spacing;
        xi5[i] = xi[i] + spacing/2.;
        f_data[i] = sin(2.*M_PI*xi[i]);
        f_exact[i] = sin(2.*M_PI*xi5[i]);
    }
    xi[nSamples-1] = 1.;
    f_data[nSamples-1] = 0.; // sin(2\pi) = 0

    // perform interpolation
    // cubic interpolation only changes at the first and last cells
    cubic_interp[0] = CBFst.interpolate(f_data, 0);
    for (int i=1; i<nSamples-2; i++) {
        cubic_interp[i] = CBMid.interpolate(f_data, i-1);
    }
    cubic_interp[nSamples-2] = CBLst.interpolate(f_data, nSamples-4); // i = nSamples-2, need to start reading buffer from i-2

    // BLi interpolation is slightly more complex
    // interpolating to xi5[0,1,2] requires us to read f_data from buffer 0
    BLi_interp[0] = BL0.interpolate(f_data);
    BLi_interp[1] = BL1.interpolate(f_data);
    BLi_interp[2] = BL2.interpolate(f_data);
    for (int i = 3; i < nSamples - 4; i++)
    {
        // need to offset now so that the correct sample points are provided
        BLi_interp[i] = BL3.interpolate(f_data, i - 3);
    }
    // running out of data points to the right means we need to fix the buffer again
    BLi_interp[nSamples - 4] = BL4.interpolate(f_data, nSamples - 8);
    BLi_interp[nSamples - 3] = BL5.interpolate(f_data, nSamples - 8);
    BLi_interp[nSamples - 2] = BL6.interpolate(f_data, nSamples - 8);

    // now we should be able to compare the interpolated values to the exact values
    for (int i=0; i<nSamples-1; i++) {
        BLi_err[i] = abs( BLi_interp[i] - f_exact[i] );
        cubic_err[i] = abs( cubic_interp[i] - f_exact[i] );
        // check that BLi beats cubic in all cases
        CHECK( BLi_err[i] < cubic_err[i] );
    }
}