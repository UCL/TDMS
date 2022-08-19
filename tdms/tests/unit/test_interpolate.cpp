# include "catch2/catch_test_macros.hpp"
# include "interpolation_methods.h"

using namespace std;

/**
 * @brief Test whether checkInterpolationPoints throws exceptions in the correct cases.
 * 
 * This function checks whether it is possible to use cubic interpolation over a given index range.
 */
TEST_CASE("checkInterpolationPoints: exceptions thrown") {

    // setup some fake field dimensions
    int I = 6, J = 7, K = 8;

    // 2D simulations have J = 0 - no width in this dimension
    // checkInterpolationPoints should not throw an exception if we provide j_l = 2 > j_u = -2, as is done for a 2D simulation (J=0), and request interpolation over the maximum i,k range
    CHECK_NOTHROW(checkInterpolationPoints(2, I - 2, 2, -2, 2, K - 2, I, 0, K));

    // conversely, it should throw an error if we try to provide j_l = 0 = j_u, which one may percieve are sensible choices, in this case
    CHECK_THROWS_AS(checkInterpolationPoints(2, I - 2, 2, 0, 0, K - 2, I, 0, K), runtime_error);

    /* 3D simulation checks
    checkInterpolationPoints should throw errors in each of the following cases:
    - i_l < 2 with text Interpolation error: i_l too small
    - j_l < 2 with text Interpolation error: j_l too small
    - k_l < 2 with text Interpolation error: k_l too small
    - i_u > I-2 with text Interpolation error: i_u too large
    - j_u > J-2 with text Interpolation error: j_u too large
    - k_u > K-2 with text Interpolation error: k_u too large */
    CHECK_THROWS_AS(checkInterpolationPoints(1, I - 2, 2, J - 2, 2, K - 2, I, J, K), runtime_error);
    CHECK_THROWS_AS(checkInterpolationPoints(2, I - 1, 2, J - 2, 2, K - 2, I, J, K), runtime_error);
    CHECK_THROWS_AS(checkInterpolationPoints(2, I - 2, 0, J - 2, 2, K - 2, I, J, K), runtime_error);
    CHECK_THROWS_AS(checkInterpolationPoints(2, I - 2, 2, J - 0, 2, K - 2, I, J, K), runtime_error);
    CHECK_THROWS_AS(checkInterpolationPoints(2, I - 1, 2, J - 2, -1, K - 2, I, J, K), runtime_error);
    CHECK_THROWS_AS(checkInterpolationPoints(2, I - 1, 2, J - 2, 2, K + 1, I, J, K), runtime_error);
}

TEST_CASE("checkInterpolationPoints: check valid inputs") {

    // setup some fake field dimensions
    int I = 6, J = 7, K = 8;

    /* checkInterpolationPoints should throw errors in each of the following cases:
    - i_l < 2 with text Interpolation error: i_l too small
    - j_l < 2 with text Interpolation error: j_l too small
    - k_l < 2 with text Interpolation error: k_l too small
    - i_u > I-2 with text Interpolation error: i_u too large
    - j_u > J-2 with text Interpolation error: j_u too large
    - k_u > K-2 with text Interpolation error: k_u too large
    When provided input arguments that do not satisfy any of the above - even cases such as i_l > i_u - should pass without throwing an exception */
    for (int k_l = 2; k_l < K - 2; k_l++)
        for (int k_u = 2; k_u < K - 2; k_u++)
            for (int j_l = 2; j_l < J - 2; j_l++)
                for (int j_u = 2; j_u < J - 2; j_u++)
                    for (int i_l = 2; i_l < I - 2; i_l++)
                        for (int i_u = 2; i_u < I - 2; i_u++)
                            CHECK_NOTHROW(checkInterpolationPoints(i_l, i_u, j_l, j_u, k_l, k_u, I, J, K));
}

/**
 * @brief Test whether determineInterpScheme correctly determines the appropriate interpolation scheme to use, given the number of Yee cells either side of cell (i,j,k)
 * 
 */
TEST_CASE("best_interp_scheme: correct interpolation chosen") {

    int N = 10;
    // should throw out_of_range exception if interpolation is impossible (<4 Yee cells in direction)
    REQUIRE_THROWS_AS(best_interp_scheme(3, 2), out_of_range);
    // should throw out_of_range exception if Yee cell of invalid index is requested
    REQUIRE_THROWS_AS(best_interp_scheme(N,0), out_of_range);
    REQUIRE_THROWS_AS(best_interp_scheme(N,N), out_of_range);

    /* Suppose we have N >= 8 Yee cells in a dimension. The program should determine:
        - cell_id == 0 : Interpolation impossible (checked previously)
        - cell_id == 1,2,3 : Use BAND_LIMITED_(0,1,2) scheme respectively
        - cell_id == 4,...,N-4 : Use BAND_LIMITED_3 scheme
        - cell_id == N-3,N-2,N-1 : Use BAND_LIMITED_(4,5,6) scheme respectively
    */
    N = 10;
    CHECK(best_interp_scheme(N, 1).get_priority() == BAND_LIMITED_0);
    CHECK(best_interp_scheme(N, 2).get_priority() == BAND_LIMITED_1);
    CHECK(best_interp_scheme(N, 3).get_priority() == BAND_LIMITED_2);
    for(int i=4; i<=N-4; i++) {CHECK(best_interp_scheme(N, i).get_priority() == BAND_LIMITED_3);}
    CHECK(best_interp_scheme(N, N-3).get_priority() == BAND_LIMITED_4);
    CHECK(best_interp_scheme(N, N-2).get_priority() == BAND_LIMITED_5);
    CHECK(best_interp_scheme(N, N-1).get_priority() == BAND_LIMITED_6);

    /* If 4 <= N < 8 we can still fall back on cubic interpolation 
        - cell_id == 0 : Interpolation impossible (checked previously)
        - cell_id == 1 : Use CUBIC_FIRST
        - cell_id == 2,...,N-2 : Use CUBIC_MIDDLE
        - cell_id == N-1 : Use CUBIC_LAST
    */
    N = 7;
    CHECK(best_interp_scheme(N, 1).get_priority() == CUBIC_INTERP_FIRST);
    for(int i=2; i<=N-2; i++) {CHECK(best_interp_scheme(N, i).get_priority() == CUBIC_INTERP_MIDDLE);}
    CHECK(best_interp_scheme(N, N-1).get_priority() == CUBIC_INTERP_LAST);
}

/**
 * @brief In the case when cubic interpolation is to be used, check that all polynomial fields up to cubic order are interpolated exactly (to within machine error)
 * 
 */
TEST_CASE("interp: cubic interpolation is exact") {

    // equidistant points
    double x[] = {0.,1.,2.,3.};
    // test acceptence tolerance. Allow for FLOP imprecision and rounding errors
    // error should be \approx 4 max(c_i) x^2 __DBL__EPSILON, so order 40 * __DBL__EPSILON__
    double tol = 4e1 * __DBL_EPSILON__;

    // constant field
    double c0 = 3.1415;
    // data points to use for interpolation
    double v[8];
    v[0] = c0, v[1] = c0, v[2] = c0, v[3] = c0;
    // exact values (interpolate f(x) = c0)
    double v12 = c0, v23 = c0, v34 = c0;
    
    CHECK(abs(v12 - CBFst.interpolate(v) <= tol));
    CHECK(abs(v23 - CBMid.interpolate(v) <= tol));
    CHECK(abs(v34 - CBLst.interpolate(v) <= tol));

    // linear, create coefficient
    double c1 = -2.7182818;
    // update interpolation data points
    v[0] += c1*x[0]; v[1] += c1*x[1]; v[2] += c1*x[2]; v[3] += c1*x[3];
    // update exact values (f(x) = c0 + c1*x)
    v12 += c1*(x[1]+x[0])/2.; v23 += c1*(x[2]+x[1])/2.; v34 += c1*(x[3]+x[2])/2.;

    CHECK(abs(v12 - CBFst.interpolate(v)) <= tol);
    CHECK(abs(v23 - CBMid.interpolate(v)) <= tol);
    CHECK(abs(v34 - CBLst.interpolate(v)) <= tol);

    // quadratic, create coefficient
    double c2 = 9.81;
    // update interpolation data points
    v[0] += c2*x[0]*x[0]; v[1] += c2*x[1]*x[1]; v[2] += c2*x[2]*x[2]; v[3] += c2*x[3]*x[3];
    // update exact values (f(x) = c0 + c1*x + c2*x^2)
    v12 += c2 * (x[1] + x[0]) * (x[1] + x[0]) / 4.;
    v23 += c2 * (x[2] + x[1]) * (x[2] + x[1]) / 4.;
    v34 += c2 * (x[3] + x[2]) * (x[3] + x[2]) / 4.;

    CHECK(abs(v12 - CBFst.interpolate(v)) <= tol);
    CHECK(abs(v23 - CBMid.interpolate(v)) <= tol);
    CHECK(abs(v34 - CBLst.interpolate(v)) <= tol);

    // cubic, create coefficient
    double c3 = 4.2;
    // update interpolation data points
    v[0] += c3*x[0]*x[0]*x[0]; v[1] += c3*x[1]*x[1]*x[1]; v[2] += c3*x[2]*x[2]*x[2]; v[3] += c3*x[3]*x[3]*x[3];
    v12 += c3 * (x[1] + x[0]) * (x[1] + x[0]) * (x[1] + x[0]) / 8.;
    v23 += c3 * (x[2] + x[1]) * (x[2] + x[1]) * (x[2] + x[1]) / 8.;
    v34 += c3 * (x[3] + x[2]) * (x[3] + x[2]) * (x[3] + x[2]) / 8.;

    CHECK(abs(v12 - CBFst.interpolate(v)) <= tol);
    CHECK(abs(v23 - CBMid.interpolate(v)) <= tol);
    CHECK(abs(v34 - CBLst.interpolate(v)) <= tol);
}

/**
 * @brief The hard-coded numerical values for the interpolation constant should all sum to the same value
 *
 * Note - the coefficients are not required to sum to unity!
 */ /*
TEST_CASE("bandlimited_interpolation: coefficient sum")
{

    double tol = __DBL_EPSILON__;
    double coeff_sums[8];
    double a[8] = {1., 1., 1., 1., 1., 1., 1., 1.};

    // fast way to sum the coefficients is to "interpolate" the constant function
    for (int i = 0; i < 8; i++)
    {
        coeff_sums[i] = bandlimited_interpolation(i, a);
    }
    // now check that the entries of coeff_sums are the same
    int n = 8;
    while (--n > 0 && abs(a[n] - a[0]) < tol)
        ;
    // we only reach the end of the while loop, IE get to n==0, when all elements in the array are the same
    REQUIRE(n == 0);
}

/**
 * @brief We will check that BLi interpolation gives comparible error to the equivalent functions in MATLAB
 *
 * For 100 sample points, we will use BLi to interpolate the following functions with 100 sample points:
 * - The constant function 1    : range 0,1
 * - sin(x)                     : range 0-2\pi
 *
 * We will then compare the maximum of the pointwise error between the interpolated values and exact values to the same quantity computed via MATLAB's interp function, off of which the bandlimited_interpolation method is based.
 *
 */ /*
TEST_CASE("bandlimited_interpolation: order of error")
{

    int N = 100;
    double x[N], x5[N];
    for (int i = 0; i < N; i++)
    {
        x[i] = ((double)i / (double)N) * 2. * M_PI;
        x5[i] = (2. * i + 1.) * M_PI / (double)N;
    }

    double sin_vals[N], interp_vals[N], diffs[N];
    double exact_vals[N];

    // sin_vals[i] = sin(2pi i/N)
    for (int i = 0; i < N; i++)
    {
        sin_vals[i] = sin(x[i]);
        exact_vals[i] = sin(x5[i]);
    }

    // now try interpolating...
    interp_vals[0] = bandlimited_interpolation(0, sin_vals);
    interp_vals[1] = bandlimited_interpolation(1, sin_vals);
    interp_vals[2] = bandlimited_interpolation(2, sin_vals);
    for (int i = 3; i < N - 4; i++)
    {
        interp_vals[i] = bandlimited_interpolation(3, sin_vals, i - 3);
    }
    interp_vals[N - 4] = bandlimited_interpolation(4, sin_vals, N - 8 - 1);
    interp_vals[N - 3] = bandlimited_interpolation(5, sin_vals, N - 8 - 1);
    interp_vals[N - 2] = bandlimited_interpolation(6, sin_vals, N - 8 - 1);
    interp_vals[N - 1] = bandlimited_interpolation(7, sin_vals, N - 8 - 1);

    // print maximum difference?
    cout << "index | difference | interp_value \n";
    for (int i = 0; i < N; i++)
    {
        diffs[i] = abs(sin_vals[i] - interp_vals[i]);
        cout << to_string(i) + " | " + to_string(diffs[i]) + " | " + to_string(interp_vals[i]) + "\n";
    }
    cout << "Maximum difference : " + to_string(*max_element(diffs, diffs + N)) + "\n";
}

/**
 * @brief Test whether the implimentation of band-limited interpolation is performing correctly.
 *
 * We will attempt to interpolate the function
 * f(x) =
 * which can be interpolated exactly by our band-limited interpolation scheme.
 */