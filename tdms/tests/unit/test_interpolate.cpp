# include <catch2/catch_test_macros.hpp>
# include <interpolation_methods.h>

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
TEST_CASE("determineInterpScheme: correct interpolation chosen") {

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

    /* If 4<=N<=8 we can still fall back on cubic interpolation 
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
    double v1 = c0, v2 = c0, v3 = c0, v4 = c0;
    double v12 = c0, v23 = c0, v34 = c0;
    
    CHECK(abs(v12 - interp2(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v23 - interp1(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v34 - interp3(v1, v2, v3, v4)) <= tol);

    // linear
    double c1 = -2.7182818;
    v1 += c1*x[0]; v2 += c1*x[1]; v3 += c1*x[2]; v4 += c1*x[3];
    v12 += c1*(x[1]+x[0])/2.; v23 += c1*(x[2]+x[1])/2.; v34 += c1*(x[3]+x[2])/2.;

    CHECK(abs(v12 - interp2(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v23 - interp1(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v34 - interp3(v1, v2, v3, v4)) <= tol);

    // quadratic
    double c2 = 9.81;
    v1 += c2*x[0]*x[0]; v2 += c2*x[1]*x[1]; v3 += c2*x[2]*x[2]; v4 += c2*x[3]*x[3];
    v12 += c2 * (x[1] + x[0]) * (x[1] + x[0]) / 4.;
    v23 += c2 * (x[2] + x[1]) * (x[2] + x[1]) / 4.;
    v34 += c2 * (x[3] + x[2]) * (x[3] + x[2]) / 4.;

    CHECK(abs(v12 - interp2(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v23 - interp1(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v34 - interp3(v1, v2, v3, v4)) <= tol);

    // cubic
    double c3 = 4.2;
    v1 += c3*x[0]*x[0]*x[0]; v2 += c3*x[1]*x[1]*x[1]; v3 += c3*x[2]*x[2]*x[2]; v4 += c3*x[3]*x[3]*x[3];
    v12 += c3 * (x[1] + x[0]) * (x[1] + x[0]) * (x[1] + x[0]) / 8.;
    v23 += c3 * (x[2] + x[1]) * (x[2] + x[1]) * (x[2] + x[1]) / 8.;
    v34 += c3 * (x[3] + x[2]) * (x[3] + x[2]) * (x[3] + x[2]) / 8.;

    CHECK(abs(v12 - interp2(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v23 - interp1(v1, v2, v3, v4)) <= tol);
    CHECK(abs(v34 - interp3(v1, v2, v3, v4)) <= tol);
}