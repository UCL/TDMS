# include "catch2/catch_test_macros.hpp"
# include "interpolate.h"

using namespace std;

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

    /* When provided input arguments that do not satisfy any of the above - even cases such as i_l > i_u - should pass without throwing an exception */
    for(int k_l=2; k_l<K-2; k_l++)
        for(int k_u=2; k_u<K-2; k_u++)
            for(int j_l=2; j_l<J-2; j_l++)
                for(int j_u=2; j_u<J-2; j_u++)
                    for(int i_l=2; i_l<I-2; i_l++)
                        for(int i_u=2; i_u++; i_u++)
                            CHECK_NOTHROW(checkInterpolationPoints(i_l, i_u, j_l, j_u, k_l, k_u, I, J, K));
}