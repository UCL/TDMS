/**
 * @file test_interpolation_determination.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests the logic that determines which interpolation schemes are appropriate.
 */
#include "interpolation_methods.h"

#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

using namespace std;

/**
 * @brief Test whether best_scheme correctly determines the appropriate interpolation scheme to use, given the number of Yee cells either side of cell (i,j,k)
 *
 */
TEST_CASE("best_interp_scheme: correct interpolation chosen") {
    SPDLOG_INFO("===== Testing interpolation scheme selection logic =====");
    int N = 10;
    // should throw out_of_range exception if interpolation is impossible (<3 Yee cells in direction)
    REQUIRE_THROWS_AS(best_scheme(1, 0), out_of_range);
    REQUIRE_THROWS_AS(best_scheme(2, 0), out_of_range);
    REQUIRE_THROWS_AS(best_scheme(2, 1), out_of_range);
    // should throw out_of_range exception if Yee cell if invalid index is requested
    REQUIRE_THROWS_AS(best_scheme(N, -1), out_of_range);
    REQUIRE_THROWS_AS(best_scheme(N, N), out_of_range);

    /* Suppose we have N >= 8 Yee cells in a dimension. The program should determine:
        - cell_id == 0 : Interpolation impossible (checked previously)
        - cell_id == 1,2,3 : Use BAND_LIMITED_(0,1,2) scheme respectively
        - cell_id == 4,...,N-4 : Use BAND_LIMITED_3 scheme
        - cell_id == N-3,N-2,N-1 : Use BAND_LIMITED_(4,5,6) scheme respectively
    */
    N = 10;
    CHECK(best_scheme(N, 0).get_priority() == BAND_LIMITED_0);
    CHECK(best_scheme(N, 1).get_priority() == BAND_LIMITED_1);
    CHECK(best_scheme(N, 2).get_priority() == BAND_LIMITED_2);
    for (int i = 3; i <= N - 4; i++) {
        CHECK(best_scheme(N, i).get_priority() == BAND_LIMITED_3);
    }
    CHECK(best_scheme(N, N - 3).get_priority() == BAND_LIMITED_4);
    CHECK(best_scheme(N, N - 2).get_priority() == BAND_LIMITED_5);
    CHECK(best_scheme(N, N - 1).get_priority() == BAND_LIMITED_6);

    /* If 4 <= N < 8 we can still fall back on cubic interpolation
        - cell_id == 0 : Interpolation impossible (checked previously)
        - cell_id == 1 : Use CUBIC_FIRST
        - cell_id == 2,...,N-2 : Use CUBIC_MIDDLE
        - cell_id == N-1 : Use CUBIC_LAST
    */
    N = 6;
    CHECK(best_scheme(N, 0).get_priority() == CUBIC_INTERP_FIRST);
    for(int i=1; i<=N-2; i++) {CHECK(best_scheme(N, i).get_priority() == CUBIC_INTERP_MIDDLE);}
    CHECK(best_scheme(N, N-1).get_priority() == CUBIC_INTERP_LAST);
}
