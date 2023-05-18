/**
 * @file test_interpolation_determination.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests the logic that determines which interpolation schemes are
 * appropriate.
 */
#include "input_flags.h"
#include "interpolation_methods.h"

#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

using namespace std;
using tdms_flags::InterpolationMethod;

/**
 * @brief Test whether best_scheme correctly determines the appropriate
 * interpolation scheme to use, given the number of Yee cells either side of
 * cell (i,j,k)
 *
 */
TEST_CASE("best_interp_scheme: correct interpolation chosen") {
  int N = 10;
  bool all_schemes_correct = true;

  // should throw out_of_range exception if interpolation is impossible (<3 Yee
  // cells in direction)
  SECTION("Too few cells to interpolate") {
    REQUIRE_THROWS_AS(best_scheme(1, 0, InterpolationMethod::BandLimited),
                      out_of_range);
    REQUIRE_THROWS_AS(best_scheme(2, 0, InterpolationMethod::BandLimited),
                      out_of_range);
    REQUIRE_THROWS_AS(best_scheme(2, 1, InterpolationMethod::BandLimited),
                      out_of_range);
  }

  /* Suppose we have N >= 8 Yee cells in a dimension. The program should
    determine:
    - cell_id <  0 : Interpolation impossible
    - cell_id == 0 : BAND_LIMITED_CELL_ZERO
    - cell_id == 1,2,3 : Use BAND_LIMITED_(0,1,2) scheme respectively
    - cell_id == 4,...,N-4 : Use BAND_LIMITED_3 scheme
    - cell_id == N-3,N-2,N-1,N : Use BAND_LIMITED_(4,5,6,7) scheme respectively
    - cell_if == N+1 : Interpolation impossible
  */
  SECTION("Bandlimited interpolation allowed") {
    SPDLOG_INFO("Interpolation scheme selection: bandlimited allowed");

    REQUIRE_THROWS_AS(best_scheme(N, -1, InterpolationMethod::BandLimited),
                      out_of_range);
    all_schemes_correct = all_schemes_correct &&
                          (best_scheme(N, 0, InterpolationMethod::BandLimited)
                                   .get_priority() == BAND_LIMITED_CELL_ZERO);
    all_schemes_correct = all_schemes_correct &&
                          (best_scheme(N, 1, InterpolationMethod::BandLimited)
                                   .get_priority() == BAND_LIMITED_0);
    all_schemes_correct = all_schemes_correct &&
                          (best_scheme(N, 2, InterpolationMethod::BandLimited)
                                   .get_priority() == BAND_LIMITED_1);
    all_schemes_correct = all_schemes_correct &&
                          (best_scheme(N, 3, InterpolationMethod::BandLimited)
                                   .get_priority() == BAND_LIMITED_2);
    for (int i = 4; i <= N - 4; i++) {
      all_schemes_correct = all_schemes_correct &&
                            (best_scheme(N, i, InterpolationMethod::BandLimited)
                                     .get_priority() == BAND_LIMITED_3);
    }
    all_schemes_correct =
            all_schemes_correct &&
            (best_scheme(N, N - 3, InterpolationMethod::BandLimited)
                     .get_priority() == BAND_LIMITED_4);
    all_schemes_correct =
            all_schemes_correct &&
            (best_scheme(N, N - 2, InterpolationMethod::BandLimited)
                     .get_priority() == BAND_LIMITED_5);
    all_schemes_correct =
            all_schemes_correct &&
            (best_scheme(N, N - 1, InterpolationMethod::BandLimited)
                     .get_priority() == BAND_LIMITED_6);
    all_schemes_correct = all_schemes_correct &&
                          (best_scheme(N, N, InterpolationMethod::BandLimited)
                                   .get_priority() == BAND_LIMITED_7);
    REQUIRE(all_schemes_correct);
    REQUIRE_THROWS_AS(best_scheme(N, N + 1, InterpolationMethod::BandLimited),
                      out_of_range);
  }

  /* Suppose we have N >= 8 Yee cells in a dimension. Then if we are restricted
    to the cubic interpolation functions, the program should determine:
    - cell_id <= 0 : Interpolation impossible
    - cell_id == 1 : CUBIC_INTERP_FIRST
    - cell_id == 2,3,...,N-2 : Use CUBIC_INTERP_MIDDLE
    - cell_id == N-1 : CUBIC_INTERP_LAST
    - cell_if >= N : Interpolation impossible
  */
  SECTION("Only cubic interpolation") {
    SPDLOG_INFO("Interpolation scheme selection: bandlimited disallowed");
    InterpolationMethod interpolation_method = InterpolationMethod::Cubic;
    REQUIRE_THROWS_AS(best_scheme(N, 0, interpolation_method), out_of_range);
    all_schemes_correct =
            all_schemes_correct &&
            (best_scheme(N, 1, interpolation_method).get_priority() ==
             CUBIC_INTERP_FIRST);
    for (int i = 2; i <= N - 2; i++) {
      all_schemes_correct =
              all_schemes_correct &&
              (best_scheme(N, i, interpolation_method).get_priority() ==
               CUBIC_INTERP_MIDDLE);
    }
    all_schemes_correct =
            all_schemes_correct &&
            (best_scheme(N, N - 1, interpolation_method).get_priority() ==
             CUBIC_INTERP_LAST);
    REQUIRE(all_schemes_correct);
    REQUIRE_THROWS_AS(best_scheme(N, N, interpolation_method), out_of_range);
  }

  /* If 4 <= N < 8 we fall back on cubic interpolation
    - cell_id == 0 : Interpolation impossible
    - cell_id == 1 : Use CUBIC_FIRST
    - cell_id == 2,...,N-2 : Use CUBIC_MIDDLE
    - cell_id == N-1 : Use CUBIC_LAST
  */
  SECTION("Interpolation scheme selection: forced cubic interpolation") {
    SPDLOG_INFO("Interpolation scheme selection: Forced cubic interpolation");
    N = 6;
    REQUIRE_THROWS_AS(best_scheme(N, 0), out_of_range);
    all_schemes_correct =
            all_schemes_correct &&
            (best_scheme(N, 1).get_priority() == CUBIC_INTERP_FIRST);
    for (int i = 2; i <= N - 2; i++) {
      all_schemes_correct =
              all_schemes_correct &&
              (best_scheme(N, i).get_priority() == CUBIC_INTERP_MIDDLE);
    }
    all_schemes_correct =
            all_schemes_correct &&
            (best_scheme(N, N - 1).get_priority() == CUBIC_INTERP_LAST);
    REQUIRE(all_schemes_correct);
    REQUIRE_THROWS_AS(best_scheme(N, N), out_of_range);
  }
}
