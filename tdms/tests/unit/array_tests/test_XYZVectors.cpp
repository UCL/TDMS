/**
 * @file test_XYZVectors.cpp
 * @author William Graham
 * @brief Tests allocation, deallocation, and methods for the XYZVectors class
 */
#include <stdexcept>

#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays/xyz_vector.h"
#include "globals.h"
#include "unit_test_utils.h"

TEST_CASE("XYZVector") {
  XYZVector xyz;
  const int n_x = 4, n_y = 8, n_z = 16;

  // Confirm assignment is possible via the [] operator
  SECTION("Access via []") {
    xyz['x'].resize(n_x);
    xyz['y'].resize(n_y);
    xyz['z'].resize(n_z);
    // [] operator throws if an unrecognised character is passed
    REQUIRE_THROWS_AS(xyz['t'], std::runtime_error);
  }

  SECTION("all_elements_less_than()") {
    xyz['x'] = {0.5, 1.0, 1.5, 2.0};
    xyz['y'] = {0.25, 0.5, 0.75};
    xyz['z'] = {0.1, 0.2, 0.3};
    REQUIRE(!xyz.all_elements_less_than(1.));
    REQUIRE(!xyz.all_elements_less_than(1., AxialDirection::X));
    REQUIRE(xyz.all_elements_less_than(1., AxialDirection::Y));
    REQUIRE(xyz.all_elements_less_than(1., AxialDirection::Z));
  }
}
