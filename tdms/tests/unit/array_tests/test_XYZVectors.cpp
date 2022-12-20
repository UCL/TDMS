/**
 * @file test_XYZVectors.cpp
 * @author William Graham
 * @brief Tests allocation, deallocation, and methods for the XYZVectors class
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>
#include <numeric>

#include "arrays.h"
#include "array_test_class.h"
#include "globals.h"
#include "unit_test_utils.h"

using std::accumulate;
using tdms_tests::is_close;

void XYZVectorsTest::test_correct_construction() {
  XYZVectors v;
  // check that, although this has been declared, its members still point to nullptrs
  bool are_nullptrs = (v.x == nullptr) && (v.y == nullptr) && (v.z == nullptr);
  REQUIRE(are_nullptrs);
}

void XYZVectorsTest::test_other_methods() {
  XYZVectors v;
  // create some arrays to assign to the members
  double x_vec[n_rows] = {0.}, y_vec[n_cols] = {2., 1.}, z_vec[n_layers] = {1., 2., 3., 4.};

  SECTION("set_ptr()") {
    // now let's try assigning
    v.set_ptr(AxialDirection::X, x_vec);
    v.set_ptr(AxialDirection::Y, y_vec);
    v.set_ptr(AxialDirection::Z, z_vec);
    bool are_not_nullptrs = (v.x != nullptr) && (v.y != nullptr) && (v.z != nullptr);
    REQUIRE(are_not_nullptrs);

    // and the components are what we expect
    double x_tot = accumulate(v.x, v.x + n_rows, 0.);
    double y_tot = accumulate(v.y, v.y + n_cols, 0.);
    double z_tot = accumulate(v.z, v.z + n_layers, 0.);
    bool correct_elements = true;
    for (int i = 0; i < n_rows; i++) {
      correct_elements = correct_elements && (is_close(v.x[i], x_vec[i]));
    }
    for (int j = 0; j < n_cols; j++) {
      correct_elements = correct_elements && (is_close(v.y[j], y_vec[j]));
    }
    for (int k = 0; k < n_layers; k++) {
      correct_elements = correct_elements && (is_close(v.z[k], z_vec[k]));
    }
    REQUIRE(correct_elements);

    // in theory, we can also swap the components by reassigning the pointers
    v.set_ptr('x', z_vec);
    bool correct_swapped_elements = true;
    for (int k = 0; k < n_layers; k++) {
      correct_elements = correct_elements && (is_close(v.x[k], z_vec[k]));
    }
    REQUIRE(correct_elements);
  }
  SECTION("all_elements_less_than()") {
    v.set_ptr(AxialDirection::X, x_vec);
    v.set_ptr(AxialDirection::Y, y_vec);
    v.set_ptr(AxialDirection::Z, z_vec);
    // v.x: max element is 0, v.y: max is 2, v.z: max is 4
    REQUIRE(v.all_elements_less_than(5., 1, 2, 4));
    // v.x has no elements less than 0
    REQUIRE(!v.all_elements_less_than(-1, 1, AxialDirection::X));
    // v.y has all elements but the first less than 2
    REQUIRE(v.all_elements_less_than(2., 1, AxialDirection::Y, 1));
    // v.z has all elements except the last less than 4
    REQUIRE(v.all_elements_less_than(4., 3, AxialDirection::Z));
  }
}

TEST_CASE("XYZVectors") {
  XYZVectorsTest().run_all_class_tests();
}
