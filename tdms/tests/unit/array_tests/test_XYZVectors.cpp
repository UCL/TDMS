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

// MacOS doesn't allow these to be a const int attribute of the base class, so here they are
const int n_layers = 4, n_cols = 8, n_rows = 16;

void XYZVectorsTest::test_correct_construction() {
  XYZVectors v;
  // check that, although this has been declared, its members still point to nullptrs
  bool are_nullptrs = (v.x == nullptr) && (v.y == nullptr) && (v.z == nullptr);
  REQUIRE(are_nullptrs);
}

void XYZVectorsTest::test_other_methods() {
  XYZVectors v;
  SECTION("set_ptr()") {
    // create some arrays to assign to the members
    double x_vec[n_rows] = {0.}, y_vec[n_cols] = {2., 1.}, z_vec[n_layers] = {4., 3., 2., 1.};
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

}

TEST_CASE("XYZVectors") {
  XYZVectorsTest().run_all_class_tests();
}
