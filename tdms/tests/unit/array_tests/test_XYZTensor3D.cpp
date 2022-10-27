/**
 * @file test_XYZTensor3D.cpp
 * @author William Graham
 * @brief Tests allocation, deallocation, and methods for the XYZTensor3D class
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"
#include "globals.h"

const double tol = 1e-8;

TEST_CASE("Test XYZTensor allocation") {
  SPDLOG_INFO("== Testing XYZTensor3D allocation/deallocation");
  // dimensions for the test-tensor
  int n_layers = 4, n_cols = 8, n_rows = 16;

  // try creating an XYZTensor3D using the default constructor
  XYZTensor3D<double> default_constructor_test;
  // check that, although this has been declared, its members still point to nullptrs
  // we should be able to index this through AxialDirection
  REQUIRE(default_constructor_test[AxialDirection::X] == nullptr);
  REQUIRE(default_constructor_test[AxialDirection::Y] == nullptr);
  REQUIRE(default_constructor_test[AxialDirection::Z] == nullptr);

  // now try to allocate the components
  default_constructor_test.allocate(n_rows, n_cols, n_layers);
  // we should be able to safely populate each of the XYZ tensors now
  for (char component : {'x', 'y', 'z'}) {
    for (int k = 0; k < n_layers; k++) {
      for (int j = 0; j < n_cols; j++) {
        for (int i = 0; i < n_rows; i++) {
          // shouldn't seg fault
          REQUIRE_NOTHROW(default_constructor_test[component][k][j][i] = 0.);
        }
      }
    }
  }
  // at end, destructor should be called as appropriate
  // but this class has no destructor --- it doesn't store the size of its arrays
  // so we should use Tensor3Ds for the components no??
}