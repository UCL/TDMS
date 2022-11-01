/**
 * @file test_XYZTensor3D.cpp
 * @author William Graham
 * @brief Tests allocation, deallocation, and methods for the XYZTensor3D class
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>
#include <complex.h>

#include "arrays.h"
#include "globals.h"

using namespace std;

TEST_CASE("Test XYZTensor allocation") {
  SPDLOG_INFO("== Testing XYZTensor3D allocation/deallocation");
  // dimensions for the test-tensor
  int n_layers = 4, n_cols = 8, n_rows = 16;

  // try creating an XYZTensor3D using the default constructor
  XYZTensor3D<complex<double>> complex_test;
  XYZTensor3D<double> double_test;
  // check that, although this has been declared, its members still point to nullptrs
  // we should be able to index this through AxialDirection
  REQUIRE(complex_test[AxialDirection::X] == nullptr);
  REQUIRE(complex_test[AxialDirection::Y] == nullptr);
  REQUIRE(complex_test[AxialDirection::Z] == nullptr);
  REQUIRE(double_test[AxialDirection::X] == nullptr);
  REQUIRE(double_test[AxialDirection::Y] == nullptr);
  REQUIRE(double_test[AxialDirection::Z] == nullptr);

  // now try to allocate the components
  complex_test.allocate(n_rows, n_cols, n_layers);
  double_test.allocate(n_rows, n_cols, n_layers);
  // we should be able to safely populate each of the XYZ tensors now
  for (char component : {'x', 'y', 'z'}) {
    for (int k = 0; k < n_layers; k++) {
      for (int j = 0; j < n_cols; j++) {
        for (int i = 0; i < n_rows; i++) {
          // shouldn't seg fault when assigning
          // cast type all at once
          REQUIRE_NOTHROW(complex_test[component][k][j][i] = complex<double>(1., 2.));
          REQUIRE_NOTHROW(double_test[component][k][j][i] = 1.);
          // check real/imag part casting
          CHECK(complex_test[component][k][j][i].real() == 1.);
          CHECK(complex_test[component][k][j][i].imag() == 2.);
          // check double assignment was correct
          CHECK(double_test[component][k][j][i] == 1.);
        }
      }
    }
  }
  // at end, destructor should be called as appropriate
  // but this class has no destructor --- it doesn't store the size of its arrays
  // hence, manually deallocate here
  for (char component : {'x', 'y', 'z'}) {
    for (int k = 0; k < n_layers; k++) {
      for (int j = 0; j < n_cols; j++) {
        free(complex_test[component][k][j]);
        free(double_test[component][k][j]);
      }
      free(complex_test[component][k]);
      free(double_test[component][k]);
    }
    free(complex_test[component]);
    free(double_test[component]);
  }
}
