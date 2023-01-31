/**
 * @file test_XYZTensor3D.cpp
 * @author William Graham
 * @brief Tests allocation, deallocation, and methods for the XYZTensor3D class
 */
#include <catch2/catch_test_macros.hpp>
#include <complex>
#include <spdlog/spdlog.h>

#include "array_test_class.h"
#include "arrays.h"
#include "globals.h"
#include "unit_test_utils.h"

using namespace std;
using tdms_tests::is_close;

void XYZTensor3DTest::test_correct_construction() {
  // try creating an XYZTensor3D using the default constructor
  XYZTensor3D<complex<double>> xyzt3d;
  // check that, although this has been declared, its members still point to nullptrs
  // we should be able to index this through AxialDirection
  bool all_nullptrs = (xyzt3d[AxialDirection::X] == nullptr) &&
                      (xyzt3d[AxialDirection::Y] == nullptr) &&
                      (xyzt3d[AxialDirection::Z] == nullptr);
  REQUIRE(all_nullptrs);
}

void XYZTensor3DTest::test_other_methods() {
  SECTION("allocate()") {
    XYZTensor3D<complex<double>> xyzt3d;
    // now try to allocate the components
    xyzt3d.allocate(n_rows, n_cols, n_layers);
    // we should be able to safely populate each of the XYZ tensors now
    bool assignments_succeeded = true;
    for (char component : {'x', 'y', 'z'}) {
      for (int k = 0; k < n_layers; k++) {
        for (int j = 0; j < n_cols; j++) {
          for (int i = 0; i < n_rows; i++) {
            // shouldn't seg fault when assigning
            // cast type all at once
            xyzt3d[component][k][j][i] = complex<double>((double) i * j, (double) i * k);
            // check real/imag part casting
            assignments_succeeded = assignments_succeeded &&
                                    is_close(xyzt3d[component][k][j][i].real(), (double) i * j) &&
                                    is_close(xyzt3d[component][k][j][i].imag(), (double) i * k);
          }
        }
      }
    }
    REQUIRE(assignments_succeeded);
    // at end, destructor should be called as appropriate
    // but this class has no destructor --- it doesn't store the size of its arrays
    // hence, manually deallocate here
    for (char component : {'x', 'y', 'z'}) {
      for (int k = 0; k < n_layers; k++) {
        for (int j = 0; j < n_cols; j++) { free(xyzt3d[component][k][j]); }
        free(xyzt3d[component][k]);
      }
      free(xyzt3d[component]);
    }
  }
}

TEST_CASE("XYZTensor") { XYZTensor3DTest().run_all_class_tests(); }
