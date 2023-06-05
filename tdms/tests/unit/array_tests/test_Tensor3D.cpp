/**
 * @file test_Tensor3D.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests the functionality of the Tensor3D class, which is the building
 * block for several further field classes
 *
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "array_test_class.h"
#include "arrays.h"
#include "unit_test_utils.h"

using tdms_tests::TOLERANCE;

void Tensor3DTest::test_correct_construction() {
  Tensor3D<int> t3d;
  SECTION("Default constructor") {
    // default constructor should assign all dimensions to 0, and the tensor
    // should have no elements as a result
    REQUIRE(!(t3d.has_elements()));
    // let us assign some free memory to this tensor via allocate()
    t3d.allocate(n_layers, n_cols, n_rows);
    // we should now "have elements", even though they are unassigned
    REQUIRE(t3d.has_elements());
  }
  SECTION("Overloaded constructor") {
    // testing initialise() as the overloaded constructor is a wrapper for this
    int ***p = nullptr;
    // assemble p such that the [i][j][k]th entry is the strided vector index,
    // i*n_layers*n_cols + j*n_layers + k
    SECTION("Reverse buffer") {
      // Test when we read in a 3D buffer of shape [n_layers][n_cols][n_rows]
      p = (int ***) malloc(n_layers * sizeof(int **));
      for (int k = 0; k < n_layers; k++) {
        p[k] = (int **) malloc(n_cols * sizeof(int *));
        for (int j = 0; j < n_cols; j++) {
          p[k][j] = (int *) malloc(n_rows * sizeof(int));
          for (int i = 0; i < n_rows; i++) {
            p[k][j][i] = i * n_layers * n_cols + j * n_layers + k;
          }
        }
      }
      t3d.initialise(p, n_layers, n_cols, n_rows, true);
      // Malloc'd values are copied, so we can tear down now
      free(p);
    }
    SECTION("Ordered buffer") {
      // Test when we read in a 3D buffer of shape [n_rows][n_cols][n_layers]
      p = (int ***) malloc(n_rows * sizeof(int **));
      for (int i = 0; i < n_rows; i++) {
        p[i] = (int **) malloc(n_cols * sizeof(int *));
        for (int j = 0; j < n_cols; j++) {
          p[i][j] = (int *) malloc(n_rows * sizeof(int));
          for (int k = 0; k < n_layers; k++) {
            p[i][j][k] = i * n_layers * n_cols + j * n_layers + k;
          }
        }
      }
      t3d.initialise(p, n_layers, n_cols, n_rows, false);
      // Malloc'd values are copied, so we can tear down now
      free(p);
    }
    // this tensor should be flagged as "having elements", since we provided a
    // pointer in the constructor
    bool elements_are_correct = t3d.has_elements();
    for (int i = 0; i < n_rows; i++) {
      for (int j = 0; j < n_cols; j++) {
        for (int k = 0; k < n_layers; k++) {
          elements_are_correct =
                  elements_are_correct &&
                  (t3d(i, j, k) == i * n_layers * n_cols + j * n_layers + k);
        }
      }
    }
    REQUIRE(t3d.has_elements());
  }
}

void Tensor3DTest::test_other_methods() {
  Tensor3D<double> t3d(n_layers, n_cols, n_rows);
  t3d.zero();
  SECTION("allocate() and zero()") {
    // we should be able to flag this tensor has elements, so the bool should be
    // set to true
    bool allocated_and_zero = t3d.has_elements();
    for (int k = 0; k < n_layers; k++) {
      for (int j = 0; j < n_cols; j++) {
        for (int i = 0; i < n_rows; i++) {
          allocated_and_zero =
                  allocated_and_zero && (abs(t3d(i, j, k)) < TOLERANCE);
        }
      }
    }
    REQUIRE(allocated_and_zero);
  }
  SECTION("frobenius()") {
    // frobenius norm should be zero after allocation and zero-ing
    REQUIRE(abs(t3d.frobenius()) < TOLERANCE);
    // assign some values to this tensor. We'll go with =0 if i+j+k is even, and
    // =1 if odd this gives us 4*8*16/2 = 4^4 = 256 non-zero entries, which are
    // 1, so the analytic norm is 16.
    for (int k = 0; k < n_layers; k++) {
      for (int j = 0; j < n_cols; j++) {
        for (int i = 0; i < n_rows; i++) { t3d(i, j, k) = (i + j + k) % 2; }
      }
    }
    // the analytic frobenuis norm of the tensor values
    double target_fro = 16.;
    // check the frobenius norms align
    REQUIRE(abs(t3d.frobenius() - target_fro) < TOLERANCE);
  }
}

TEST_CASE("Tensor3D: zero, frobenius") { Tensor3DTest().run_all_class_tests(); }
