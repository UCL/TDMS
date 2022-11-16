/**
 * @file test_Matrix.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests for the Matrix class and its subclasses
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

TEST_CASE("Matrix: allocation and deallocation") {

  // mock-up dimensions
  int n_rows = 4, n_cols = 8;

  // create a Matrix via the default constructor
  SECTION("Default constructor") {
    Matrix<int> M_int;
    // although created, the matrix should not have any elements, and should be flagged as such
    REQUIRE(!M_int.has_elements());
    // we need to use allocate to provide space for elements
    M_int.allocate(n_rows, n_cols);
    REQUIRE(M_int.has_elements());

    // should be able to assign to these values without seg faults now
    for (int i = 0; i < n_rows; i++) {
      for (int j = 0; j < n_cols; j++) {
        CHECK_NOTHROW(M_int[i][j] = 0);
      }
    }
  }

  // create a Matrix using the overloaded constructor, and fill it with 1s
  SECTION("Overloaded constructor") {
    Matrix<double> M_double(n_rows, n_cols);
    // because we used the overloaded constructor, the matrix should have elements
    REQUIRE(M_double.has_elements());

    // should be able to assign to these values without seg faults now
    for (int i = 0; i < n_rows; i++) {
      for (int j = 0; j < n_cols; j++) {
        CHECK_NOTHROW(M_double[i][j] = 1.);
      }
    }
  }
}
