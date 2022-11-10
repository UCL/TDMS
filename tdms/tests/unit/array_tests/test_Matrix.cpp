/**
 * @file test_Matrix.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests for the Matrix class and its subclasses
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

TEST_CASE("Matrix: allocation and deallocation") {
  SPDLOG_INFO("== Testing Matrix allocation/deallocation");

  // mock-up dimensions
  int n_rows = 4, n_cols = 8;

  // create a Matrix via the default constructor
  Matrix<int> M_int;
  // although created, the matrix should not have any elements, and should be flagged as such
  REQUIRE(!M_int.has_elements());
  // we need to use allocate to provide space for elements
  M_int.allocate(n_rows, n_cols);
  REQUIRE(M_int.has_elements());

  // create a Matrix using the overloaded constructor, and fill it with 1s
  Matrix<double> M_double(n_rows, n_cols);
  // because we used the overloaded constructor, the matrix should have elements
  REQUIRE(M_double.has_elements());

  // should be able to assign to these values without seg faults now
  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < n_cols; j++) {
      CHECK_NOTHROW(M_int[i][j] = 0);
      CHECK_NOTHROW(M_double[i][j] = 1.);
    }
  }
}

TEST_CASE("Vertices: allocation and deallocation") {
  SPDLOG_INFO("== Testing Vertices allocation/deallocation");

  // initialise the struct, it needs the fieldname "vertices"
  const int n_fields = 1;
  const int n_vertex_elements = 8;
  const char *fieldnames[n_fields] = {"vertices"};
  mxArray *struct_pointer = mxCreateStructMatrix(1, 1, n_fields, fieldnames);
  mxArray *vertices_array = mxCreateNumericMatrix(n_vertex_elements, 3, mxINT32_CLASS, mxREAL);
  mxSetField(struct_pointer, 0, fieldnames[0], vertices_array);
  // create object
  Vertices v;
  // initialise
  v.initialise(struct_pointer);
  // we should have n_vertex_elements number of vertices stored
  REQUIRE(v.n_vertices() == n_vertex_elements);
  // what's more, we should have decremented all the elements of v from 0 to -1
  for (int i = 0; i < n_vertex_elements; i++) {
   for (int j = 0; j < 3; j++) {
     CHECK(v[j][i] == -1);
   }
  }

  // destroy MATLAB object
  mxDestroyArray(struct_pointer);
}
