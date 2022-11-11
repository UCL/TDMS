/**
 * @file test_Matrix.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests for the Matrix class and its subclasses (Vertices, GratingStructure)
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

TEST_CASE("GratingStructure: allocation and deallocation") {
  SPDLOG_INFO("== Testing GratingStructure allocation/deallocation");

  // non-empty input must be a pointer to a 2D matlab array (of ints, although non-interleaved API does not grant us the luxury of enforcing this),
  // dimensions must be 2 by I_tot+1
  const int I_tot = 4;

  // empty pointer input doesn't throw an error, but also doesn't initialise anything meaningful
  mxArray *empty_array = mxCreateNumericMatrix(0, 4, mxINT32_CLASS, mxREAL);
  GratingStructure *created_with_empty_ptr;
  REQUIRE_NOTHROW(created_with_empty_ptr = new GratingStructure(empty_array, 4));
  REQUIRE(!created_with_empty_ptr->has_elements());
  // memory cleanup
  delete created_with_empty_ptr;
  mxDestroyArray(empty_array);

  // try sending in a 3D array instead
  int dims_3d[3] = {2, I_tot+1, 3};
  mxArray *array_3by3 = mxCreateNumericArray(3, (mwSize *) dims_3d, mxINT32_CLASS, mxREAL);
  REQUIRE_THROWS_AS(GratingStructure(array_3by3, I_tot), std::runtime_error);
  mxDestroyArray(array_3by3);

  // try sending in 2D arrays with the wrong dimensions
  // wrong dimension 1
  int dims_2d[2] = {2, I_tot};
  mxArray *array_2d_wrong_dim1 = mxCreateNumericArray(2, (mwSize *) dims_2d, mxINT32_CLASS, mxREAL);
  REQUIRE_THROWS_AS(GratingStructure(array_2d_wrong_dim1, I_tot), std::runtime_error);
  // wrong dimension 0
  dims_2d[0] = 3; dims_2d[1] = I_tot + 1;
  mxArray *array_2d_wrong_dim0 = mxCreateNumericArray(2, (mwSize *) dims_2d, mxINT32_CLASS, mxREAL);
  REQUIRE_THROWS_AS(GratingStructure(array_2d_wrong_dim0, I_tot), std::runtime_error);
  // memory cleanup
  mxDestroyArray(array_2d_wrong_dim0);
  mxDestroyArray(array_2d_wrong_dim1);

  // now try sending in something that's actually useful
  dims_2d[0] = 2; dims_2d[1] = I_tot + 1; // reassign, but better to be safe
  mxArray *useful_array = mxCreateNumericMatrix(dims_2d[0], dims_2d[1], mxINT32_CLASS, mxREAL);
  GratingStructure *gs;
  REQUIRE_NOTHROW(gs = new GratingStructure(useful_array, I_tot));
  REQUIRE(gs->has_elements());

  // memory cleanup
  /* Disassociate the GratingStructure from the MATLAB array that it was cast to.
  In iterator.cpp, this ensures that gs.matrix is freed and set to nullptr, but doesn't free the memory that matrix[j] was pointing to, since this is returned through the MEX function.
  Since we are unit testing here, we have dynamically created this memory for the MATLAB array, and no longer require it, thus we also need to cleanup our MATLAB array manually.
  */
  delete gs; // disassociate from MATLAB array and free malloced space
  mxDestroyArray(useful_array); // clear MATLAB array and free mxMalloc'd space
}
