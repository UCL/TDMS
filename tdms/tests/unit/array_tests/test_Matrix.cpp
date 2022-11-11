/**
 * @file test_Matrix.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests for the Matrix class and its subclasses (Vertices, GratingStructure, Pupil, EHVec)
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

TEST_CASE("Matrix") {
  SPDLOG_INFO("== Matrix");

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

TEST_CASE("Vertices") {
  SPDLOG_INFO("== Vertices");

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

TEST_CASE("GratingStructure") {
  SPDLOG_INFO("== GratingStructure");

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

TEST_CASE("Pupil") {
  SPDLOG_INFO("== Pupil");

  // only default constructor exists, which doesn't even assign memory
  Pupil p;
  REQUIRE(!p.has_elements());

  // we'll use these as the target dimensions
  const int n_rows = 4, n_cols = 8;

  // passing in an empty array to initialise() doesn't error, but also doesn't assign
  // additionally, the rows and columns arguments aren't even used, so can be garbage
  mxArray *empty_array = mxCreateNumericMatrix(0, n_cols, mxDOUBLE_CLASS, mxREAL);
  p.initialise(empty_array, 1, 1);
  REQUIRE(!p.has_elements()); // shouldn't have assigned any memory or pointers
  mxDestroyArray(empty_array);

  // wrong dimensions or wrong number of dimensions will cause an error, and also not assign
  // wrong number of dimensions
  const int dims_3d[3] = {n_rows, n_cols, 2};
  mxArray *array_3d = mxCreateNumericArray(3, (mwSize *) dims_3d, mxDOUBLE_CLASS, mxREAL);
  REQUIRE_THROWS_AS(p.initialise(array_3d, n_rows, n_cols), std::runtime_error);
  CHECK(!p.has_elements());
  mxDestroyArray(array_3d);
  // wrong number of rows/cols
  mxArray *array_wrong_dims = mxCreateNumericMatrix(2 * n_rows, n_cols + 1, mxDOUBLE_CLASS, mxREAL);
  REQUIRE_THROWS_AS(p.initialise(array_wrong_dims, n_rows, n_cols), std::runtime_error);
  CHECK(!p.has_elements());
  mxDestroyArray(array_wrong_dims);

  // correct size should successfully assign memory
  mxArray *array_2d = mxCreateNumericMatrix(n_rows, n_cols, mxDOUBLE_CLASS, mxREAL);
  REQUIRE_NOTHROW(p.initialise(array_2d, n_rows, n_cols));
  REQUIRE(p.has_elements());
  mxDestroyArray(array_2d); // need to manually cleanup our MATLAB array, since we are not preserving the data via MEX function
}

TEST_CASE("EHVec") {
  SPDLOG_INFO("== EHVec");

  // because we're storing fftw_complex variables, this class has a custom destructor but nothing else
  // as such, we should just be able to initialise it using allocate as per
  const int n_rows = 4, n_cols = 8;
  const int REAL = 0, IMAG = 1;
  EHVec eh;
  CHECK(!eh.has_elements());

  // allocate memory
  eh.allocate(n_rows, n_cols);
  CHECK(eh.has_elements());

  // check that we an assign fftw_complexes to the elements
  eh[0][0][REAL] = 1.; eh[0][0][IMAG] = 0.;
  eh[0][1][REAL] = 0.; eh[0][1][IMAG] = 1.;
  fftw_complex fftw_unit {1.,0.};
  fftw_complex fftw_imag_unit {0.,1.};

  CHECK(eh[0][0][REAL] == fftw_unit[REAL]);
  CHECK(eh[0][0][IMAG] == fftw_unit[IMAG]);
  CHECK(eh[0][1][REAL] == fftw_imag_unit[REAL]);
  CHECK(eh[0][1][IMAG] == fftw_imag_unit[IMAG]);
}
