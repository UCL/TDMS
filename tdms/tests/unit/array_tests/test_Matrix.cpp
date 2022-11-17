/**
 * @file test_Matrix.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests for the Matrix class and its subclasses (Vertices, GratingStructure, Pupil, EHVec)
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

TEST_CASE("Matrix") {

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

TEST_CASE("Vertices") {

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

  // tear down - destroy MATLAB object
  mxDestroyArray(struct_pointer);
}

TEST_CASE("GratingStructure") {

  // non-empty input must be a pointer to a 2D matlab array (of ints, although non-interleaved API does not grant us the luxury of enforcing this),
  // dimensions must be 2 by I_tot+1
  const int I_tot = 4;
  mxArray *matlab_input;

  // empty pointer input doesn't throw an error, but also doesn't initialise anything meaningful
  SECTION("Empty input") {
    matlab_input = mxCreateNumericMatrix(0, 4, mxINT32_CLASS, mxREAL);
    GratingStructure *gs;
    REQUIRE_NOTHROW(gs = new GratingStructure(matlab_input, 4));
    REQUIRE(!gs->has_elements());
    // memory cleanup
    delete gs;
  }

  // try sending in a 3D array instead
  SECTION("Wrong number of dimensions (3D)") {
    int dimensions_3d[3] = {2, I_tot + 1, 3};
    matlab_input = mxCreateNumericArray(3, (mwSize *) dimensions_3d, mxINT32_CLASS, mxREAL);
    REQUIRE_THROWS_AS(GratingStructure(matlab_input, I_tot), std::runtime_error);
  }
  SECTION("Wrong dimension size") {
    int dimensions_2d[2];
    SECTION("(axis 0)") {
      dimensions_2d[0] = 3;
      dimensions_2d[1] = I_tot + 1;
    }
    SECTION("(axis 1)") {
      dimensions_2d[0] = 2;
      dimensions_2d[1] = I_tot;
    }
    matlab_input = mxCreateNumericArray(2, (mwSize *) dimensions_2d, mxINT32_CLASS, mxREAL);
    REQUIRE_THROWS_AS(GratingStructure(matlab_input, I_tot), std::runtime_error);
  }

  // now try sending in something that's actually useful
  SECTION("Expected input") {
    int dimensions_2d[2] = {2, I_tot + 1};
    matlab_input = mxCreateNumericMatrix(dimensions_2d[0], dimensions_2d[1], mxINT32_CLASS, mxREAL);
    GratingStructure *gs;
    REQUIRE_NOTHROW(gs = new GratingStructure(matlab_input, I_tot));
    REQUIRE(gs->has_elements());
    // memory cleanup
    delete gs;
  }

  mxDestroyArray(matlab_input);
}

TEST_CASE("Pupil") {

  // only default constructor exists, which doesn't even assign memory
  Pupil p;
  REQUIRE(!p.has_elements());
  mxArray *matlab_input;
  // we'll use these as the target dimensions
  const int n_rows = 4, n_cols = 8;

  // passing in an empty array to initialise() doesn't error, but also doesn't assign
  // additionally, the rows and columns arguments aren't even used, so can be garbage
  SECTION("Empty input") {
    matlab_input = mxCreateNumericMatrix(0, n_cols, mxDOUBLE_CLASS, mxREAL);
    p.initialise(matlab_input, 1, 1);
    REQUIRE(!p.has_elements());// shouldn't have assigned any memory or pointers
  }

  // wrong dimensions or wrong number of dimensions will cause an error, and also not assign
  SECTION("Wrong number of dimensions (3D)") {
    const int dimensions_3d[3] = {n_rows, n_cols, 2};
    matlab_input = mxCreateNumericArray(3, (mwSize *) dimensions_3d, mxDOUBLE_CLASS, mxREAL);
    REQUIRE_THROWS_AS(p.initialise(matlab_input, n_rows, n_cols), std::runtime_error);
    CHECK(!p.has_elements());
  }
  SECTION("Wrong dimension size") {
    matlab_input = mxCreateNumericMatrix(2 * n_rows, n_cols + 1, mxDOUBLE_CLASS, mxREAL);
    REQUIRE_THROWS_AS(p.initialise(matlab_input, n_rows, n_cols), std::runtime_error);
    CHECK(!p.has_elements());
  }

  // correct size should successfully assign memory
  SECTION("Expected input") {
    matlab_input = mxCreateNumericMatrix(n_rows, n_cols, mxDOUBLE_CLASS, mxREAL);
    REQUIRE_NOTHROW(p.initialise(matlab_input, n_rows, n_cols));
    REQUIRE(p.has_elements());
  }

  mxDestroyArray(matlab_input);
}

TEST_CASE("EHVec") {

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
