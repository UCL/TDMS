/**
 * @file test_Matrix.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests for the Matrix class and its subclasses (Vertices,
 * GratingStructure, Pupil)
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "array_test_class.h"
#include "arrays/tdms_matrix.h"

void MatrixTest::test_correct_construction() {
  // create a Matrix via the default constructor
  SECTION("Default constructor") {
    Matrix<int> M_int;
    // although created, the matrix should not have any elements, and should be
    // flagged as such
    REQUIRE(!M_int.has_elements());
  }
  // create a Matrix using the overloaded constructor, and fill it with 1s
  SECTION("Overloaded constructor") {
    Matrix<double> M_double(n_rows, n_cols);
    // because we used the overloaded constructor, the matrix should have
    // elements
    REQUIRE(M_double.has_elements());

    // should be able to assign to these values without seg faults now
    for (int i = 0; i < n_rows; i++) {
      for (int j = 0; j < n_cols; j++) { CHECK_NOTHROW(M_double[i][j] = 1.); }
    }
  }
}

void MatrixTest::test_other_methods() {
  SECTION("allocate()") {
    Matrix<int> M_int;
    // we need to use allocate to provide space for elements
    M_int.allocate(n_rows, n_cols);
    REQUIRE(M_int.has_elements());

    // should be able to assign to these values without seg faults now
    for (int i = 0; i < n_rows; i++) {
      for (int j = 0; j < n_cols; j++) { CHECK_NOTHROW(M_int[i][j] = 0); }
    }
  }
}

void VerticesTest::test_correct_construction() {
  // initialise the struct, it needs the fieldname "vertices"
  create_1by1_struct(n_fields, fieldnames);
  mxArray *vertices_array =
          mxCreateNumericMatrix(n_numeric_elements, 3, mxINT32_CLASS, mxREAL);
  mxSetField(matlab_input, 0, fieldnames[0], vertices_array);
  // create object
  Vertices v;
  // initialise
  v.initialise_from_matlab(matlab_input);
  // we should have n_vertex_elements number of vertices stored
  REQUIRE(v.n_vertices() == n_numeric_elements);
  // what's more, we should have decremented all the elements of v from 0 to -1
  for (int i = 0; i < n_numeric_elements; i++) {
    for (int j = 0; j < 3; j++) { CHECK(v(i, j) == -1); }
  }
}

void GratingStructureTest::test_empty_construction() {
  dimensions_2d[0] = 0;
  dimensions_2d[1] = I_tot;
  create_numeric_array(2, dimensions_2d, mxINT32_CLASS);
  // note: "new" used here since we need to delete gs to then safely delete
  // matlab_input AFTERWARDS we cannot delete matlab_array before gs
  GratingStructure *gs;
  REQUIRE_NOTHROW(gs = new GratingStructure(matlab_input, I_tot));
  REQUIRE(!gs->has_elements());
  // tear down
  delete gs;
}

void GratingStructureTest::test_wrong_input_dimensions() {
  SECTION("Wrong number of dimensions (3D)") {
    dimensions_3d[0] = 2;
    dimensions_3d[1] = I_tot;
    dimensions_3d[2] = 3;
    create_numeric_array(3, dimensions_3d, mxINT32_CLASS);
    REQUIRE_THROWS_AS(GratingStructure(matlab_input, I_tot),
                      std::runtime_error);
  }
  SECTION("Wrong dimension size") {
    SECTION("(axis 0)") {
      dimensions_2d[0] = 3;
      dimensions_2d[1] = I_tot + 1;
    }
    SECTION("(axis 1)") {
      dimensions_2d[0] = 2;
      dimensions_2d[1] = I_tot;
    }
    create_numeric_array(2, dimensions_2d, mxINT32_CLASS);
    REQUIRE_THROWS_AS(GratingStructure(matlab_input, I_tot),
                      std::runtime_error);
  }
}

void GratingStructureTest::test_correct_construction() {
  dimensions_2d[0] = 2;
  dimensions_2d[1] = I_tot + 1;
  create_numeric_array(2, dimensions_2d, mxINT32_CLASS, mxREAL);
  // note: "new" used here since we need to delete gs to then safely delete
  // matlab_input AFTERWARDS we cannot delete matlab_array before gs
  GratingStructure *gs;
  REQUIRE_NOTHROW(gs = new GratingStructure(matlab_input, I_tot));
  REQUIRE(gs->has_elements());
  // tear down
  delete gs;
}

void PupilTest::test_empty_construction() {
  Pupil p;
  dimensions_2d[0] = 0;
  dimensions_2d[1] = n_cols;
  create_numeric_array(2, dimensions_2d);
  // passing in an empty array to initialise() doesn't error, but also doesn't
  // assign additionally, the rows and columns arguments aren't even used, so
  // can be garbage
  p.initialise_from_matlab(matlab_input, 1, 1);
  REQUIRE(!p.has_elements());// shouldn't have assigned any memory or pointers
}

void PupilTest::test_wrong_input_dimensions() {
  Pupil p;
  SECTION("Wrong number of dimensions (3D)") {
    dimensions_3d[0] = n_rows;
    dimensions_3d[1] = n_cols;
    dimensions_3d[2] = 2;
    create_numeric_array(3, dimensions_3d);
    REQUIRE_THROWS_AS(p.initialise_from_matlab(matlab_input, n_rows, n_cols),
                      std::runtime_error);
    REQUIRE(!p.has_elements());
  }
  SECTION("Wrong dimension size") {
    dimensions_2d[0] = 2 * n_rows;
    dimensions_2d[1] = n_cols + 1;
    create_numeric_array(2, dimensions_2d);
    REQUIRE_THROWS_AS(p.initialise_from_matlab(matlab_input, n_rows, n_cols),
                      std::runtime_error);
    REQUIRE(!p.has_elements());
  }
}

void PupilTest::test_correct_construction() {
  Pupil p;
  SECTION("Default constructor") { REQUIRE(!p.has_elements()); }
  SECTION("Overloaded constructor") {
    dimensions_2d[0] = n_rows;
    dimensions_2d[1] = n_cols;
    create_numeric_array(2, dimensions_2d);
    REQUIRE_NOTHROW(p.initialise_from_matlab(matlab_input, n_rows, n_cols));
    REQUIRE(p.has_elements());
  }
}

TEST_CASE("Matrix") { MatrixTest().run_all_class_tests(); }

TEST_CASE("Vertices") { VerticesTest().run_all_class_tests(); }

TEST_CASE("GratingStructure") { GratingStructureTest().run_all_class_tests(); }

TEST_CASE("Pupil") { PupilTest().run_all_class_tests(); }
