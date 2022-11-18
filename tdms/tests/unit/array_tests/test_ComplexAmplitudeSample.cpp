/**
 * @file test_ComplexAmplitudeSample.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the ComplexAmplitudeSample class and its subclasses
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "array_test_class.h"
#include "arrays.h"
#include "globals.h"

using namespace std;

bool ComplexAmplitudeSampleTest::test_empty_construction() {
  create_empty_struct(2, fieldnames);
  CHECK_NOTHROW(ComplexAmplitudeSample(matlab_input));
  ComplexAmplitudeSample empty_test(matlab_input);
  // n_vertices should return 0, since Vector has not been set so should default initialise to 0-vertex matrix
  // similarly, vector & components should be flagged as having no elements
  CHECK(empty_test.n_vertices() == 0);
  CHECK(!empty_test.vertices.has_elements());
  CHECK(!empty_test.components.has_elements());
  return true;
}

bool ComplexAmplitudeSampleTest::test_correct_construction() {
  create_struct_array(2, dimensions_2d, 2, fieldnames);
  // each entry is a numeric array, setup for vertices and components
  // array for "vertices" field
  mxArray *vertices_array = mxCreateNumericMatrix(n_numeric_elements, 3, mxINT32_CLASS, mxREAL);
  // array for "components" field
  mxArray *components_vector = mxCreateNumericMatrix(1, n_numeric_elements, mxINT16_CLASS, mxREAL);
  // add fields to struct that will be passed to CAS constructor function
  mxSetField(matlab_input, 0, fieldnames[0], vertices_array);
  mxSetField(matlab_input, 0, fieldnames[1], components_vector);
  // create an instance using this struct
  REQUIRE_NOTHROW(ComplexAmplitudeSample(matlab_input));
  return true;
}

TEST_CASE("ComplexAmplitudeSample") { ComplexAmplitudeSampleTest().run_all_class_tests(); }
