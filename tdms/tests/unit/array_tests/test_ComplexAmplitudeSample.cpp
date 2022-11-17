/**
 * @file test_ComplexAmplitudeSample.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the ComplexAmplitudeSample class and its subclasses
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"
#include "globals.h"

using namespace std;

TEST_CASE("ComplexAmplitudeSample") {

  mxArray *matlab_input;
  int dimensions[2] = {1, 1};
  const int n_fields = 2;
  const char *fieldnames[n_fields] = {"vertices", "components"};

  // Constructor should just exit if recieving an empty struct
  SECTION("Empty input") {
    dimensions[0] = 0;
    matlab_input = mxCreateStructArray(2, (const mwSize *) dimensions, n_fields, fieldnames);
    CHECK_NOTHROW(ComplexAmplitudeSample(matlab_input));
    ComplexAmplitudeSample empty_test(matlab_input);
    // n_vertices should return 0, since Vector has not been set so should default initialise to 0-vertex matrix
    // similarly, vector & components should be flagged as having no elements
    CHECK(empty_test.n_vertices() == 0);
    CHECK(!empty_test.vertices.has_elements());
    CHECK(!empty_test.components.has_elements());
  }

  // For successful construction, we need to build a MATLAB struct with 2 fields
  // these are the fieldnames that are expected
  SECTION("Expected input") {
    matlab_input = mxCreateStructArray(2, (const mwSize *) dimensions, n_fields, fieldnames);
    // each entry is a numeric array, setup for vertices and components
    // copy code from Vertices and FieldComponentsVector once those tests have been added
    const int n_elements = 8;
    // array for "vertices" field
    mxArray *vertices_array = mxCreateNumericMatrix(n_elements, 3, mxINT32_CLASS, mxREAL);
    // array for "components" field
    mxArray *components_vector = mxCreateNumericMatrix(1, n_elements, mxINT16_CLASS, mxREAL);
    // add fields to struct that will be passed to CAS constructor function
    mxSetField(matlab_input, 0, fieldnames[0], vertices_array);
    mxSetField(matlab_input, 0, fieldnames[1], components_vector);
    // create an instance using this struct
    REQUIRE_NOTHROW(ComplexAmplitudeSample(matlab_input));
  }

  // tear down
  mxDestroyArray(matlab_input);
}
