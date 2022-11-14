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

  // Constructor should just exit if recieving an empty struct
  const char *empty_fields[] = {};
  const int empty_dims[2] = {0, 1};
  mxArray *empty_struct_pointer =
          mxCreateStructArray(2, (const mwSize *) empty_dims, 0, empty_fields);
  CHECK_NOTHROW(ComplexAmplitudeSample(empty_struct_pointer));
  ComplexAmplitudeSample empty_test(empty_struct_pointer);
  // n_vertices should return 0, since Vector has not been set so should default initialise to 0-vertex matrix
  // similarly, vector & components should be flagged as having no elements
  CHECK(empty_test.n_vertices()==0);
  CHECK(!empty_test.vertices.has_elements());
  CHECK(!empty_test.components.has_elements());
  mxDestroyArray(empty_struct_pointer);

  // For successful construction, we need to build a MATLAB struct with 2 fields
  // these are the fieldnames that are expected
  const int n_fields = 2;
  const char *fieldnames[n_fields] = {"vertices", "components"};
  // build our struct
  const int dims[2] = {1, 2};
  mxArray *cas_struct_pointer = mxCreateStructArray(2, (const mwSize *) dims, n_fields, fieldnames);
  // each entry is a numeric array, setup for vertices and components
  // copy code from Vertices and FieldComponentsVector once those tests have been added
  const int n_elements = 8;
  // array for "vertices" field
  mxArray *vertices_array = mxCreateNumericMatrix(n_elements, 3, mxINT32_CLASS, mxREAL);
  // array for "components" field
  mxArray *components_vector = mxCreateNumericMatrix(1, n_elements, mxINT16_CLASS, mxREAL);
  // add fields to struct that will be passed to CAS constructor function
  mxSetField(cas_struct_pointer, 0, fieldnames[0], vertices_array);
  mxSetField(cas_struct_pointer, 0, fieldnames[1], components_vector);
  // create an instance using this struct
  REQUIRE_NOTHROW(ComplexAmplitudeSample(cas_struct_pointer));
  // do whatever checks we need

  // destroy what we created
  mxDestroyArray(cas_struct_pointer);
}
