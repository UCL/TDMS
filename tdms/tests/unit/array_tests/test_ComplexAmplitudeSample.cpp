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

TEST_CASE("ComplexAmplitudeSample: allocation and deallocation") {
  SPDLOG_INFO("== Testing ComplexAmplitudeSample class");

  // Constructor should just exit if recieving an empty struct
  // POTENTIALLY DANGEROUS BEHAVIOUR AS WE RE THEN ABLE TO ACCESS OTHER METHODS!!!
  const char *empty_fields[] = {};
  const int empty_dims[2] = {1, 1};
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
  const int n_field_elements = 5;
  // build our struct
  const int dims[2] = {1, 2};
  mxArray *cas_struct_pointer = mxCreateStructArray(2, (const mwSize *) dims, n_fields, fieldnames);
  // each entry is another structure array, setup for vertices and components
  // copy code from Vertices and FieldComponentsVector once those tests have been added
  mxDestroyArray(cas_struct_pointer);
}
