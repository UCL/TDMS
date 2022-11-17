/**
 * @file test_FieldSample.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the FieldSample class and its subclasses
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

using namespace std;

const double TOLERANCE = 1e-16;

TEST_CASE("FieldSample") {

  mxArray *matlab_input;
  int dimensions[2] = {1, 1};
  const char *fieldnames[] = {"i", "j", "k", "n"};

  // attempting to construct with an empty array should exit with default values
  SECTION("Empty array input") {
    dimensions[0] = 0;
    matlab_input = mxCreateNumericArray(2, (const mwSize *) dimensions, mxUINT8_CLASS, mxREAL);
    // attempt assignment (2nd and 3rd args don't matter)
    FieldSample empty_fs(matlab_input);
    // should still be unassigned vectors
    bool no_vectors_set = (!empty_fs.i.has_elements()) && (!empty_fs.j.has_elements()) &&
                          (!empty_fs.k.has_elements()) && (!empty_fs.n.has_elements());
    REQUIRE(no_vectors_set);
  }

  // construction fails if provided a non-empty, non-struct array
  SECTION("Non-struct input") {
    dimensions[0] = 2;
    dimensions[1] = 3;
    matlab_input = mxCreateNumericArray(2, (const mwSize *) dimensions, mxUINT8_CLASS, mxREAL);
    CHECK_THROWS_AS(FieldSample(matlab_input), runtime_error);
  }

  // construction fails if provided a struct array that doesn't have 4 fields
  SECTION("Struct with too many fields") {
    const char *too_many_names[5] = {"field1", "field2", "field3", "field4", "field5"};
    matlab_input = mxCreateStructArray(2, (const mwSize *) dimensions, 5, too_many_names);
    CHECK_THROWS_AS(FieldSample(matlab_input), runtime_error);
  }
  SECTION("Struct with too few fields") {
    const char *too_few_names[1] = {"field1"};
    matlab_input = mxCreateStructArray(2, (const mwSize *) dimensions, 1, too_few_names);
    CHECK_THROWS_AS(FieldSample(matlab_input), runtime_error);
  }

  // expecting a struct with fields i, j, k, and n
  SECTION("Expected input (populated)") {
    const int number_of_vector_elements = 8;
    const int vector_dimensions[2] = {1, number_of_vector_elements};
    // create the struct
    matlab_input = mxCreateStructArray(2, (const mwSize *) dimensions, 4, fieldnames);
    // create vectors to place in the fields of the struct
    mxArray *i_vector =
            mxCreateNumericArray(2, (const mwSize *) vector_dimensions, mxINT32_CLASS, mxREAL);
    mxArray *j_vector =
            mxCreateNumericArray(2, (const mwSize *) vector_dimensions, mxINT32_CLASS, mxREAL);
    mxArray *k_vector =
            mxCreateNumericArray(2, (const mwSize *) vector_dimensions, mxINT32_CLASS, mxREAL);
    mxArray *n_vector =
            mxCreateNumericArray(2, (const mwSize *) vector_dimensions, mxDOUBLE_CLASS, mxREAL);
    // populate the vectors with some information
    mxDouble *place_i_data = mxGetPr(i_vector);//< 0,1,2,...
    mxDouble *place_j_data = mxGetPr(j_vector);//< number_of_vector_elements-1, -2, -3...
    mxDouble *place_k_data = mxGetPr(k_vector);//< 0,1,2,...,N/2,0,1,2,...
    mxDouble *place_n_data = mxGetPr(n_vector);//< 1, 1/2, 1/3, ...
    for (int i = 0; i < 4; i++) {
      place_i_data[i] = i;
      place_j_data[i] = number_of_vector_elements - (i + 1);
      place_k_data[i] = i % (number_of_vector_elements / 2);
      place_n_data[i] = 1. / (double) (i + 1);
    }
    // set the vectors to be the field values
    mxSetField(matlab_input, 0, fieldnames[0], i_vector);
    mxSetField(matlab_input, 0, fieldnames[1], j_vector);
    mxSetField(matlab_input, 0, fieldnames[2], k_vector);
    mxSetField(matlab_input, 0, fieldnames[3], n_vector);
    // we should be able to construct things now
    REQUIRE_NOTHROW(FieldSample(matlab_input));
    FieldSample fs(matlab_input);
    // check that the tensor attribute is the correct size
    const mwSize *dimensions = mxGetDimensions(fs.mx);
    for (int i = 0; i < 4; i++) { CHECK(dimensions[i] == number_of_vector_elements); }
    // check that we did indeed copy the data across
    for (int i = 0; i < 4; i++) {
      CHECK(abs(fs.i[i] - i) < TOLERANCE);
      CHECK(abs(fs.j[i] - (number_of_vector_elements - (i + 1))) < TOLERANCE);
      CHECK(abs(fs.k[i] - (i % (number_of_vector_elements / 2))) < TOLERANCE);
      CHECK(abs(fs.n[i] - 1. / (double) (i + 1)) < TOLERANCE);
    }
  }

  // alternatively, we can pass in empty vectors to the constructor and it should create an empty 4D tensor
  SECTION("Expected inputs (empty vectors)") {
    const int empty_vector_dimensions[2] = {1, 0};
    // create the struct
    matlab_input = mxCreateStructArray(2, (const mwSize *) dimensions, 4, fieldnames);
    // create vectors to place in the fields of the struct
    mxArray *empty_i_vector =
            mxCreateNumericArray(2, (const mwSize *) empty_vector_dimensions, mxINT32_CLASS, mxREAL);
    mxArray *empty_j_vector =
            mxCreateNumericArray(2, (const mwSize *) empty_vector_dimensions, mxINT32_CLASS, mxREAL);
    mxArray *empty_k_vector =
            mxCreateNumericArray(2, (const mwSize *) empty_vector_dimensions, mxINT32_CLASS, mxREAL);
    mxArray *empty_n_vector =
            mxCreateNumericArray(2, (const mwSize *) empty_vector_dimensions, mxDOUBLE_CLASS, mxREAL);
    mxSetField(matlab_input, 0, fieldnames[0], empty_i_vector);
    mxSetField(matlab_input, 0, fieldnames[1], empty_j_vector);
    mxSetField(matlab_input, 0, fieldnames[2], empty_k_vector);
    mxSetField(matlab_input, 0, fieldnames[3], empty_n_vector);
    // we should be able to construct things now
    REQUIRE_NOTHROW(FieldSample(matlab_input));
    FieldSample fs_from_empty_vectors(matlab_input);
    // check that all vectors are empty flags as being true
    CHECK(!fs_from_empty_vectors.all_vectors_are_non_empty());
    // check that the tensor attribute is the correct size
    const mwSize *empty_tensor_dimensions = mxGetDimensions(fs_from_empty_vectors.mx);
    for (int i = 0; i < 4; i++) { CHECK(empty_tensor_dimensions[i] == 0); }
  }

  // cleanup
  mxDestroyArray(matlab_input);
}
