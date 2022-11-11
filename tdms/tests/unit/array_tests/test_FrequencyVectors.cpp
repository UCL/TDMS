/**
 * @file test_FrequencyVectors.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests for the FrequencyVectors class and its subclasses
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

using namespace std;

const double tol = 1e-16;

TEST_CASE("FrequencyVectors") {
  SPDLOG_INFO("== FrequencyVectors");
  // there is no custom constructor for this class
  FrequencyVectors fv;
  // but the members should start as nullptrs
  bool not_assigned = (!fv.x.has_elements() && !fv.y.has_elements());
  REQUIRE(not_assigned);
  // initialise() method should exit without assignment if we pass in a pointer to an empty array
  // create empty array
  const int empty_dims[2] = {0, 1};
  mxArray *empty_array =
          mxCreateNumericArray(2, (const mwSize *) empty_dims, mxUINT8_CLASS, mxREAL);
  // attempt assignment
  fv.initialise(empty_array);
  // should still be unassigned vectors
  not_assigned = (!fv.x.has_elements() && !fv.y.has_elements());
  REQUIRE(not_assigned);
  // cleanup the empty array mess we just made
  mxDestroyArray(empty_array);

  // assignment will throw error if we attempt to provide a non-empty, non-struct array
  int dims[2] = {2, 3};
  mxArray *not_a_struct = mxCreateNumericArray(2, (const mwSize *) dims, mxUINT8_CLASS, mxREAL);
  CHECK_THROWS_AS(fv.initialise(not_a_struct), runtime_error);
  // cleanup
  mxDestroyArray(not_a_struct);

  // assignment will throw error if we attempt to provide a struct array that doesn't have two fields
  const char *too_few_names[1] = {"field1"};
  const char *too_mny_names[3] = {"field1", "field2", "field3"};
  mxArray *struct_with_too_few_fields =
          mxCreateStructArray(2, (const mwSize *) dims, 1, too_few_names);
  mxArray *struct_with_too_mny_fields =
          mxCreateStructArray(2, (const mwSize *) dims, 3, too_mny_names);
  CHECK_THROWS_AS(fv.initialise(struct_with_too_few_fields), runtime_error);
  CHECK_THROWS_AS(fv.initialise(struct_with_too_mny_fields), runtime_error);
  // cleanup
  mxDestroyArray(struct_with_too_few_fields);
  mxDestroyArray(struct_with_too_mny_fields);

  // otherwise, we need to provide a struct, whose fields are vectors that can be assigned to
  // setup for our struct
  const char *fieldnames[] = {"fx_vec", "fy_vec"};
  dims[0] = 1;
  dims[1] = 1;
  const int n_field_array_elements = 10;
  const int field_array_dims[2] = {1, n_field_array_elements};
  // create the struct
  mxArray *example_struct = mxCreateStructArray(2, (const mwSize *) dims, 2, fieldnames);
  // create the data for the fields of our struct
  mxArray *field_array_ptrs[2];
  for (int i = 0; i < 2; i++) {
    field_array_ptrs[i] =
            mxCreateNumericArray(2, (const mwSize *) field_array_dims, mxDOUBLE_CLASS, mxREAL);
    mxDouble *where_to_place_data = mxGetPr(field_array_ptrs[i]);
    // 0th field, fx_vec[i], will be 1/(i+1)
    // 1st field, fy_vec[i], will be -1/(i+1)
    for (int j = 0; j < n_field_array_elements; j++) {
      where_to_place_data[j] = pow(-1., (double) i) / ((double) (j + 1));
    }
    mxSetField(example_struct, 0, fieldnames[i], field_array_ptrs[i]);
  }
  // attempt to create a vector from this struct
  REQUIRE_NOTHROW(fv.initialise(example_struct));
  // check that we actually assigned values to the Vectors under the hood
  not_assigned = (!fv.x.has_elements() && !fv.y.has_elements());
  bool expected_size = (fv.x.size()==n_field_array_elements && fv.y.size()==n_field_array_elements);
  bool assigned_and_correct_size = ((!not_assigned) && expected_size);
  REQUIRE(assigned_and_correct_size);
  // and the values themselves are what we expect
  for (int i = 0; i < n_field_array_elements; i++) {
    CHECK(abs(fv.x[i] - 1. / ((double) (i + 1))) < tol);
    CHECK(abs(fv.y[i] + 1. / ((double) (i + 1))) < tol);
  }
  // cleanup
  mxDestroyArray(example_struct);
}
