/**
 * @file test_DTilde.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the DTilde class
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

using namespace std;

TEST_CASE("DTilde: allocation and deallocation") {
  SPDLOG_INFO("== Testing DTilde class");
  // no custom constructor
  DTilde dt;
  // check no information is contained in the DTilde object
  bool no_information_stored =
          (dt.num_det_modes() == 0) && (!dt.x.has_elements()) && (!dt.y.has_elements());
  REQUIRE(no_information_stored);

  // initialise() method should exit without assignment if we pass in a pointer to an empty array
  // create empty array
  const int empty_dims[2] = {0, 1};
  mxArray *empty_array =
          mxCreateNumericArray(2, (const mwSize *) empty_dims, mxUINT8_CLASS, mxREAL);
  // attempt assignment (2nd and 3rd args don't matter)
  dt.initialise(empty_array, 0, 0);
  // should still be unassigned vectors
  no_information_stored =
          (dt.num_det_modes() == 0) && (!dt.x.has_elements()) && (!dt.y.has_elements());
  REQUIRE(no_information_stored);
  // cleanup the empty array mess we just made
  mxDestroyArray(empty_array);

  // assignment will throw error if we attempt to provide a non-empty, non-struct array
  int dims[2] = {2, 3};
  mxArray *not_a_struct = mxCreateNumericArray(2, (const mwSize *) dims, mxUINT8_CLASS, mxREAL);
  CHECK_THROWS_AS(dt.initialise(not_a_struct, 0, 0), runtime_error);
  // cleanup
  mxDestroyArray(not_a_struct);

  // assignment will throw error if we attempt to provide a struct array that doesn't have two fields
  const char *too_few_names[1] = {"field1"};
  const char *too_mny_names[3] = {"field1", "field2", "field3"};
  mxArray *struct_with_too_few_fields =
          mxCreateStructArray(2, (const mwSize *) dims, 1, too_few_names);
  mxArray *struct_with_too_mny_fields =
          mxCreateStructArray(2, (const mwSize *) dims, 3, too_mny_names);
  CHECK_THROWS_AS(dt.initialise(struct_with_too_few_fields, 0, 0), runtime_error);
  CHECK_THROWS_AS(dt.initialise(struct_with_too_mny_fields, 0, 0), runtime_error);
  // cleanup
  mxDestroyArray(struct_with_too_few_fields);
  mxDestroyArray(struct_with_too_mny_fields);

  // otherwise, we need to provide a struct with two fields, Dx_tilde and Dy_tilde
  // these fields must contain (n_det_modes, n_rows, n_cols) arrays of doubles/complex
  const char *fieldnames[] = {"Dx_tilde", "Dy_tilde"};
  dims[0] = 1;
  dims[1] = 1;
  const int target_n_det_modes = 5, n_rows = 6, n_cols = 4;
  const int field_array_dims[3] = {target_n_det_modes, n_rows, n_cols};
  // create the struct
  mxArray *example_struct = mxCreateStructArray(2, (const mwSize *) dims, 2, fieldnames);
  // create the data for the fields of our struct
  mxArray *field_array_ptrs[2];
  for (int i = 0; i < 2; i++) {
    field_array_ptrs[i] =
            mxCreateNumericArray(3, (const mwSize *) field_array_dims, mxDOUBLE_CLASS, mxCOMPLEX);
    mxSetField(example_struct, 0, fieldnames[i], field_array_ptrs[i]);
  }
  // attempt to create a vector from this struct
  REQUIRE_NOTHROW(dt.initialise(example_struct, n_rows, n_cols));
  // check that we actually assigned values to the Vectors under the hood
  bool information_stored = (dt.num_det_modes() == target_n_det_modes) && (dt.x.has_elements()) &&
                            (dt.y.has_elements());
  CHECK(information_stored);
  // cleanup
  mxDestroyArray(example_struct);
}
