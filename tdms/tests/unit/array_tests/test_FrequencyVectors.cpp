/**
 * @file test_FrequencyVectors.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests for the FrequencyVectors class and its subclasses
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"
#include "array_test_class.h"
#include "unit_test_utils.h"

using namespace std;
using tdms_tests::TOLERANCE;

void FrequencyVectorsTest::test_empty_construction() {
  // initialise() method should exit without assignment if we pass in a pointer to an empty array (regardless of whether this is a struct or not)
  FrequencyVectors fv;
  dimensions_2d[0] = 0;
  create_numeric_array(2, dimensions_2d, mxUINT8_CLASS);
  // attempt assignment
  fv.initialise(matlab_input);
  // should still be unassigned vectors
  bool not_assigned = (!fv.x.has_elements() && !fv.y.has_elements());
  REQUIRE(not_assigned);
}

void FrequencyVectorsTest::test_wrong_input_type() {
  FrequencyVectors fv;
  // throw error if we attempt to provide a non-empty, non-struct array
  dimensions_2d[0] = 2;
  dimensions_2d[1] = 3;
  create_numeric_array(2, dimensions_2d, mxUINT8_CLASS);
  REQUIRE_THROWS_AS(fv.initialise(matlab_input), runtime_error);
}

void FrequencyVectorsTest::test_incorrect_number_of_fields() {
  // assignment will throw error if we attempt to provide a struct array that doesn't have two fields
  FrequencyVectors fv;
  SECTION("Struct with too many inputs") {
    const char *too_many_names[3] = {"field1", "field2", "field3"};
    create_1by1_struct(3, too_many_names);
    CHECK_THROWS_AS(fv.initialise(matlab_input), runtime_error);
  }
  SECTION("Struct with too few inputs") {
    create_1by1_struct(n_fields - 1, fieldnames);
    CHECK_THROWS_AS(fv.initialise(matlab_input), runtime_error);
  }
}

void FrequencyVectorsTest::test_correct_construction() {
  FrequencyVectors fv;
  // members should start unassigned
  bool not_assigned = (!fv.x.has_elements() && !fv.y.has_elements());
  REQUIRE(not_assigned);
}

void FrequencyVectorsTest::test_initialise_method() {
  FrequencyVectors fv;
  const int field_array_dimensions[2] = {1, n_numeric_elements};
  // create the struct
  create_1by1_struct(n_fields, fieldnames);
  // create the data for the fields of our struct
  mxArray *field_array_ptrs[2];
  for (int i = 0; i < 2; i++) {
    field_array_ptrs[i] = mxCreateNumericArray(2, (const mwSize *) field_array_dimensions,
                                               mxDOUBLE_CLASS, mxREAL);
    mxDouble *where_to_place_data = mxGetPr(field_array_ptrs[i]);
    // 0th field, fx_vec[i], will be 1/(i+1)
    // 1st field, fy_vec[i], will be -1/(i+1)
    for (int j = 0; j < n_numeric_elements; j++) {
      where_to_place_data[j] = pow(-1., (double) i) / ((double) (j + 1));
    }
    mxSetField(matlab_input, 0, fieldnames[i], field_array_ptrs[i]);
  }
  // attempt to create a vector from this struct
  REQUIRE_NOTHROW(fv.initialise(matlab_input));
  // check that we actually assigned values to the Vectors under the hood
  bool not_assigned = (!fv.x.has_elements() && !fv.y.has_elements());
  bool expected_size =
          (fv.x.size() == n_numeric_elements && fv.y.size() == n_numeric_elements);
  bool assigned_and_correct_size = ((!not_assigned) && expected_size);
  REQUIRE(assigned_and_correct_size);
  // and the values themselves are what we expect
  bool values_are_correct = true;
  for (int i = 0; i < n_numeric_elements; i++) {
    values_are_correct = values_are_correct &&
                         (abs(fv.x[i] - 1. / ((double) (i + 1))) < TOLERANCE) &&
                         (abs(fv.y[i] + 1. / ((double) (i + 1))) < TOLERANCE);
  }
  REQUIRE(values_are_correct);
}

TEST_CASE("FrequencyVectors") {
  FrequencyVectorsTest().run_all_class_tests();
}
