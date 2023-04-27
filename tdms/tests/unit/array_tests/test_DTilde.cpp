/**
 * @file test_DTilde.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the DTilde class
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "array_test_class.h"
#include "arrays.h"

using namespace std;

void DTildeTest::test_correct_construction() {
  // check no information is contained in the DTilde object when declared with
  // default constructor
  DTilde dt;
  bool no_information_stored = (dt.num_det_modes() == 0) &&
                               (!dt.x.has_elements()) && (!dt.y.has_elements());
  REQUIRE(no_information_stored);
}

void DTildeTest::test_empty_construction() {
  DTilde dt;
  bool no_information_stored;
  // initialise() method should exit without assignment
  create_empty_struct();
  // attempt assignment (2nd and 3rd args don't matter)
  dt.initialise(matlab_input, 0, 0);
  // should still be unassigned vectors
  no_information_stored = (dt.num_det_modes() == 0) && (!dt.x.has_elements()) &&
                          (!dt.y.has_elements());
  REQUIRE(no_information_stored);
}

void DTildeTest::test_wrong_input_type() {
  DTilde dt;
  // initialise() will throw error if we attempt to provide a non-empty,
  // non-struct array
  dimensions_2d[0] = 2;
  dimensions_2d[1] = 3;
  create_numeric_array(2, dimensions_2d, mxUINT8_CLASS);
  REQUIRE_THROWS_AS(dt.initialise(matlab_input, 0, 0), runtime_error);
}

void DTildeTest::test_incorrect_number_of_fields() {
  DTilde dt;
  // assignment will throw error if we attempt to provide a struct array that
  // doesn't have two fields
  SECTION("Struct with too many fields") {
    const char *too_mny_names[3] = {"field1", "field2", "field3"};
    create_struct_array(2, dimensions_2d, 3, too_mny_names);
    REQUIRE_THROWS_AS(dt.initialise(matlab_input, 0, 0), runtime_error);
  }
  SECTION("Struct with too few fields") {
    SKIP("Causes segmentation violation.");
    const char *too_few_names[1] = {"field1"};
    create_struct_array(2, dimensions_2d, 3, too_few_names);
    REQUIRE_THROWS_AS(dt.initialise(matlab_input, 0, 0), runtime_error);
  }
}

void DTildeTest::test_initialise_method() {
  DTilde dt;
  // otherwise, we need to provide a struct with two fields, Dx_tilde and
  // Dy_tilde these fields must contain (n_det_modes, n_rows, n_cols) arrays of
  // doubles/complex
  const int target_n_det_modes = 5, n_rows = 6, n_cols = 4;
  const int field_array_dimensions[3] = {target_n_det_modes, n_rows, n_cols};
  // create the struct
  create_struct_array(2, dimensions_2d, n_fields, fieldnames);
  // create the data for the fields of our struct
  mxArray *field_array_ptrs[n_fields];
  for (int i = 0; i < n_fields; i++) {
    field_array_ptrs[i] =
            mxCreateNumericArray(3, (const mwSize *) field_array_dimensions,
                                 mxDOUBLE_CLASS, mxCOMPLEX);
    mxSetField(matlab_input, 0, fieldnames[i], field_array_ptrs[i]);
  }
  // attempt to create a vector from this struct
  REQUIRE_NOTHROW(dt.initialise(matlab_input, n_rows, n_cols));
  // check that we actually assigned values to the Vectors under the hood
  bool information_stored = (dt.num_det_modes() == target_n_det_modes) &&
                            (dt.x.has_elements()) && (dt.y.has_elements());
  CHECK(information_stored);
}

TEST_CASE("DTilde") { DTildeTest().run_all_class_tests(); }
