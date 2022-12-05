/**
 * @file test_DTilde.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the DTilde class
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

using namespace std;

TEST_CASE("DTilde") {

  DTilde dt;
  bool no_information_stored;
  int dimensions[2];
  mxArray *test_input;

  // check no information is contained in the DTilde object
  SECTION("No information stored on declaration") {
    no_information_stored =
            (dt.num_det_modes() == 0) && (!dt.x.has_elements()) && (!dt.y.has_elements());
    REQUIRE(no_information_stored);
  }

  // initialise() method should exit without assignment if we pass in a pointer to an empty array
  // create empty array
  SECTION("Empty array input") {
    dimensions[0] = 0;
    dimensions[1] = 1;
    test_input = mxCreateNumericArray(2, (const mwSize *) dimensions, mxUINT8_CLASS, mxREAL);
    // attempt assignment (2nd and 3rd args don't matter)
    dt.initialise(test_input, 0, 0);
    // should still be unassigned vectors
    no_information_stored =
            (dt.num_det_modes() == 0) && (!dt.x.has_elements()) && (!dt.y.has_elements());
    REQUIRE(no_information_stored);
  }

  // assignment will throw error if we attempt to provide a non-empty, non-struct array
  SECTION("Non-struct input") {
    dimensions[0] = 2;
    dimensions[1] = 3;
    test_input = mxCreateNumericArray(2, (const mwSize *) dimensions, mxUINT8_CLASS, mxREAL);
    CHECK_THROWS_AS(dt.initialise(test_input, 0, 0), runtime_error);
  }

  // assignment will throw error if we attempt to provide a struct array that doesn't have two fields
  SECTION("Struct with too many fields") {
    dimensions[0] = 1;
    dimensions[1] = 1;
    const char *too_mny_names[3] = {"field1", "field2", "field3"};
    test_input = mxCreateStructArray(2, (const mwSize *) dimensions, 3, too_mny_names);
    CHECK_THROWS_AS(dt.initialise(test_input, 0, 0), runtime_error);
  }
  SECTION("Struct with too few fields") {
    dimensions[0] = 1;
    dimensions[1] = 1;
    const char *too_few_names[1] = {"field1"};
    test_input = mxCreateStructArray(2, (const mwSize *) dimensions, 1, too_few_names);
    CHECK_THROWS_AS(dt.initialise(test_input, 0, 0), runtime_error);
  }

  // otherwise, we need to provide a struct with two fields, Dx_tilde and Dy_tilde
  // these fields must contain (n_det_modes, n_rows, n_cols) arrays of doubles/complex
  SECTION("Correct input") {
    const char *fieldnames[] = {"Dx_tilde", "Dy_tilde"};
    dimensions[0] = 1;
    dimensions[1] = 1;
    const int target_n_det_modes = 5, n_rows = 6, n_cols = 4;
    const int field_array_dimensions[3] = {target_n_det_modes, n_rows, n_cols};
    // create the struct
    test_input = mxCreateStructArray(2, (const mwSize *) dimensions, 2, fieldnames);
    // create the data for the fields of our struct
    mxArray *field_array_ptrs[2];
    for (int i = 0; i < 2; i++) {
      field_array_ptrs[i] =
              mxCreateNumericArray(3, (const mwSize *) field_array_dimensions, mxDOUBLE_CLASS, mxCOMPLEX);
      mxSetField(test_input, 0, fieldnames[i], field_array_ptrs[i]);
    }
    // attempt to create a vector from this struct
    REQUIRE_NOTHROW(dt.initialise(test_input, n_rows, n_cols));
    // check that we actually assigned values to the Vectors under the hood
    bool information_stored = (dt.num_det_modes() == target_n_det_modes) && (dt.x.has_elements()) &&
                              (dt.y.has_elements());
    CHECK(information_stored);
  }

  // tear down
  mxDestroyArray(test_input);
}
