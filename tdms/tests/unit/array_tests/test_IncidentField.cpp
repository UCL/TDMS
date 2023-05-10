/**
 * @file test_IncidentField.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the IncidentField class
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "array_test_class.h"
#include "arrays.h"
#include "unit_test_utils.h"

using namespace std;
using tdms_tests::TOLERANCE;

void IncidentFieldTest::test_empty_construction() {
  // construction should fail if given an empty array
  dimensions_2d[0] = 0;
  create_empty_struct();
  // attempt construction (2nd and 3rd args don't matter)
  REQUIRE_THROWS_AS(IncidentField(matlab_input), runtime_error);
}

void IncidentFieldTest::test_wrong_input_type() {
  dimensions_2d[0] = 2;
  dimensions_2d[1] = 3;
  create_numeric_array(2, dimensions_2d, mxUINT8_CLASS);
  CHECK_THROWS_AS(IncidentField(matlab_input), runtime_error);
}

void IncidentFieldTest::test_incorrect_number_of_fields() {
  SECTION("Struct with too many fields") {
    const char *too_many_names[3] = {"field1", "field2", "field3"};
    create_1by1_struct(3, too_many_names);
    CHECK_THROWS_AS(IncidentField(matlab_input), runtime_error);
  }
  SECTION("Struct with too few fields") {
    create_1by1_struct(n_fields - 1, fieldnames);
    CHECK_THROWS_AS(IncidentField(matlab_input), runtime_error);
  }
}
void IncidentFieldTest::test_correct_construction() {
  dimensions_3d[0] = n_rows;
  dimensions_3d[1] = n_cols;
  dimensions_3d[2] = n_layers;
  // create the struct
  create_1by1_struct(n_fields, fieldnames);
  // create the data for the fields of our struct
  mxArray *field_array_ptrs[2];
  for (int i = 0; i < 2; i++) {
    field_array_ptrs[i] = mxCreateNumericArray(
            3, (const mwSize *) dimensions_3d, mxDOUBLE_CLASS, mxREAL);
    mxDouble *place_data = mxGetPr(field_array_ptrs[i]);
    for (int ii = 0; ii < n_rows; ii++) {
      for (int jj = 0; jj < n_cols; jj++) {
        for (int kk = 0; kk < n_layers; kk++) {
          *(place_data + kk * n_cols * n_rows + jj * n_rows + ii) =
                  1. / ((double) (kk + jj + ii + 1));
        }
      }
    }
    mxSetField(matlab_input, 0, fieldnames[i], field_array_ptrs[i]);
  }
  // attempt to create a vector from this struct
  IncidentField i_field(matlab_input);
  // check that we actually assigned values to the Vectors under the hood
  bool information_stored =
          (i_field.x.has_elements()) && (i_field.y.has_elements());
  REQUIRE(information_stored);
  bool elements_set_correctly = true;
  for (int ii = 0; ii < n_rows; ii++) {
    for (int jj = 0; jj < n_cols; jj++) {
      for (int kk = 0; kk < n_layers; kk++) {
        elements_set_correctly =
                elements_set_correctly &&
                (abs(i_field.x(ii, jj, kk) -
                     1. / ((double) (kk + jj + ii + 1))) < TOLERANCE) &&
                (abs(i_field.y(ii, jj, kk) -
                     1. / ((double) (kk + jj + ii + 1))) < TOLERANCE);
      }
    }
  }
  REQUIRE(elements_set_correctly);
}

TEST_CASE("IncidentField") { IncidentFieldTest().run_all_class_tests(); }
