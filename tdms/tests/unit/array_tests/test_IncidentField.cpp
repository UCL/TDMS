/**
 * @file test_IncidentField.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the IncidentField class
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

using namespace std;

const double TOLERANCE = 1e-16;

TEST_CASE("IncidentField") {

  mxArray *matlab_input;
  int dims[2] = {1, 1};

  // construction should fail if given an empty array
  SECTION("Empty array input") {
    dims[0] = 0;
    matlab_input = mxCreateNumericArray(2, (const mwSize *) dims, mxUINT8_CLASS, mxREAL);
    // attempt construction (2nd and 3rd args don't matter)
    REQUIRE_THROWS_AS(IncidentField(matlab_input), runtime_error);
  }

  // construction fails if provided a non-empty, non-struct array
  SECTION("Non-struct input") {
    dims[0] = 2; dims[1] = 3;
    matlab_input = mxCreateNumericArray(2, (const mwSize *) dims, mxUINT8_CLASS, mxREAL);
    CHECK_THROWS_AS(IncidentField(matlab_input), runtime_error);
  }

  // construction fails if provided a struct array that doesn't have two fields
  SECTION("Struct with too many fields") {
    const char *too_mny_names[3] = {"field1", "field2", "field3"};
    matlab_input = mxCreateStructArray(2, (const mwSize *) dims, 3, too_mny_names);
    CHECK_THROWS_AS(IncidentField(matlab_input), runtime_error);
  }
  SECTION("Struct with too few fields") {
    const char *too_few_names[1] = {"field1"};
    matlab_input = mxCreateStructArray(2, (const mwSize *) dims, 1, too_few_names);
    CHECK_THROWS_AS(IncidentField(matlab_input), runtime_error);
  }

  // otherwise, we need to provide a struct with two fields, Dx_tilde and Dy_tilde
  // these fields must contain (n_det_modes, n_rows, n_cols) arrays of doubles/complex
  SECTION("Expected input") {
    const char *fieldnames[] = {"exi", "eyi"};
    const int n_rows = 6, n_cols = 4, n_layers = 5;
    const int field_array_dims[3] = {n_rows, n_cols, n_layers};
    // create the struct
    matlab_input = mxCreateStructArray(2, (const mwSize *) dims, 2, fieldnames);
    // create the data for the fields of our struct
    mxArray *field_array_ptrs[2];
    for (int i = 0; i < 2; i++) {
      field_array_ptrs[i] =
              mxCreateNumericArray(3, (const mwSize *) field_array_dims, mxDOUBLE_CLASS, mxREAL);
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
    bool information_stored = (i_field.x.has_elements()) && (i_field.y.has_elements());
    REQUIRE(information_stored);
    for (int ii = 0; ii < n_rows; ii++) {
      for (int jj = 0; jj < n_cols; jj++) {
        for (int kk = 0; kk < n_layers; kk++) {
          bool element_set_correctly =
                  (abs(i_field.x[kk][jj][ii] - 1. / ((double) (kk + jj + ii + 1))) < tol) &&
                  (abs(i_field.y[kk][jj][ii] - 1. / ((double) (kk + jj + ii + 1))) < tol);
          CHECK(element_set_correctly);
        }
      }
    }
  }

  // cleanup
  mxDestroyArray(matlab_input);
}
