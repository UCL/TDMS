/**
 * @file test_DispersiveMultiLayer.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests for the DispersiveMultiLayer class and its subclasses
 */
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

using namespace std;
using Catch::Approx;

const double TOLERANCE = 1e-16;

TEST_CASE("DispersiveMultiLayer") {

  mxArray *matlab_input;
  // Constructor should throw runtime_error at not recieving struct
  // create a MATLAB pointer to something that isn't an array
  SECTION("Non-struct input") {
    int n_elements[2] = {2, 3};
    matlab_input = mxCreateNumericArray(2, (const mwSize *) n_elements, mxUINT8_CLASS, mxREAL);
    CHECK_THROWS_AS(DispersiveMultiLayer(matlab_input), runtime_error);
  }

  // Constructor should error if recieving an empty struct
  // create a MATLAB pointer to an empty struct
  SECTION("Empty struct input") {
    const char* empty_fields[] = {};
    const int empty_dims[2] = {1,1};
    matlab_input = mxCreateStructArray(2, (const mwSize *) empty_dims, 0, empty_fields);
    CHECK_THROWS_AS(DispersiveMultiLayer(matlab_input), runtime_error);
  }

  // For successful construction, we need to build a MATLAB struct with 9 fields
  // these are the fieldnames that are expected
  SECTION("Successful construction") {
    const int n_fields = 9;
    const char *fieldnames[n_fields] = {"alpha", "beta", "gamma", "kappa_x", "kappa_y",
                                        "kappa_z", "sigma_x", "sigma_y", "sigma_z"};
    const int n_field_elements = 5;
    // build our struct
    const int dims[2] = {1, 1};
    matlab_input = mxCreateStructArray(2, (const mwSize *) dims, n_fields, fieldnames);
    // build "data" for each of the fields, which is going to be the same array filled with consecutive integers
    const int array_size[2] = {1, n_field_elements};
    mxArray *field_array_ptrs[n_fields];
    for (int i = 0; i < n_fields; i++) {
      field_array_ptrs[i] =
              mxCreateNumericArray(2, (const mwSize *) array_size, mxDOUBLE_CLASS, mxREAL);
      mxDouble *where_to_place_data = mxGetPr(field_array_ptrs[i]);
      for (int i = 0; i < 5; i++) { where_to_place_data[i] = (double) i; }
      mxSetField(matlab_input, 0, fieldnames[i], field_array_ptrs[i]);
    }
    // we should now be able to create a DispersiveMultiLayer object
    REQUIRE_NOTHROW(DispersiveMultiLayer(matlab_input));
    DispersiveMultiLayer dml(matlab_input);
    // now check that the data has been correctly assigned
    for (int i = 0; i < 5; i++) {
      CHECK(dml.alpha[i] == Approx(i).epsilon(tol));
      CHECK(dml.beta[i] == Approx(i).epsilon(tol));
      CHECK(dml.gamma[i] == Approx(i).epsilon(tol));
      CHECK(dml.kappa.x[i] == Approx(i).epsilon(tol));
      CHECK(dml.kappa.y[i] == Approx(i).epsilon(tol));
      CHECK(dml.kappa.z[i] == Approx(i).epsilon(tol));
      CHECK(dml.sigma.x[i] == Approx(i).epsilon(tol));
      CHECK(dml.sigma.y[i] == Approx(i).epsilon(tol));
      CHECK(dml.sigma.z[i] == Approx(i).epsilon(tol));
    }
  }

  // cleanup MATLAB array
  mxDestroyArray(matlab_input);
}
