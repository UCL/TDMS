/**
 * @file test_DispersiveMultiLayer.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests for the DispersiveMultiLayer class and its subclasses
 */
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"
#include "array_test_class.h"
#include "unit_test_utils.h"

using namespace std;
using Catch::Approx;
using tdms_tests::TOLERANCE;

void DispersiveMultilayerTest::test_empty_construction() {
  // Constructor should error if recieving a struct with no fields
  create_1by1_struct(0, {});
  REQUIRE_THROWS_AS(DispersiveMultiLayer(matlab_input), runtime_error);
}

void DispersiveMultilayerTest::test_wrong_input_type() {
  // Constructor should throw runtime_error at not recieving struct
  dimensions_2d[0] = 2; dimensions_2d[1] = 3;
  create_numeric_array(2, dimensions_2d, mxUINT16_CLASS);
  REQUIRE_THROWS_AS(DispersiveMultiLayer(matlab_input), runtime_error);
}

void DispersiveMultilayerTest::test_correct_construction() {
  create_1by1_struct(n_fields, fieldnames);
  // build "data" for each of the fields, which is going to be the same array filled with consecutive integers
  const int array_size[2] = {1, n_numeric_elements};
  mxArray *field_array_ptrs[n_fields];
  for (int i = 0; i < n_fields; i++) {
    field_array_ptrs[i] =
            mxCreateNumericArray(2, (const mwSize *) array_size, mxDOUBLE_CLASS, mxREAL);
    mxDouble *where_to_place_data = mxGetPr(field_array_ptrs[i]);
    for (int i = 0; i < n_numeric_elements; i++) { where_to_place_data[i] = (double) i; }
    mxSetField(matlab_input, 0, fieldnames[i], field_array_ptrs[i]);
  }
  // we should now be able to create a DispersiveMultiLayer object
  REQUIRE_NOTHROW(DispersiveMultiLayer(matlab_input));
  DispersiveMultiLayer dml(matlab_input);
  // now check that the data has been correctly assigned
  for (int i = 0; i < n_numeric_elements; i++) {
    CHECK(dml.alpha[i] == Approx(i).epsilon(TOLERANCE));
    CHECK(dml.beta[i] == Approx(i).epsilon(TOLERANCE));
    CHECK(dml.gamma[i] == Approx(i).epsilon(TOLERANCE));
    CHECK(dml.kappa.x[i] == Approx(i).epsilon(TOLERANCE));
    CHECK(dml.kappa.y[i] == Approx(i).epsilon(TOLERANCE));
    CHECK(dml.kappa.z[i] == Approx(i).epsilon(TOLERANCE));
    CHECK(dml.sigma.x[i] == Approx(i).epsilon(TOLERANCE));
    CHECK(dml.sigma.y[i] == Approx(i).epsilon(TOLERANCE));
    CHECK(dml.sigma.z[i] == Approx(i).epsilon(TOLERANCE));
  }
}

TEST_CASE("DispersiveMultiLayer") {
  DispersiveMultilayerTest().run_all_class_tests();
}
