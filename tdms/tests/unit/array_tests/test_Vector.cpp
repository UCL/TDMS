/**
 * @file test_Vector.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the Vector class and its subclasses
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"
#include "matlabio.h"

TEST_CASE("Vector") {

  // default constructor should not assign any elements or pointers
  Vector<double> v_double;
  REQUIRE(!v_double.has_elements());
  REQUIRE(v_double.size() == 0);

  // create a MATLAB array to cast to
  const int n_elements = 8;
  const int dims[2] = {1, n_elements};
  mxArray *array = mxCreateNumericArray(2, (const mwSize *) dims, mxDOUBLE_CLASS, mxREAL);
  mxDouble *where_to_place_data = mxGetPr(array);
  // data is just the integers starting from 0
  for (int i = 0; i < n_elements; i++) { where_to_place_data[i] = (double) i; }

  // assign this MATLAB array to a vector
  Vector<double> v_from_matlab(array);
  bool has_elements_and_right_size =
          (v_from_matlab.has_elements() && v_from_matlab.size() == n_elements);
  REQUIRE(has_elements_and_right_size);
  // check that the data is the _right_ data too
  for (int i = 0; i < n_elements; i++) { CHECK(v_from_matlab[i] == i); }

  // cleanup our MATLAB array
  mxDestroyArray(array);
}
