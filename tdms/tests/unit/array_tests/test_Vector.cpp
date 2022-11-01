/**
 * @file test_Vector.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the Vector class and its subclasses
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"
#include "matlabio.h"

TEST_CASE("Vector: allocation and deallocation") {
  SPDLOG_INFO("== Testing Vector allocation/deallocation");
  // default constructor should not assign any elements or pointers
  Vector<double> v_double;
  REQUIRE(!v_double.has_elements());
  REQUIRE(v_double.size() == 0);

  // create a MATLAB array to cast to
  const int n_elements[1] = {8}; // MATLAB'S LOGIC IS BEYOND MY UNDERSTANDING, CAN'T have n_elements as an int, must be int*
  mxArray *array = mxCreateNumericArray(1, (const mwSize *) n_elements, mxUINT8_CLASS, mxREAL);
  REQUIRE( array != NULL );

  // assign this MATLAB array to a vector
  Vector<double> v_int(array);
  REQUIRE(v_int.has_elements());
  REQUIRE(v_int.size()==n_elements[0]);

  // cleanup our MATLAB array
  mxDestroyArray(array);
}
