/**
 * @file test_Vector.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the Vector class
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "array_test_class.h"
#include "arrays.h"
#include "globals.h"
#include "unit_test_utils.h"

using namespace tdms_math_constants;
using tdms_tests::is_close;
using tdms_tests::TOLERANCE;

void VectorTest::test_correct_construction() {
  // default constructor should not assign any elements or pointers
  SECTION("Default constructor") {
    Vector<double> v_double;
    REQUIRE(!v_double.has_elements());
    REQUIRE(v_double.size() == 0);
  }
  SECTION("Overloaded constructor") {
    dimensions_2d[1] = n_numeric_elements;
    create_numeric_array(2, dimensions_2d);
    mxDouble *where_to_place_data = mxGetPr(matlab_input);
    // data is just the integers starting from 0
    for (int i = 0; i < n_numeric_elements; i++) {
      where_to_place_data[i] = (double) i;
    }
    // assign this MATLAB array to a vector
    Vector<double> v_from_matlab(matlab_input);
    bool has_elements_and_right_size =
            (v_from_matlab.has_elements() &&
             v_from_matlab.size() == n_numeric_elements);
    REQUIRE(has_elements_and_right_size);
    // check that the data is the _right_ data too
    bool correct_data_accessed = true;
    for (int i = 0; i < n_numeric_elements; i++) {
      correct_data_accessed =
              correct_data_accessed && is_close(v_from_matlab[i], (double) i);
    }
    REQUIRE(correct_data_accessed);
  }
}

TEST_CASE("Vector") { VectorTest().run_all_class_tests(); }
