/**
 * @file test_Vector.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the Vector class and its subclasses (FieldComponentsVector, FrequencyExtractVector)
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
    for (int i = 0; i < n_numeric_elements; i++) { where_to_place_data[i] = (double) i; }
    // assign this MATLAB array to a vector
    Vector<double> v_from_matlab(matlab_input);
    bool has_elements_and_right_size =
            (v_from_matlab.has_elements() && v_from_matlab.size() == n_numeric_elements);
    REQUIRE(has_elements_and_right_size);
    // check that the data is the _right_ data too
    bool correct_data_accessed = true;
    for (int i = 0; i < n_numeric_elements; i++) {
      correct_data_accessed = correct_data_accessed && is_close(v_from_matlab[i], (double) i);
    }
    REQUIRE(correct_data_accessed);
  }
}

void FieldComponentsVectorTest::test_correct_construction() {
  // FieldComponentsVector can be initialised with the default constructor, then have initialise() called to assign values from a pre-existing MATLAB array
  FieldComponentsVector fcv;
  bool fcv_ready = ((!fcv.has_elements()) && (fcv.size() == 0));
  REQUIRE(fcv_ready);
}

void FieldComponentsVectorTest::test_initialise_method() {
  FieldComponentsVector fcv;
  // create a MATLAB struct with a "components" field to cast to
  // the constructor shouldn't care whether or not the components field contains an n*1 or 1*n vector either
  // struct with horizontal vector
  create_1by1_struct(n_fields, fieldnames);
  mxArray *vector_array;
  SECTION("(horz)") {
    vector_array = mxCreateNumericMatrix(1, n_numeric_elements, mxINT16_CLASS, mxREAL);
  }
  SECTION("(vert)") {
    vector_array = mxCreateNumericMatrix(n_numeric_elements, 1, mxINT16_CLASS, mxREAL);
  }
  mxSetField(matlab_input, 0, fieldnames[0], vector_array);
  // initialise
  REQUIRE_NOTHROW(fcv.initialise(matlab_input));
  // number of elements should be n_vector_elements
  REQUIRE(fcv.size() == n_numeric_elements);
}

void FrequencyExtractVectorTest::test_empty_construction() {
  // if passed in an empty pointer, creates a length-1 array with a single, predetermined element
  dimensions_2d[0] = 0;
  dimensions_2d[1] = n_numeric_elements;
  create_numeric_array(2, dimensions_2d);
  FrequencyExtractVector fev_empty(matlab_input, omega_an);
  CHECK(fev_empty.has_elements());
  CHECK(fev_empty.size() == 1);
  CHECK(fev_empty[0] == omega_an / 2. / DCPI);
}

void FrequencyExtractVectorTest::test_wrong_input_dimensions() {
  SECTION("(2D array)") {
    dimensions_2d[0] = 2;
    dimensions_2d[1] = 4;
    create_numeric_array(2, dimensions_2d);
  }
  SECTION("(3D array)") {
    dimensions_3d[0] = 1;
    dimensions_3d[1] = n_numeric_elements;
    dimensions_3d[2] = 3;
    create_numeric_array(3, dimensions_3d);
  }
  REQUIRE_THROWS_AS(FrequencyExtractVector(matlab_input, omega_an), std::runtime_error);
}

void FrequencyExtractVectorTest::test_correct_construction() {
  SECTION("(horz)") { dimensions_2d[1] = n_numeric_elements; }
  SECTION("(vert)") { dimensions_2d[0] = n_numeric_elements; }
  create_numeric_array(2, dimensions_2d);
  FrequencyExtractVector fev(matlab_input, omega_an);
  CHECK(fev.size() == n_numeric_elements);
}

void FrequencyExtractVectorTest::test_other_methods() {
  SECTION("max()") {
    dimensions_2d[1] = n_numeric_elements;
    create_numeric_array(2, dimensions_2d);
    double *data_placement = (double *) mxGetPr(matlab_input);
    int target_max_index;
    SECTION("Find element n_numeric_elements-1") {
      for (int i = 0; i < n_numeric_elements; i++) {
        data_placement[i] = ((double) i / (double) (n_numeric_elements - 1)) * omega_an * DCPI;
      }
      target_max_index = n_numeric_elements - 1;
    }
    SECTION("Find element 0") {
      for (int i = 0; i < n_numeric_elements; i++) {
        data_placement[i] =
                ((double) (n_numeric_elements - i) / (double) (n_numeric_elements - 1)) * omega_an *
                DCPI;
      }
      target_max_index = 0;
    }
    // create array
    FrequencyExtractVector fev(matlab_input, omega_an);
    // check max identifies the correct value in fev
    REQUIRE(is_close(fev.max(), fev[target_max_index]));
  }
}

TEST_CASE("Vector") { VectorTest().run_all_class_tests(); }

TEST_CASE("FieldComponentsVector") { FieldComponentsVectorTest().run_all_class_tests(); }

TEST_CASE("FrequencyExtractVector") { FrequencyExtractVectorTest().run_all_class_tests(); }
