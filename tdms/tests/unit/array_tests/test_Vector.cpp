/**
 * @file test_Vector.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the Vector class and its subclasses (FieldComponentsVector, FrequencyExtractVector)
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"
#include "globals.h"

using namespace tdms_math_constants;

const double TOLERANCE = 1e-16;

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

TEST_CASE("FieldComponentsVector") {

  // FieldComponentsVector can be initialised with the default constructor, then have initialise() called to assign values from a pre-existing MATLAB array
  FieldComponentsVector fcv;
  bool fcv_ready = ((!fcv.has_elements()) && (fcv.size() == 0));
  REQUIRE(fcv_ready);

  const int n_vector_elements = 8;
  const char *fieldnames[1] = {"components"};
  mxArray *matlab_input;

  // create a MATLAB struct with a "components" field to cast to
  // the constructor shouldn't care whether or not the components field contains an n*1 or 1*n vector either
  // struct with horizontal vector
  SECTION("Expected input") {
    matlab_input = mxCreateStructMatrix(1, 1, 1, fieldnames);
    mxArray *vector_array;
    SECTION("(horz)") {
      vector_array = mxCreateNumericMatrix(1, n_vector_elements, mxINT16_CLASS, mxREAL);
    }
    SECTION("(vert)") {
      vector_array = mxCreateNumericMatrix(n_vector_elements, 1, mxINT16_CLASS, mxREAL);
    }
    mxSetField(matlab_input, 0, fieldnames[0], vector_array);
    // initialise
    REQUIRE_NOTHROW(fcv.initialise(matlab_input));
    // number of elements should be n_vector_elements
    REQUIRE(fcv.size() == n_vector_elements);
  }

  // destroy what we created
  mxDestroyArray(matlab_input);
}

TEST_CASE("FrequencyExtractVector") {

  const double omega_an = 1;
  const int n_elements = 8;
  mxArray *matlab_input;

  // if passed in an empty pointer, creates a length-1 array with a single, predetermined element
  SECTION("Empty input") {
    matlab_input = mxCreateNumericMatrix(0, n_elements, mxDOUBLE_CLASS, mxREAL);
    FrequencyExtractVector fev_empty(matlab_input, omega_an);
    CHECK(fev_empty.has_elements());
    CHECK(fev_empty.size() == 1);
    CHECK(fev_empty[0] == omega_an / 2. / DCPI);
  }

  // catch wrong-dimension inputs
  SECTION("Wrong number of dimensions") {
    SECTION("(2D array)") {
      matlab_input = mxCreateNumericMatrix(2, 4, mxDOUBLE_CLASS, mxREAL);
      REQUIRE_THROWS_AS(FrequencyExtractVector(matlab_input, omega_an), std::runtime_error);
    }
    SECTION("(3D array)") {
      int dims_3d[3] = {1, n_elements, 3};
      matlab_input = mxCreateNumericArray(3, (mwSize *) dims_3d, mxDOUBLE_CLASS, mxREAL);
      REQUIRE_THROWS_AS(FrequencyExtractVector(matlab_input, omega_an), std::runtime_error);
    }
  }

  // otherwise, it shouldn't matter whether we pass in a vertical or horizontal array in
  SECTION("Expected input") {
    SECTION("(horz)") {
      matlab_input = mxCreateNumericMatrix(1, n_elements, mxDOUBLE_CLASS, mxREAL);
      double *horz_data_placement = (double *) mxGetPr(matlab_input);
      for (int i = 0; i < n_elements; i++) {
        horz_data_placement[i] = ((double) i / (double) (n_elements - 1)) * omega_an * DCPI;
      }
        FrequencyExtractVector fev_horz(matlab_input, omega_an);
        CHECK(fev_horz.size() == n_elements);
        CHECK(abs(fev_horz.max() - fev_horz[n_elements - 1]) < TOLERANCE);
    }
    SECTION("(vert)") {
      matlab_input = mxCreateNumericMatrix(n_elements, 1, mxDOUBLE_CLASS, mxREAL);
      // insert some data into the arrays
      double *vert_data_placement = (double *) mxGetPr(matlab_input);
      for (int i = 0; i < n_elements; i++) {
        vert_data_placement[i] =
                ((double) (n_elements - i) / (double) (n_elements - 1)) * omega_an * DCPI;
      }
      // create the arrays
      FrequencyExtractVector fev_vert(matlab_input, omega_an);
      CHECK(fev_vert.size() == n_elements);
      // max should return the element at index n_elements-1 for horz and 0 for vert
      CHECK(abs(fev_vert.max() - fev_vert[0]) < TOLERANCE);
    }
  }

  // cleanup memory
  mxDestroyArray(matlab_input);
}
