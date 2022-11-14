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

const double tol = 1e-16;

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
  SPDLOG_INFO("== FieldComponentsVector");

  // FieldComponentsVector can be initialised with the default constructor, then have initialise() called to assign values from a pre-existing MATLAB array
  FieldComponentsVector v_horz;
  bool v_horz_ready = ((!v_horz.has_elements()) && (v_horz.size() == 0));
  REQUIRE(v_horz_ready);

  // create a MATLAB struct with a "components" field to cast to
  // the constructor shouldn't care whether or not the components field contains an n*1 or 1*n vector either
  // struct with horizontal vector
  const int n_fields_horz = 1;
  const int n_vector_elements_horz = 8;
  const int horz_struct_dims[2] = {1, 1};
  const char *horz_struct_fieldnames[n_fields_horz] = {"components"};
  mxArray *fcv_struct_pointer_horz = mxCreateStructArray(2, (const mwSize *) horz_struct_dims,
                                                         n_fields_horz, horz_struct_fieldnames);
  const int vector_dims_horz[2] = {1, n_vector_elements_horz};
  mxArray *horz_vector =
          mxCreateNumericArray(2, (const mwSize *) vector_dims_horz, mxINT16_CLASS, mxREAL);
  mxSetField(fcv_struct_pointer_horz, 0, horz_struct_fieldnames[0], horz_vector);
  // initialise
  v_horz.initialise(fcv_struct_pointer_horz);
  // in both cases, the number of elements should be n_vector_elements
  REQUIRE(v_horz.size() == n_vector_elements_horz);

  // destroy what we created
  mxDestroyArray(fcv_struct_pointer_horz);

  // now do the same for a vertical array
  FieldComponentsVector v_vert;
  bool v_vert_ready = ((!v_vert.has_elements()) && (v_vert.size() == 0));
  REQUIRE(v_vert_ready);

  // struct with vertical vector
  const int n_fields_vert = 1;
  const int n_vector_elements_vert = 8;
  const int vert_struct_dims[2] = {1, 1};
  const char *vert_struct_fieldnames[n_fields_vert] = {"components"};
  mxArray *fcv_struct_pointer_vert = mxCreateStructArray(2, (const mwSize *) vert_struct_dims,
                                                         n_fields_vert, vert_struct_fieldnames);
  const int vector_dims_vert[2] = {n_vector_elements_vert, 1};
  mxArray *vert_vector =
          mxCreateNumericArray(2, (const mwSize *) vector_dims_vert, mxINT16_CLASS, mxREAL);
  mxSetField(fcv_struct_pointer_vert, 0, vert_struct_fieldnames[0], vert_vector);
  // initialise
  v_vert.initialise(fcv_struct_pointer_vert);
  // in both cases, the number of elements should be n_vector_elements
  REQUIRE(v_vert.size() == n_vector_elements_vert);

  // destroy what we created
  mxDestroyArray(fcv_struct_pointer_vert);
}

TEST_CASE("FrequencyExtractVector") {
  SPDLOG_INFO("== FrequencyExtractVector");

  const double omega_an = 1;
  const int n_elements = 8;
  // if passed in an empty pointer, creates a length-1 array with a single, predetermined element
  mxArray *empty_array = mxCreateNumericMatrix(0, n_elements, mxDOUBLE_CLASS, mxREAL);
  FrequencyExtractVector fev_empty(empty_array, omega_an);
  CHECK(fev_empty.has_elements());
  CHECK(fev_empty.size()==1);
  CHECK(fev_empty[0] == omega_an / 2. / DCPI);
  mxDestroyArray(empty_array);

  // catch wrong-dimension inputs
  // input is not a vector
  mxArray *not_1D = mxCreateNumericMatrix(2, 4, mxDOUBLE_CLASS, mxREAL);
  REQUIRE_THROWS_AS(FrequencyExtractVector(not_1D, omega_an), std::runtime_error);
  mxDestroyArray(not_1D);
  // input is 3D array
  int dims_3d[3] = {1, n_elements, 3};
  mxArray *too_many_dims = mxCreateNumericArray(3, (mwSize *) dims_3d, mxDOUBLE_CLASS, mxREAL);
  REQUIRE_THROWS_AS(FrequencyExtractVector(too_many_dims, omega_an), std::runtime_error);
  mxDestroyArray(too_many_dims);

  // otherwise, it shouldn't matter whether we pass in a vertical or horizontal array in
  mxArray *horz_array = mxCreateNumericMatrix(1, n_elements, mxDOUBLE_CLASS, mxREAL);
  mxArray *vert_array = mxCreateNumericMatrix(n_elements, 1, mxDOUBLE_CLASS, mxREAL);
  // insert some data into the arrays
  double *horz_data_placement = (double *) mxGetPr(horz_array);
  double *vert_data_placement = (double *) mxGetPr(vert_array);
  for (int i = 0; i < n_elements; i++) {
    horz_data_placement[i] = ((double) i / (double) (n_elements - 1)) * omega_an * DCPI;
    vert_data_placement[i] = ((double) (n_elements - i) / (double) (n_elements - 1)) * omega_an * DCPI;
  }
  // create the arrays
  FrequencyExtractVector fev_horz(horz_array, omega_an);
  CHECK(fev_horz.size() == n_elements);
  FrequencyExtractVector fev_vert(vert_array, omega_an);
  CHECK(fev_vert.size() == n_elements);
  // max should return the element at index n_elements-1 for horz and 0 for vert
  CHECK(abs(fev_horz.max() - fev_horz[n_elements - 1]) < tol);
  CHECK(abs(fev_vert.max() - fev_vert[0]) < tol);
  // cleanup memory
  mxDestroyArray(horz_array);
  mxDestroyArray(vert_array);
}
