/**
 * @file test_VertexPhasors.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the VertexPhasors class
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "array_test_class.h"
#include "globals.h"
#include "vertex_phasors.h"

using namespace std;

void VertexPhasorsTest::test_empty_construction() {
  create_empty_struct();
  CHECK_NOTHROW(VertexPhasors(matlab_input));
  VertexPhasors empty_test(matlab_input);
  // n_vertices should return 0, since Vector has not been set so should default
  // initialise to 0-vertex matrix vector & components should be flagged as
  // having no elements the object should also flag that there are no vertices
  // to extract at (since there are 0 vertices listed)
  CHECK(empty_test.n_vertices() == 0);
  CHECK(!empty_test.there_are_vertices_to_extract_at());
  CHECK(!empty_test.there_are_elements_in_arrays());
}

void VertexPhasorsTest::test_correct_construction() {
  // dummy array dimensions
  int n_components = 4;
  int f_ex_vec_size = 5;

  create_struct_array(2, dimensions_2d, 2, fieldnames);
  // each entry is a numeric array, setup for vertices and components
  // array for "vertices" field
  mxArray *vertices_array =
          mxCreateNumericMatrix(n_numeric_elements, 3, mxINT32_CLASS, mxREAL);
  // array for "components" field
  mxArray *components_vector =
          mxCreateNumericMatrix(1, n_components, mxINT16_CLASS, mxREAL);
  // add fields to struct that will be passed to CAS constructor function
  mxSetField(matlab_input, 0, fieldnames[0], vertices_array);
  mxSetField(matlab_input, 0, fieldnames[1], components_vector);
  // create an instance using this struct
  VertexPhasors vp(matlab_input);
  // check that we are flagged as having elements to extract at
  CHECK(vp.there_are_vertices_to_extract_at());
  // attempt to initialise camplitude{R,I} arrays
  vp.setup_complex_amplitude_arrays(f_ex_vec_size);
}

TEST_CASE("VertexPhasors") { VertexPhasorsTest().run_all_class_tests(); }
