/**
 * @file test_FieldSample.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the FieldSample class and its subclasses
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

using namespace std;

TEST_CASE("FieldSample: allocation and deallocation") {
  SPDLOG_INFO("== Testing FieldSample class");

  // attempting to construct with an empty array should exit with default values
  const int empty_dims[2] = {0, 1};
  mxArray *empty_array =
          mxCreateNumericArray(2, (const mwSize *) empty_dims, mxUINT8_CLASS, mxREAL);
  // attempt assignment (2nd and 3rd args don't matter)
  FieldSample empty_fs(empty_array);
  // should still be unassigned vectors
  bool no_vectors_set = (!empty_fs.i.has_elements()) && (!empty_fs.j.has_elements()) && (!empty_fs.k.has_elements()) && (!empty_fs.n.has_elements());
  REQUIRE(no_vectors_set);
  // cleanup the empty array mess we just made
  mxDestroyArray(empty_array);

  // should really test whether the inputs being the wrong size/shape/etc are bad, but equally this is pointless, as we're just checking against the same things
  // as such, let's just skip to the correct constructor

  // expecting a struct with fields i, j, k, and n
  const char *fieldnames[] = {"i", "j", "k", "n"};
  const int struct_dims[2] = {1, 1};
  const int number_of_vector_elements = 8;
  const int vector_dims[2] = {1, number_of_vector_elements};
  // create the struct
  mxArray *populated_struct = mxCreateStructArray(2, (const mwSize *) struct_dims, 4, fieldnames);
  // create vectors to place in the fields of the struct
  mxArray *i_vector = mxCreateNumericArray(2, (const mwSize *) vector_dims, mxINT32_CLASS, mxREAL);
  mxSetField(populated_struct, 0, fieldnames[0], i_vector);
  mxArray *j_vector = mxCreateNumericArray(2, (const mwSize *) vector_dims, mxINT32_CLASS, mxREAL);
  mxSetField(populated_struct, 0, fieldnames[1], j_vector);
  mxArray *k_vector = mxCreateNumericArray(2, (const mwSize *) vector_dims, mxINT32_CLASS, mxREAL);
  mxSetField(populated_struct, 0, fieldnames[2], k_vector);
  mxArray *n_vector = mxCreateNumericArray(2, (const mwSize *) vector_dims, mxDOUBLE_CLASS, mxREAL);
  mxSetField(populated_struct, 0, fieldnames[3], n_vector);
  // populate the vectors with some information
  mxDouble *place_i_data = mxGetPr(i_vector); //< 0,1,2,...
  mxDouble *place_j_data = mxGetPr(j_vector); //< number_of_vector_elements-1, -2, -3...
  mxDouble *place_k_data = mxGetPr(k_vector); //< 0,1,2,...,N/2,0,1,2,...
  mxDouble *place_n_data = mxGetPr(n_vector); //< 1, 1/2, 1/3, ...
  for (int i = 0; i < 4; i++) {
    place_i_data[i] = i;
    place_j_data[i] = number_of_vector_elements - (i+1);
    place_k_data[i] = i%(number_of_vector_elements/2);
    place_n_data[i] = 1./(double)(i+1);
  }
  // we should be able to construct things now
  REQUIRE_NOTHROW(FieldSample(populated_struct));
  FieldSample fs(populated_struct);
  // check that the tensor attribute is the correct size
  const mwSize *dims = mxGetDimensions(fs.mx);
  for(int i = 0; i < 4; i++) {
    CHECK(dims[i] == number_of_vector_elements);
  }
  // cleanup
  mxDestroyArray(populated_struct);
  // FieldSample destructor will deallocate fs.tensor correctly

  // alternatively, we can pass in empty vectors to the constructor and it should create an empty 4D tensor
  const int empty_vector_dims[2] = {1, 0};
  // create the struct
  mxArray *struct_with_empty_vectors = mxCreateStructArray(2, (const mwSize *) struct_dims, 4, fieldnames);
  // create vectors to place in the fields of the struct
  mxArray *empty_i_vector = mxCreateNumericArray(2, (const mwSize *) empty_vector_dims, mxINT32_CLASS, mxREAL);
  mxArray *empty_j_vector = mxCreateNumericArray(2, (const mwSize *) empty_vector_dims, mxINT32_CLASS, mxREAL);
  mxArray *empty_k_vector = mxCreateNumericArray(2, (const mwSize *) empty_vector_dims, mxINT32_CLASS, mxREAL);
  mxArray *empty_n_vector = mxCreateNumericArray(2, (const mwSize *) empty_vector_dims, mxDOUBLE_CLASS, mxREAL);
  mxSetField(struct_with_empty_vectors, 0, fieldnames[0], empty_i_vector);
  mxSetField(struct_with_empty_vectors, 0, fieldnames[1], empty_j_vector);
  mxSetField(struct_with_empty_vectors, 0, fieldnames[2], empty_k_vector);
  mxSetField(struct_with_empty_vectors, 0, fieldnames[3], empty_n_vector);
  // we should be able to construct things now
  REQUIRE_NOTHROW(FieldSample(struct_with_empty_vectors));
  FieldSample fs_from_empty_vectors(struct_with_empty_vectors);
  // check that all vectors are empty flags as being true
  CHECK(!fs_from_empty_vectors.all_vectors_are_non_empty());
  // check that the tensor attribute is the correct size
  const mwSize *empty_tensor_dims = mxGetDimensions(fs_from_empty_vectors.mx);
  for (int i = 0; i < 4; i++) {
    CHECK(empty_tensor_dims[i] == 0);
  }
  // cleanup
  mxDestroyArray(struct_with_empty_vectors);
}
