/**
 * @file test_CollectionBase.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for CollectionBase class and its subclasses (CCollection, DCollection)
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

TEST_CASE("CCollection") {

  // setup - an empty mxArray and the names of the field components
  mxArray *matlab_input;
  const int n_numeric_elements = 8;
  const char *fieldnames[] = {"Cax", "Cay", "Caz", "Cbx", "Cby", "Cbz", "Ccx", "Ccy", "Ccz"};

  // we need to provide a struct array with either 6 or 9 fields
  // incorrect number of fields should throw an error
  SECTION("Incorrect number of fields") {
    matlab_input = mxCreateStructMatrix(1, 1, 3, fieldnames);
    REQUIRE_THROWS_AS(CCollection(matlab_input), std::runtime_error);
  }

  // now create an array with the correct fieldnames
  SECTION("Expected input (6 fields)") {
    matlab_input = mxCreateStructMatrix(1, 1, 6, fieldnames);
    mxArray *elements6[6];
    for (int i = 0; i < 6; i++) {
      elements6[i] = mxCreateNumericMatrix(1, n_numeric_elements, mxDOUBLE_CLASS, mxREAL);
      mxSetField(matlab_input, 0, fieldnames[i], elements6[i]);
    }
    CCollection cc6(matlab_input);
    REQUIRE(!cc6.is_disp_ml);
    REQUIRE(!cc6.is_multilayer);
  }
  SECTION("Expected input (9 fields)") {
    matlab_input = mxCreateStructMatrix(1, 1, 9, fieldnames);
    mxArray *elements9[9];
    for (int i = 0; i < 9; i++) {
      elements9[i] = mxCreateNumericMatrix(2, n_numeric_elements, mxDOUBLE_CLASS, mxREAL);
      mxSetField(matlab_input, 0, fieldnames[i], elements9[i]);
    }
    CCollection cc9(matlab_input);
    REQUIRE(cc9.is_disp_ml);
    REQUIRE(cc9.is_multilayer);
  }

  // tear down
  mxDestroyArray(matlab_input);
}

TEST_CASE("DCollection") {

  mxArray *matlab_input;
  const char *fieldnames[] = {"Dax", "Day", "Daz", "Dbx", "Dby", "Dbz"};
  const int n_numeric_elements = 8;

  // we need to provide a struct array with either 6 or 9 fields
  // incorrect number of fields should throw an error
  SECTION("Incorrect number of fields") {
    matlab_input = mxCreateStructMatrix(1, 1, 3, fieldnames);
    REQUIRE_THROWS_AS(CCollection(matlab_input), std::runtime_error);
  }

  // now create an array with the correct fieldnames
  SECTION("Expected input") {
    matlab_input = mxCreateStructMatrix(1, 1, 6, fieldnames);
    mxArray *elements6[6];
    for (int i = 0; i < 6; i++) {
      elements6[i] = mxCreateNumericMatrix(1, n_numeric_elements, mxDOUBLE_CLASS, mxREAL);
      mxSetField(matlab_input, 0, fieldnames[i], elements6[i]);
    }
    REQUIRE_NOTHROW(DCollection(matlab_input));
  }

  // tear down - cleanup MATLAB memory
  mxDestroyArray(matlab_input);
};
