/**
 * @file test_CollectionBase.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for CollectionBase class and its subclasses (CCollection,
 * DCollection)
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "array_test_class.h"
#include "arrays.h"

void CCollectionTest::test_incorrect_number_of_fields() {
  create_1by1_struct(3, fieldnames);
  REQUIRE_THROWS_AS(CCollection(matlab_input), std::runtime_error);
}

void CCollectionTest::test_correct_construction() {
  SECTION("(6 inputs)") {
    create_1by1_struct(6, fieldnames);
    mxArray *field_arrays[6];
    for (int i = 0; i < 6; i++) {
      field_arrays[i] = mxCreateNumericMatrix(1, n_numeric_elements,
                                              mxDOUBLE_CLASS, mxREAL);
      mxSetField(matlab_input, 0, fieldnames[i], field_arrays[i]);
    }
    CCollection cc6(matlab_input);
    REQUIRE(!cc6.is_disp_ml);
    REQUIRE(!cc6.is_multilayer);
  }
  SECTION("(9 inputs)") {
    create_1by1_struct(9, fieldnames);
    mxArray *field_arrays[9];
    for (int i = 0; i < 9; i++) {
      field_arrays[i] = mxCreateNumericMatrix(2, n_numeric_elements,
                                              mxDOUBLE_CLASS, mxREAL);
      mxSetField(matlab_input, 0, fieldnames[i], field_arrays[i]);
    }
    CCollection cc9(matlab_input);
    REQUIRE(cc9.is_disp_ml);
    REQUIRE(cc9.is_multilayer);
  }
}

void DCollectionTest::test_incorrect_number_of_fields() {
  create_1by1_struct(3, fieldnames);
  REQUIRE_THROWS_AS(DCollection(matlab_input), std::runtime_error);
}

void DCollectionTest::test_correct_construction() {
  create_1by1_struct(6, fieldnames);
  mxArray *elements6[6];
  for (int i = 0; i < 6; i++) {
    elements6[i] = mxCreateNumericMatrix(1, n_numeric_elements, mxDOUBLE_CLASS,
                                         mxREAL);
    mxSetField(matlab_input, 0, fieldnames[i], elements6[i]);
  }
  REQUIRE_NOTHROW(DCollection(matlab_input));
}

TEST_CASE("CCollection") { CCollectionTest().run_all_class_tests(); }

TEST_CASE("DCollection") { DCollectionTest().run_all_class_tests(); };
