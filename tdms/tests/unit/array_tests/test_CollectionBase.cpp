/**
 * @file test_CollectionBase.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for CollectionBase class and its subclasses ()
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

TEST_CASE("CCollection"){
    SPDLOG_INFO("== Testing CCollection class");

    // we need to provide a struct array with either 6 or 9 fields
    // incorrect number of fields should throw an error
    const char *bad_fieldnames[] = {"1", "2", "3"};
    mxArray *bad_struct = mxCreateStructMatrix(1, 1, 3, bad_fieldnames);
    REQUIRE_THROWS_AS(CCollection(bad_struct), std::runtime_error);
    mxDestroyArray(bad_struct);

    // now create an array with the correct fieldnames
    const char *fieldnames[] = {"Cax", "Cay", "Caz", "Cbx", "Cby", "Cbz", "Ccx", "Ccy", "Ccz"};
    const int n_numeric_elements = 8;
    mxArray *struct_6 = mxCreateStructMatrix(1, 1, 6, fieldnames);
    mxArray *struct_9 = mxCreateStructMatrix(1, 1, 9, fieldnames);
    mxArray *elements9[9];
    mxArray *elements6[6];
    for (int i = 0; i < 6; i++) {
      elements6[i] = mxCreateNumericMatrix(1, n_numeric_elements, mxDOUBLE_CLASS, mxREAL);
      mxSetField(struct_6, 0, fieldnames[i], elements6[i]);
      elements9[i] = mxCreateNumericMatrix(1, n_numeric_elements, mxDOUBLE_CLASS, mxREAL);
      mxSetField(struct_9, 0, fieldnames[i], elements9[i]);
    }
    // make the last set 2 by n_numeric_elements so as to check is_multilayer is flagged
    for (int i = 6; i < 9; i++) {
      elements9[i] = mxCreateNumericMatrix(2, n_numeric_elements, mxDOUBLE_CLASS, mxREAL);
      mxSetField(struct_9, 0, fieldnames[i], elements9[i]);
    }
    CCollection cc6(struct_6), cc9(struct_9);
    REQUIRE(!cc6.is_disp_ml);
    REQUIRE(cc9.is_disp_ml);
    REQUIRE(!cc6.is_multilayer);
    REQUIRE(cc9.is_multilayer);

    // cleanup MATLAB memory (these are recursive so elements{6,9} also destroyed)
    mxDestroyArray(struct_6);
    mxDestroyArray(struct_9);
};
