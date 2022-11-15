/**
 * @file test_Material.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the Material class and it's subclasses (CMaterial, DMaterial)
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

TEST_CASE("CMaterial") {

  mxArray *matlab_input;
  // the constructor expects a struct of 9 fields
  const int n_fields = 9;
  const char *fieldnames[n_fields] = {"Cax", "Cay", "Caz", "Cbx", "Cby", "Cbz", "Ccx", "Ccy", "Ccz"};
  const int n_numeric_elements = 8;

  SECTION("Bad fieldname") {
    const char *wrong_fieldnames[] = {"Cax", "Cay", "BAD", "Cbx", "Cby",
                                      "Cbz", "Ccx", "Ccy", "Ccz"};
    matlab_input = mxCreateStructMatrix(1, 1, 9, wrong_fieldnames);
    REQUIRE_THROWS_AS(CMaterial(matlab_input), std::runtime_error);
  }
  SECTION("Wrong number of fields") {
    matlab_input = mxCreateStructMatrix(1, 1, n_fields - 1, fieldnames);
    REQUIRE_THROWS_AS(CMaterial(matlab_input), std::runtime_error);
  }

  // try and construct something meaningful
  SECTION("Expected input") {
    matlab_input = mxCreateStructMatrix(1, 1, 9, fieldnames);
    mxArray *elements9[9];
    for (int i = 0; i < 9; i++) {
      elements9[i] = mxCreateNumericMatrix(1, n_numeric_elements, mxDOUBLE_CLASS, mxREAL);
      mxSetField(matlab_input, 0, fieldnames[i], elements9[i]);
    }
    // construction should succeed
    CMaterial cm(matlab_input);
  }

  mxDestroyArray(matlab_input);
}

TEST_CASE("DMaterial") {

  mxArray *matlab_input;
  // the constructor expects a struct of 6 fields
  const int n_fields = 6;
  const char *fieldnames[n_fields] = {"Dax", "Day", "Daz", "Dbx", "Dby", "Dbz"};
  const int n_numeric_elements = 8;

  SECTION("Bad fieldname") {
    const char *wrong_fieldnames[] = {"Cax", "Day", "Daz", "Dbx", "Cby", "Dbz"};
    matlab_input = mxCreateStructMatrix(1, 1, 6, wrong_fieldnames);
    REQUIRE_THROWS_AS(DMaterial(matlab_input), std::runtime_error);
  }
  SECTION("Wrong number of fields") {
    matlab_input = mxCreateStructMatrix(1, 1, n_fields - 1, fieldnames);
    REQUIRE_THROWS_AS(DMaterial(matlab_input), std::runtime_error);
  }

  SECTION("Expected input") {
    matlab_input = mxCreateStructMatrix(1, 1, 6, fieldnames);
    mxArray *elements6[6];
    for (int i = 0; i < 6; i++) {
      elements6[i] = mxCreateNumericMatrix(1, n_numeric_elements, mxDOUBLE_CLASS, mxREAL);
      mxSetField(matlab_input, 0, fieldnames[i], elements6[i]);
    }
    // construction should succeed
    DMaterial dm(matlab_input);
  }

  mxDestroyArray(matlab_input);
}
