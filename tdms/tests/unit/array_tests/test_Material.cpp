/**
 * @file test_Material.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the Material class and it's subclasses (CMaterial, DMaterial)
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "array_test_class.h"
#include "arrays.h"

void CMaterialTest::test_incorrect_number_of_fields() {
  SECTION("Too few fields") {
    create_1by1_struct(n_fields - 1, fieldnames);
    REQUIRE_THROWS_AS(CMaterial(matlab_input), std::runtime_error);
  }
  SECTION("Too many fields") {
    const char *too_many_fields[] = {"Cax", "Cay", "Caz", "Cbx", "Cby",
                                     "Cbz", "Ccx", "Ccy", "Ccz", "extra"};
    create_1by1_struct(n_fields + 1, too_many_fields);
    REQUIRE_THROWS_AS(CMaterial(matlab_input), std::runtime_error);
  }
}

void CMaterialTest::test_incorrect_fieldname() {
  create_1by1_struct(n_fields, wrong_fieldnames);
  REQUIRE_THROWS_AS(CMaterial(matlab_input), std::runtime_error);
}

void CMaterialTest::test_correct_construction() {
  create_1by1_struct(n_fields, fieldnames);
  mxArray *elements[n_fields];
  for (int i = 0; i < n_fields; i++) {
    elements[i] = mxCreateNumericMatrix(1, n_numeric_elements, mxDOUBLE_CLASS, mxREAL);
    mxSetField(matlab_input, 0, fieldnames[i], elements[i]);
  }
  // construction should succeed
  CMaterial cm(matlab_input);
}

void DMaterialTest::test_incorrect_number_of_fields() {
  SECTION("Too few fields") {
    create_1by1_struct(n_fields - 1, fieldnames);
    REQUIRE_THROWS_AS(DMaterial(matlab_input), std::runtime_error);
  }
  SECTION("Too many fields") {
    const char *too_many_fields[] = {"Dax", "Day", "Daz", "Dbx", "Dby", "Dbz", "extra"};
    create_1by1_struct(n_fields + 1, too_many_fields);
    REQUIRE_THROWS_AS(DMaterial(matlab_input), std::runtime_error);
  }
}

void DMaterialTest::test_incorrect_fieldname() {
  create_1by1_struct(n_fields, wrong_fieldnames);
  REQUIRE_THROWS_AS(DMaterial(matlab_input), std::runtime_error);
}

void DMaterialTest::test_correct_construction() {
  create_1by1_struct(n_fields, fieldnames);
  mxArray *elements[n_fields];
  for (int i = 0; i < n_fields; i++) {
    elements[i] = mxCreateNumericMatrix(1, n_numeric_elements, mxDOUBLE_CLASS, mxREAL);
    mxSetField(matlab_input, 0, fieldnames[i], elements[i]);
  }
  // construction should succeed
  DMaterial dm(matlab_input);
}

TEST_CASE("CMaterial") { CMaterialTest().run_all_class_tests(); }

TEST_CASE("DMaterial") { DMaterialTest().run_all_class_tests(); }
