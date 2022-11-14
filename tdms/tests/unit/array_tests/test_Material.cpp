/**
 * @file test_Material.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for the Material class and it's subclasses (CMaterial, DMaterial)
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

TEST_CASE("CMaterial") {
  SPDLOG_INFO("== CMaterial");
  // the constructor expects a struct of 9 fields
  const char *fieldnames[] = {"Cax", "Cay", "Caz", "Cbx", "Cby", "Cbz", "Ccx", "Ccy", "Ccz"};
  const int n_numeric_elements = 8;

  // providing incorrect names or sizes results in an error
  // bad fieldname
  const char *wrong_fieldnames[] = {"Cax", "Cay", "BAD", "Cbx", "Cby", "Cbz", "Ccx", "Ccy", "Ccz"};
  mxArray *wrong_struct = mxCreateStructMatrix(1, 1, 9, wrong_fieldnames);
  REQUIRE_THROWS_AS(CMaterial(wrong_struct), std::runtime_error);
  mxDestroyArray(wrong_struct);
  // bad number of fields
  const char *too_few_names[] = {"1", "2", "3", "4"};
  mxArray *too_few_fields_struct = mxCreateStructMatrix(1, 1, 4, too_few_names);
  REQUIRE_THROWS_AS(CMaterial(too_few_fields_struct), std::runtime_error);
  mxDestroyArray(too_few_fields_struct);

  // try and construct something meaningful
  mxArray *struct_9 = mxCreateStructMatrix(1, 1, 9, fieldnames);
  mxArray *elements9[9];
  for (int i = 0; i < 9; i++) {
    elements9[i] = mxCreateNumericMatrix(1, n_numeric_elements, mxDOUBLE_CLASS, mxREAL);
    mxSetField(struct_9, 0, fieldnames[i], elements9[i]);
  }
  // construction should succeed
  CMaterial cm(struct_9);
  mxDestroyArray(struct_9);
}
