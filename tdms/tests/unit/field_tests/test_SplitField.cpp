/**
 * @file test_SplitField.cpp
 * @brief Unit tests for the SplitField class and its subclasses
 * (ElectricSplitField, MagneticSplitField)
 */
#include "field.h"

#include <catch2/catch_test_macros.hpp>

#include "unit_test_utils.h"

using tdms_tests::is_close;

// NOTE: MagneticSplitField and ElectricSplitField inherit these methods from
// SplitField, so we can use either here

TEST_CASE("SplitField: allocate and zero") {

  ElectricSplitField field(2, 1, 0);
  field.allocate_and_zero();

  bool all_components_have_elements =
          field.xy.has_elements() && field.xz.has_elements() &&
          field.yx.has_elements() && field.yz.has_elements() &&
          field.zx.has_elements() && field.zy.has_elements();
  REQUIRE(all_components_have_elements);

  // NOTE: the allocated field is actually [I_tot+1, J_tot+1, K_tot+1] in size
  for (int k = 0; k < 1; k++) {
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 3; i++) {
        REQUIRE(is_close(field.xy(i, j, k), 0.0));
        REQUIRE(is_close(field.zy(i, j, k), 0.0));
        REQUIRE(is_close(field.yz(i, j, k), 0.0));
      }
    }
  }
}

TEST_CASE("SplitField: allocate without initialising") {

  auto field = ElectricSplitField();
  REQUIRE(!field.xy.has_elements());
}
