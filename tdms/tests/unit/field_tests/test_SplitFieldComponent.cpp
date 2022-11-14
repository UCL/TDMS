/**
 * @file test_SplitFieldComponent.cpp
 * @brief Unit tests for SplitFieldComponent
 */
#include "field.h"

#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "unit_test_utils.h"

TEST_CASE("SplitFieldComponent: set component") {

  auto tmp = SplitFieldComponent();// initialisation != allocation
  REQUIRE(!tmp.has_elements());

  tmp.allocate(2, 2, 2);// 2x2x2 tensor
  tmp.zero();

  REQUIRE(is_close(tmp[0][1][0], 0.0));

  tmp[0][1][0] = 1.0;
  REQUIRE(is_close(tmp[0][1][0], 1.0));
}
