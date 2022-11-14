/**
 * @file test_SplitFieldComponent.cpp
 * @brief Unit tests for SplitFieldComponent
 */
#include "field.h"

#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

template<typename T>
inline bool is_close(T a, T b) {

  auto max_norm = max(abs(a), abs(b));

  if (max_norm < 1E-30) {// Prevent dividing by zero
    return true;
  }

  return abs(a - b) / max(abs(a), abs(b)) < 1E-10;
}

TEST_CASE("SplitFieldComponent: set component") {

  auto tmp = SplitFieldComponent();// initialisation != allocation
  REQUIRE(!tmp.has_elements());

  tmp.allocate(2, 2, 2);// 2x2x2 tensor
  tmp.zero();

  REQUIRE(is_close(tmp[0][1][0], 0.0));

  tmp[0][1][0] = 1.0;
  REQUIRE(is_close(tmp[0][1][0], 1.0));
}
