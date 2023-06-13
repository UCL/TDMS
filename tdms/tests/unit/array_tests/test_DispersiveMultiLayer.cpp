/**
 * @file test_DispersiveMultiLayer.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests for the DispersiveMultiLayer struct and its methods
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "array_test_class.h"
#include "arrays/dispersive_multilayer.h"
#include "unit_test_utils.h"

using namespace std;

/** @brief Test the methods associated to the DispersiveMultiLayer struct */
TEST_CASE("DispersiveMultiLayer") {
  spdlog::debug("[Unit] DispersiveMultiLayer");

  int N_ELEMENTS = 5;
  DispersiveMultiLayer dml;

  // Test the is_dispersive() method
  SECTION("is_dispersive()") {
    SECTION("non-empty vector") {
      /* Initialise the gamma component and set to 0 */
      dml.gamma.resize(N_ELEMENTS);
      // All entries being zero should result in a non-dispersive medium
      REQUIRE(!dml.is_dispersive());

      // Populate an element with the value 1.
      dml.gamma[N_ELEMENTS / 2] = 1.;
      /* A tolerance of 1.5 should still flag the dml as not dispersive */
      REQUIRE(!dml.is_dispersive(1.5));
      /* Yet a tolerance of 0.5 should flag it as dispersive */
      REQUIRE(dml.is_dispersive(0.5));
    }
    SECTION("empty vector") {
      /* When empty, is_dispersive() should not error and return false */
      REQUIRE(!dml.is_dispersive());
    }
  }
}
