/**
 * @file test_TDFieldExporter2D.cpp
 * @brief Unit tests for the TDFieldExporter2D class.
 */
#include "field.h"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("TDFieldExporter2D") {

  // create a 0-field for us to export
  const int I_tot = 8, J_tot = 1, K_tot = 8;
  ElectricSplitField Es(I_tot, J_tot, K_tot);
  Es.allocate_and_zero();

  // create our field exporter
  TDFieldExporter2D tdfe2d;
  tdfe2d.folder_name = ".";
  int stride;

  SECTION("Stride too small") {
    int nI = 4, nK = 4;
    stride = 1;
    REQUIRE_NOTHROW(tdfe2d.allocate(nI, nK));
    REQUIRE_THROWS_AS(tdfe2d.export_field(Es, stride, 0), std::runtime_error);
  }
  SECTION("Stride able to write out") {
    int nI = 8, nK = 8;
    stride = 2;
    REQUIRE_NOTHROW(tdfe2d.allocate(nI, nK));
    REQUIRE_NOTHROW(tdfe2d.export_field(Es, stride, 0));
  }
}
