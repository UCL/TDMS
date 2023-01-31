/**
 * @file test_FullFieldSnapshot.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Unit tests for FullFieldSnapshot class
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "array_test_class.h"
#include "arrays.h"
#include "globals.h"
#include "unit_test_utils.h"

using namespace std;
using tdms_math_constants::IMAGINARY_UNIT;
using tdms_tests::is_close;

void FullFieldSnapshotTest::test_other_methods() {
  FullFieldSnapshot F;
  F.Ex = IMAGINARY_UNIT;
  F.Ey = 1.;
  F.Ez = 0.;
  F.Hx = 1. + IMAGINARY_UNIT;
  F.Hy = 1.;
  F.Hz = 0.;
  SECTION("multiply_{E,H}_by") {
    bool multiply_E_by_correct = true;
    // multiply E-field by imaginary unit
    F.multiply_E_by(IMAGINARY_UNIT);
    multiply_E_by_correct = multiply_E_by_correct && is_close(F.Ex, complex<double>(-1., 0.));
    multiply_E_by_correct = multiply_E_by_correct && is_close(F.Ey, IMAGINARY_UNIT);
    multiply_E_by_correct = multiply_E_by_correct && is_close(F.Ez, complex<double>(0., 0.));
    REQUIRE(multiply_E_by_correct);

    bool multiply_H_by_correct = true;
    // multiply H-field by the value of the Hx field
    F.multiply_H_by(F.Hx);
    multiply_H_by_correct = multiply_H_by_correct && is_close(F.Hx, 2. * IMAGINARY_UNIT);
    multiply_H_by_correct = multiply_H_by_correct && is_close(F.Hy, 1. + IMAGINARY_UNIT);
    multiply_H_by_correct = multiply_H_by_correct && is_close(F.Hz, complex<double>(0., 0.));
    REQUIRE(multiply_H_by_correct);
  }
}

TEST_CASE("FullFieldSnapshot") { FullFieldSnapshotTest().run_all_class_tests(); }
