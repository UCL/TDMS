/**
 * @file test_Field.cpp
 * @brief Unit tests for the Field class and its subclasses (ElectricField, MagneticField)
 */
#include "field.h"

#include <complex>

#include <catch2/catch_test_macros.hpp>

#include "globals.h"
#include "unit_test_utils.h"

using namespace std;
using namespace tdms_math_constants;

const double tol = 1e-16;

TEST_CASE("Field: normalise_volume") {

  const int I_tot = 2, J_tot = 2, K_tot = 2;
  // create a field (note electric or magnetic should be sufficient here as we are testing base class method)
  ElectricField F(I_tot, J_tot, K_tot);
  F.allocate_and_zero();

  // set the angular norm to be 1 + i
  F.angular_norm = 1. + IMAGINARY_UNIT * 1.;
  // set components
  for (int i = 0; i < I_tot; i++) {
    for (int j = 0; j < J_tot; j++) {
      for (int k = 0; k < K_tot; k++) {
        // set x component to be 1 + i, so should normalise to 1.
        F.real.x[k][j][i] = 1.;
        F.imag.x[k][j][i] = 1.;
        // set y component to be 2, so should normalise to be 1 - i
        F.real.y[k][j][i] = 2.;
        F.imag.y[k][j][i] = 0.;
        // keep z components as 0
      }
    }
  }

  // now normalise
  F.normalise_volume();
  // check normalisation was correct
  bool x_normalised_correctly = true;
  bool y_normalised_correctly = true;
  bool z_normalised_correctly = true;
  for (int i = 0; i < I_tot; i++) {
    for (int j = 0; j < J_tot; j++) {
      for (int k = 0; k < K_tot; k++) {
        // set x component to be 1 + i, so should normalise to 1.
        x_normalised_correctly =
                (x_normalised_correctly &&
                 (abs(F.real.x[k][j][i] - 1. + IMAGINARY_UNIT * (F.imag.x[k][j][i] - 0.)) < tol));
        // set y component to be 2, so should normalise to be 1 - i
        y_normalised_correctly =
                (y_normalised_correctly &&
                 (abs(F.real.y[k][j][i] - 1. + IMAGINARY_UNIT * (F.imag.y[k][j][i] + 1.)) < tol));
        // keep z components as 0
        z_normalised_correctly =
                (z_normalised_correctly &&
                 (abs(F.real.z[k][j][i] - 0. + IMAGINARY_UNIT * (F.imag.z[k][j][i] - 0.)) < tol));
      }
    }
  }
  REQUIRE(x_normalised_correctly);
  REQUIRE(y_normalised_correctly);
  REQUIRE(z_normalised_correctly);
}

TEST_CASE("ElectricField: angular norm addition") {

  double OMEGA = 0.7;
  double DT = 0.1;
  int N = 3;
  int N_T = 5;

  auto E = ElectricField();
  E.ft = 1.0;
  auto params = SimulationParameters();
  params.omega_an = OMEGA;
  params.dt = DT;

  // z = e^(iω(n+1)dt) / N_t
  auto expected = exp(OMEGA * ((double) N + 1) * DT * IMAGINARY_UNIT) / ((double) N_T);

  E.add_to_angular_norm(N, N_T, params);
  REQUIRE(is_close(E.angular_norm, expected));
}

TEST_CASE("MagneticField: angular norm addition") {

  double OMEGA = 0.3;
  double DT = 0.2;
  int N = 7;
  int N_T = 5;

  auto H = MagneticField();
  H.ft = 1.0;
  auto params = SimulationParameters();
  params.omega_an = OMEGA;
  params.dt = DT;

  // z = e^(iω(n+1/2)dt) / N_t
  auto expected = exp(OMEGA * ((double) N + 0.5) * DT * IMAGINARY_UNIT) / ((double) N_T);

  H.add_to_angular_norm(N, N_T, params);
  REQUIRE(is_close(H.angular_norm, expected));
}
