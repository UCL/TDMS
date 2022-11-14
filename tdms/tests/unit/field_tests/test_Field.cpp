/**
 * @file test_Field.cpp
 * @brief Unit tests for the Field class and its subclasses (ElectricField, MagneticField)
 */
#include "field.h"

#include <complex>

#include <catch2/catch_test_macros.hpp>

#include "globals.h"

using namespace std;
using namespace tdms_math_constants;

template<typename T>
inline bool is_close(T a, T b) {

  auto max_norm = max(abs(a), abs(b));

  if (max_norm < 1E-30) {// Prevent dividing by zero
    return true;
  }

  return abs(a - b) / max(abs(a), abs(b)) < 1E-10;
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
