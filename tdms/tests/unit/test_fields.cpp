/**
 * @file test_fields.cpp
 * @brief Test of the Field class and subclasses.
 */
#include <catch2/catch_test_macros.hpp>
#include <complex>
#include "field.h"
#include "globals.h"

using namespace std;

template<typename T>
inline bool is_close(T a, T b){
  return abs(a - b) / max(abs(a), abs(b)) < 1E-10;
}

/**
 * @brief Test that the electric field's angular normalisation is as expected.
 */
TEST_CASE("Test electric field angular norm addition") {

  double OMEGA = 0.7;
  double DT = 0.1;
  int N = 3;
  int N_T = 5;
  auto I = TDMS_MATH_CONSTANTS::IMAGINARY_UNIT;

  auto E = ElectricField();
  auto params = SimulationParameters();
  params.omega_an = OMEGA;
  params.dt = DT;

  // z = e^(iω(n+1)dt) / N_t
  auto expected = exp(OMEGA * ((double)N + 1) * DT * I) / ((double) N_T);

  E.add_to_angular_norm(1., N, N_T, params);
  REQUIRE(is_close(E.angular_norm, expected));
}

/**
 * @brief Test that the magnetic field's angular normalisation is as expected.
 */
TEST_CASE("Test magnetic field angular norm addition") {

  double OMEGA = 0.3;
  double DT = 0.2;
  int N = 7;
  int N_T = 5;
  auto I = TDMS_MATH_CONSTANTS::IMAGINARY_UNIT;

  auto H = MagneticField();
  auto params = SimulationParameters();
  params.omega_an = OMEGA;
  params.dt = DT;

  // z = e^(iω(n+1/2)dt) / N_t
  auto expected = exp(OMEGA * ((double)N + 0.5) * DT * I) / ((double) N_T);

  H.add_to_angular_norm(1., N, N_T, params);
  REQUIRE(is_close(H.angular_norm, expected));
}

