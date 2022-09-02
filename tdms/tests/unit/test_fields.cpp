#include <catch2/catch_test_macros.hpp>
#include <complex>
#include "field.h"

using namespace std;


template<typename T>
inline bool is_close(T a, T b){
  return abs(a - b) / max(abs(a), abs(b)) < 1E-10;
}

TEST_CASE("Test electric field angular norm addition") {

  double OMEGA = 0.7;
  double DT = 0.1;
  int N = 3;
  int N_T = 5;
  auto I = complex<double>(0., 1.);

  auto E = ElectricField();
  E.ft = 1.0;
  auto params = SimulationParameters();
  params.omega_an = OMEGA;
  params.dt = DT;

  // z = e^(iω(n+1)dt) / N_t
  auto expected = exp(OMEGA * ((double)N + 1) * DT * I) / ((double) N_T);

  E.add_to_angular_norm(N, N_T, params);
  REQUIRE(is_close(E.angular_norm, expected));
}

TEST_CASE("Test magnetic field angular norm addition") {

  double OMEGA = 0.3;
  double DT = 0.2;
  int N = 7;
  int N_T = 5;
  auto I = complex<double>(0., 1.);

  auto H = MagneticField();
  H.ft = 1.0;
  
  auto params = SimulationParameters();
  params.omega_an = OMEGA;
  params.dt = DT;

  // z = e^(iω(n+1/2)dt) / N_t
  auto expected = exp(OMEGA * ((double)N + 0.5) * DT * I) / ((double) N_T);

  H.add_to_angular_norm(N, N_T, params);
  REQUIRE(is_close(H.angular_norm, expected));
}

