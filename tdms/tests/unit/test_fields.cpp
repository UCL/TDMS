#include <catch2/catch_test_macros.hpp>
#include <complex>
#include "field.h"

using namespace std;


template<typename T>
inline bool is_close(T a, T b){

  auto max_norm = max(abs(a), abs(b));

  if (max_norm < 1E-30){  // Prevent dividing by zero
    return true;
  }

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

TEST_CASE("Test that a split field can be allocated and zeroed"){

  auto field = ElectricSplitField(2, 1, 0);
  field.allocate_and_zero();

  REQUIRE(field.xy.has_elements());

  // NOTE: the allocated field is actually [I_tot+1, J_tot+1, K_tot+1] in size
  for (int k = 0; k < 1; k++){
    for (int j = 0; j < 2; j++){
      for (int i = 0; i < 3; i++){
        REQUIRE(is_close(field.xy[k][j][i], 0.0));
        REQUIRE(is_close(field.zy[k][j][i], 0.0));
        REQUIRE(is_close(field.yz[k][j][i], 0.0));
      }
    }
  }
}

TEST_CASE("Test that a split field can be constructed without allocation"){

  auto field = ElectricSplitField();
  REQUIRE(!field.xy.has_elements());
}

TEST_CASE("Test setting a component of an allocated vector field"){

  auto tmp = SplitFieldComponent();    // initialisation != allocation
  REQUIRE(!tmp.has_elements());

  tmp.allocate(2, 2, 2);  // 2x2x2 tensor
  tmp.zero();

  REQUIRE(is_close(tmp[0][1][0], 0.0));

  tmp[0][1][0] = 1.0;
  REQUIRE(is_close(tmp[0][1][0], 1.0));
}
