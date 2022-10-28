#include "field.h"

#include "matlabio.h"
#include "globals.h"

using namespace std;
using namespace tdms_math_constants;

void Field::add_to_angular_norm(int n, int Nt, SimulationParameters &params) {
  angular_norm += phasor_norm(ft, n, params.omega_an, params.dt, Nt);
}



void Field::normalise_volume() {

  double norm_r = std::real(angular_norm);
  double norm_i = std::imag(angular_norm);
  double denom = norm_r * norm_r + norm_i * norm_i;

  for (char c : {'x', 'y', 'z'})
    for (int k = 0; k < K_tot; k++)
      for (int j = 0; j < J_tot; j++)
        for (int i = 0; i < I_tot; i++){

          auto temp_r = real[c][k][j][i];
          auto temp_i = imag[c][k][j][i];
          real[c][k][j][i] = (norm_r * temp_r + norm_i * temp_i) / denom;
          imag[c][k][j][i] = (norm_r * temp_i - norm_i * temp_r) / denom;
        }
}

Field::Field(int I_total, int J_total, int K_total) {
  I_tot = I_total;
  J_tot = J_total;
  K_tot = K_total;
}

void Field::allocate() {
  for (auto arr : {&real, &imag}) {
    arr->allocate(I_tot+1, J_tot+1, K_tot+1);
  }
};

void Field::zero() {
  for (auto &arr : {real, imag})
    for (char c : {'x', 'y', 'z'})
      for (int k = 0; k < K_tot; k++)
        for (int j = 0; j < J_tot; j++)
          for (int i = 0; i < I_tot; i++)
              arr[c][k][j][i] = 0.;
}

complex<double> Field::phasor_norm(double f, int n, double omega, double dt, int Nt){
  return f
         * exp( fmod(phase(n, omega, dt), 2*DCPI) * IMAGINARY_UNIT)
         * 1./((double) Nt);
}



Field::~Field() {

  for (auto &arr : {real, imag})
    for (char c : {'x', 'y', 'z'}){
      if (arr[c] != nullptr){
        free_cast_matlab_3D_array(arr[c], K_tot);
      }
    }
}

void Field::set_phasors(SplitField &F, int n, double omega, double dt, int Nt) {

  complex<double> subResult;
  double x_m, y_m, z_m;

  auto phaseTerm = exp(phase(n, omega, dt) * IMAGINARY_UNIT) * 1. / ((double) Nt);

#pragma omp parallel default(shared) private(x_m, y_m, z_m, subResult)
  {
#pragma omp for
    for (int k = kl; k <= ku; k++)
      for (int j = jl; j <= ju; j++)
        for (int i = il; i <= iu; i++) {

          x_m = F.xy[k][j][i] + F.xz[k][j][i];
          y_m = F.yx[k][j][i] + F.yz[k][j][i];
          z_m = F.zx[k][j][i] + F.zy[k][j][i];

          int di = i - il;
          int dj = j - jl;
          int dk = k - kl;

          subResult = x_m * phaseTerm;
          real.x[dk][dj][di] += std::real(subResult);
          imag.x[dk][dj][di] += std::imag(subResult);

          subResult = y_m * phaseTerm;
          real.y[dk][dj][di] += std::real(subResult);
          imag.y[dk][dj][di] += std::imag(subResult);

          subResult = z_m * phaseTerm;
          real.z[dk][dj][di] += std::real(subResult);
          imag.z[dk][dj][di] += std::imag(subResult);
        }
  }
}

void Field::set_values_from(Field &other) {

  if (other.I_tot != I_tot || other.J_tot != J_tot || other.K_tot != K_tot){
    throw runtime_error("Cannot set the values of a field with a different size");
  }

  for (auto &component : {"real", "imag"}) {

    auto this_component = are_equal(component, "real") ? real : imag;
    auto other_component = are_equal(component, "real") ? other.real : other.imag;

    for (char c : {'x', 'y', 'z'})
      for (int k = 0; k < K_tot; k++)
        for (int j = 0; j < J_tot; j++)
          for (int i = 0; i < I_tot; i++) {
            this_component[c][k][j][i] = other_component[c][k][j][i];
          }
  }
}
