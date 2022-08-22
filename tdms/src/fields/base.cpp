#include "matlabio.h"
#include "field.h"

using namespace std;


void Field::normalise_volume() {

  double norm_r = std::real(angular_norm);
  double norm_i = std::imag(angular_norm);
  double denom = norm_r * norm_r + norm_i * norm_i;

  for (int k = 0; k < K_tot; k++)
    for (int j = 0; j < J_tot; j++)
      for (int i = 0; i < I_tot; i++)
        for (char c : {'x', 'y', 'z'}){

          auto temp_r = real(c)[k][j][i];
          auto temp_i = imag(c)[k][j][i];
          real(c)[k][j][i] = (norm_r * temp_r + norm_i * temp_i) / denom;
          imag(c)[k][j][i] = (norm_r * temp_i - norm_i * temp_r) / denom;
        }
}

void Field::zero() {

  for (auto &arr : {real, imag})
    for (char c : {'x', 'y', 'z'})
      for (int k = 0; k < K_tot; k++)
        for (int j = 0; j < J_tot; j++)
          for (int i = 0; i < I_tot; i++)
              arr(c)[k][j][i] = 0.;
}

Field::~Field() {

  for (auto &arr : {real, imag})
    for (char c : {'x', 'y', 'z'})
      freeCastMatlab3DArray(arr(c), K_tot);
}
