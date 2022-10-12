#include "field.h"
#include "globals.h"
#include "interpolation_methods.h"

using namespace std;
using namespace tdms_math_constants;

double ElectricField::phase(int n, double omega, double dt){
  return omega * ((double) n + 1) * dt;
}

complex<double> ElectricField::interpolate_x_to_centre(int i, int j, int k) {
  // determine the interpolation scheme to use
  const InterpolationScheme &scheme = best_scheme(I_tot, i);
  // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
  complex<double> interp_data[8];

  // now fill the interpolation data
  // i - (scheme.number_of_datapoints_to_left) is the index of the Yee cell that plays the role of v0 in the interpolation
  for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
    interp_data[ind] = real.x[k][j][i - scheme.number_of_datapoints_to_left + ind] +
                       IMAGINARY_UNIT * imag.x[k][j][i - scheme.number_of_datapoints_to_left + ind];
  }
  // now run the interpolation scheme and place the result into the output
  return scheme.interpolate(interp_data);
}
complex<double> ElectricField::interpolate_y_to_centre(int i, int j, int k) {
  // determine the interpolation scheme to use
  const InterpolationScheme &scheme = best_scheme(J_tot, j);
  // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
  complex<double> interp_data[8];

  // now fill the interpolation data
  // j - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
  for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
    interp_data[ind] = real.y[k][j - scheme.number_of_datapoints_to_left + ind][i] +
                       IMAGINARY_UNIT * imag.y[k][j - scheme.number_of_datapoints_to_left + ind][i];
  }

  // now run the interpolation scheme and place the result into the output
  return scheme.interpolate(interp_data);
}
complex<double> ElectricField::interpolate_z_to_centre(int i, int j, int k) {
  // determine the interpolation scheme to use
  const InterpolationScheme &scheme = best_scheme(K_tot, k);
  // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
  complex<double> interp_data[8];

  // now fill the interpolation data
  // k - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
  for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
    interp_data[ind] = real.z[k - scheme.number_of_datapoints_to_left + ind][j][i] +
                       IMAGINARY_UNIT * imag.z[k - scheme.number_of_datapoints_to_left + ind][j][i];
  }

  // now run the interpolation scheme and place the result into the output
  return scheme.interpolate(interp_data);
}

double ElectricSplitField::interpolate_x_to_centre(int i, int j, int k) {
  // determine the interpolation scheme to use
  const InterpolationScheme &scheme = best_scheme(I_tot, i);
  // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
  double interp_data[8];

  // now fill the interpolation data
  // i - (scheme.number_of_datapoints_to_left) is the index of the Yee cell that plays the role of v0 in the interpolation
  for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
    interp_data[ind] = xy[k][j][i - scheme.number_of_datapoints_to_left + ind] +
                       xz[k][j][i - scheme.number_of_datapoints_to_left + ind];
  }

  // now run the interpolation scheme and place the result into the output
  return scheme.interpolate(interp_data);
}
double ElectricSplitField::interpolate_y_to_centre(int i, int j, int k) {
  // determine the interpolation scheme to use
  const InterpolationScheme &scheme = best_scheme(J_tot, j);
  // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
  double interp_data[8];

  // now fill the interpolation data
  // j - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
  for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
    interp_data[ind] = yx[k][j - scheme.number_of_datapoints_to_left + ind][i] +
                       yz[k][j - scheme.number_of_datapoints_to_left + ind][i];
  }

  // now run the interpolation scheme and place the result into the output
  return scheme.interpolate(interp_data);
}
double ElectricSplitField::interpolate_z_to_centre(int i, int j, int k) {
  // determine the interpolation scheme to use
  const InterpolationScheme &scheme = best_scheme(K_tot, k);
  // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
  double interp_data[8];

  // now fill the interpolation data
  // k - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
  for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
    interp_data[ind] = zx[k - scheme.number_of_datapoints_to_left + ind][j][i] +
                       zy[k - scheme.number_of_datapoints_to_left + ind][j][i];
  }

  // now run the interpolation scheme and place the result into the output
  return scheme.interpolate(interp_data);
}
