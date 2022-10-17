#include "field.h"
#include "globals.h"
#include "interpolation_methods.h"

using namespace std;
using namespace tdms_math_constants;

double ElectricField::phase(int n, double omega, double dt){
  return omega * ((double) n + 1) * dt;
}

/*
You may notice several +1's appearing when the interpolation data (interp_data) arrays are constructed from the field values.
This is because E.x[k][j][i] is, by construction of the Grid class, not the value associated to Yee cell (i,j,k) but rather cell (i-1,j,k). Similarly Ey[k][j][i] is the value of cell (i,j-1,k) and Ez[k][j][i] of cell (i,j,k-1).
*/

complex<double> ElectricField::interpolate_to_centre_of(AxialDirection d, int i, int j, int k) {
  const InterpolationScheme *scheme;
  // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
  complex<double> interp_data[8];
  switch (d) {
    case X:
      // determine the interpolation scheme to use
      scheme = &(best_scheme(I_tot, i));
      // now fill the interpolation data
      // i - (scheme.number_of_datapoints_to_left) is the index of the Yee cell that plays the role of v0 in the interpolation
      for (int ind = scheme->first_nonzero_coeff; ind <= scheme->last_nonzero_coeff; ind++) {
        interp_data[ind] =
                real.x[k][j][i + 1 - scheme->number_of_datapoints_to_left + ind] +
                IMAGINARY_UNIT * imag.x[k][j][i + 1 - scheme->number_of_datapoints_to_left + ind];
      }
      break;
    case Y:
      // determine the interpolation scheme to use
      scheme = &(best_scheme(J_tot, j));

      // now fill the interpolation data
      // j - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
      for (int ind = scheme->first_nonzero_coeff; ind <= scheme->last_nonzero_coeff; ind++) {
        interp_data[ind] =
                real.y[k][j + 1 - scheme->number_of_datapoints_to_left + ind][i] +
                IMAGINARY_UNIT * imag.y[k][j + 1 - scheme->number_of_datapoints_to_left + ind][i];
      }
      break;
    case Z:
      // determine the interpolation scheme to use
      scheme = &(best_scheme(K_tot, k));

      // now fill the interpolation data
      // k - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
      for (int ind = scheme->first_nonzero_coeff; ind <= scheme->last_nonzero_coeff; ind++) {
        interp_data[ind] =
                real.z[k + 1 - scheme->number_of_datapoints_to_left + ind][j][i] +
                IMAGINARY_UNIT * imag.z[k + 1 - scheme->number_of_datapoints_to_left + ind][j][i];
      }
      break;
    default:
      throw runtime_error("Invalid axial direction selected for interpolation!\n");
      break;
  }
  // now run the interpolation scheme and place the result into the output
  return scheme->interpolate(interp_data);
}

double ElectricSplitField::interpolate_to_centre_of(AxialDirection d, int i, int j, int k) {
  const InterpolationScheme *scheme;
  // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
  double interp_data[8];

  switch (d) {
    case X:
      scheme = &(best_scheme(I_tot, i));
      // now fill the interpolation data
      // i - (scheme.number_of_datapoints_to_left) is the index of the Yee cell that plays the role of v0 in the interpolation
      for (int ind = scheme->first_nonzero_coeff; ind <= scheme->last_nonzero_coeff; ind++) {
        interp_data[ind] = xy[k][j][i + 1 - scheme->number_of_datapoints_to_left + ind] +
                           xz[k][j][i + 1 - scheme->number_of_datapoints_to_left + ind];
      }
      // now run the interpolation scheme and place the result into the output
      return scheme->interpolate(interp_data);
      break;
    case Y:
      scheme = &(best_scheme(J_tot, j));
      // now fill the interpolation data
      // j - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
      for (int ind = scheme->first_nonzero_coeff; ind <= scheme->last_nonzero_coeff; ind++) {
        interp_data[ind] = yx[k][j + 1 - scheme->number_of_datapoints_to_left + ind][i] +
                           yz[k][j + 1 - scheme->number_of_datapoints_to_left + ind][i];
      }
      // now run the interpolation scheme and place the result into the output
      return scheme->interpolate(interp_data);
      break;
    case Z:
      scheme = &(best_scheme(K_tot, k));
      // now fill the interpolation data
      // k - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
      for (int ind = scheme->first_nonzero_coeff; ind <= scheme->last_nonzero_coeff; ind++) {
        interp_data[ind] = zx[k + 1 - scheme->number_of_datapoints_to_left + ind][j][i] +
                           zy[k + 1 - scheme->number_of_datapoints_to_left + ind][j][i];
      }
      // now run the interpolation scheme and place the result into the output
      return scheme->interpolate(interp_data);
      break;
    default:
      throw runtime_error("Invalid axial direction selected for interpolation!\n");
      break;
  }
}
