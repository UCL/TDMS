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


void ElectricField::interpolate_across_range(mxArray **Ex_out, mxArray **Ey_out, mxArray **Ez_out) {
  interpolate_across_range(Ex_out, Ey_out, Ez_out, 0, I_tot, 0, J_tot, 0, K_tot);
}
void ElectricField::interpolate_across_range(mxArray **Ex_out, mxArray **Ey_out, mxArray **Ez_out, int i_lower_cutoff, int i_upper_cutoff, int j_lower_cutoff, int j_upper_cutoff, int k_lower_cutoff, int k_upper_cutoff) {
  // check we have 3 dimensional arrays - this is superflous thanks to our new classes?
  if ((int)mxGetNumberOfDimensions( (const mxArray *)real.x) < 3) {
    throw runtime_error("real.x is not 3D!\n");
  }
  // these are the dimensions of the field arrays. We assume they all have the same dimensions. We don't actually know this, unless we added an error check to the input array pointers when we read them in and used them to construct the field
  const int *indims = (int *)mxGetDimensions( (mxArray *)real.x );
  const int ndims = 3;

  // construct the output arrays
  int outdims[ndims] = {iu-il+1, jl-ju+1, ku-kl+1};
  if (outdims[1]<1) {
    // if simulation is 2D, allow one cell in the y-direction so the output is not NULL
    outdims[1] = 1;
  }
  *Ex_out = mxCreateNumericArray(ndims, (const mwSize *) outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *Ey_out = mxCreateNumericArray(ndims, (const mwSize *) outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *Ez_out = mxCreateNumericArray(ndims, (const mwSize *) outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  XYZTensor3D real_out, imag_out;
  real_out.x = cast_matlab_3D_array(mxGetPr(*Ex_out), outdims[0], outdims[1], outdims[2]);
  imag_out.x = cast_matlab_3D_array(mxGetPi(*Ex_out), outdims[0], outdims[1], outdims[2]);
  real_out.y = cast_matlab_3D_array(mxGetPr(*Ey_out), outdims[0], outdims[1], outdims[2]);
  imag_out.y = cast_matlab_3D_array(mxGetPi(*Ey_out), outdims[0], outdims[1], outdims[2]);
  real_out.z = cast_matlab_3D_array(mxGetPr(*Ez_out), outdims[0], outdims[1], outdims[2]);
  imag_out.z = cast_matlab_3D_array(mxGetPi(*Ez_out), outdims[0], outdims[1], outdims[2]);

  // now interpolate the fields, placing the values into the outputs
  if (J_tot==0) {
    // in a 2D simulation, interpolate across only two dimensions
    for(int i=i_lower_cutoff; i<=i_upper_cutoff; i++) {
      for(int k=k_lower_cutoff; k<=k_upper_cutoff; k++) {
        complex<double> x_at_centre = interpolate_x_to_centre(i, 0, k),
                        z_at_centre = interpolate_z_to_centre(i, 0, k);
        real_out.x[k][0][i] = x_at_centre.real();
        imag_out.x[k][0][i] = x_at_centre.imag();
        real_out.y[k][0][i] = real.y[k][0][i];
        imag_out.y[k][0][i] = imag.y[k][0][i];
        real_out.z[k][0][i] = z_at_centre.real();
        imag_out.z[k][0][i] = z_at_centre.imag();
      }
    }
  } else {
    for(int i=i_lower_cutoff; i<=i_upper_cutoff; i++) {
      for(int j=j_lower_cutoff; j<=j_upper_cutoff; j++) {
        for(int k=k_lower_cutoff; k<=k_upper_cutoff; k++) {
          complex<double> x_at_centre = interpolate_x_to_centre(i, j, k),
                          y_at_centre = interpolate_y_to_centre(i, j, k),
                          z_at_centre = interpolate_z_to_centre(i, j, k);
          real_out.x[k][j][i] = x_at_centre.real();
          imag_out.x[k][j][i] = x_at_centre.imag();
          real_out.y[k][j][i] = y_at_centre.real();
          imag_out.y[k][j][i] = y_at_centre.imag();
          real_out.z[k][j][i] = z_at_centre.real();
          imag_out.z[k][j][i] = z_at_centre.imag();
        }
      }
    }
  }
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
