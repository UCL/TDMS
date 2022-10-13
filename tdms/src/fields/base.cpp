#include "matlabio.h"
#include "field.h"
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

void Field::interpolate_across_range(mxArray **x_out, mxArray **y_out, mxArray **z_out) {
  interpolate_across_range(x_out, y_out, z_out, 0, I_tot, 0, J_tot, 0, K_tot);
}
void Field::interpolate_across_range(mxArray **x_out, mxArray **y_out, mxArray **z_out,
                                     int i_lower_cutoff, int i_upper_cutoff, int j_lower_cutoff,
                                     int j_upper_cutoff, int k_lower_cutoff, int k_upper_cutoff) {
  // check we have 3 dimensional arrays - this is superflous thanks to our new classes?
  if ((int) mxGetNumberOfDimensions((const mxArray *) real.x) < 3) {
    throw runtime_error("real.x is not 3D!\n");
  }
  // these are the dimensions of the field arrays. We assume they all have the same dimensions. We don't actually know this, unless we added an error check to the input array pointers when we read them in and used them to construct the field
  const int *indims = (int *) mxGetDimensions((mxArray *) real.x);
  const int ndims = 3;

  // construct the output arrays
  int outdims[ndims] = {iu - il + 1, jl - ju + 1, ku - kl + 1};
  if (outdims[1] < 1) {
    // if simulation is 2D, allow one cell in the y-direction so the output is not NULL
    outdims[1] = 1;
  }
  *x_out = mxCreateNumericArray(ndims, (const mwSize *) outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *y_out = mxCreateNumericArray(ndims, (const mwSize *) outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *z_out = mxCreateNumericArray(ndims, (const mwSize *) outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  XYZTensor3D real_out, imag_out;
  real_out.x = cast_matlab_3D_array(mxGetPr(*x_out), outdims[0], outdims[1], outdims[2]);
  imag_out.x = cast_matlab_3D_array(mxGetPi(*x_out), outdims[0], outdims[1], outdims[2]);
  real_out.y = cast_matlab_3D_array(mxGetPr(*y_out), outdims[0], outdims[1], outdims[2]);
  imag_out.y = cast_matlab_3D_array(mxGetPi(*y_out), outdims[0], outdims[1], outdims[2]);
  real_out.z = cast_matlab_3D_array(mxGetPr(*z_out), outdims[0], outdims[1], outdims[2]);
  imag_out.z = cast_matlab_3D_array(mxGetPi(*z_out), outdims[0], outdims[1], outdims[2]);

  // now interpolate the fields, placing the values into the outputs
  if (J_tot == 0) {
    // in a 2D simulation, interpolate across only two dimensions
    for (int i = i_lower_cutoff; i <= i_upper_cutoff; i++) {
      for (int k = k_lower_cutoff; k <= k_upper_cutoff; k++) {
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
    for (int i = i_lower_cutoff; i <= i_upper_cutoff; i++) {
      for (int j = j_lower_cutoff; j <= j_upper_cutoff; j++) {
        for (int k = k_lower_cutoff; k <= k_upper_cutoff; k++) {
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
