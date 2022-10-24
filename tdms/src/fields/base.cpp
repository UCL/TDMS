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
    arr->allocate(I_tot, J_tot, K_tot);
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

void Field::interpolate_over_range(mxArray **x_out, mxArray **y_out, mxArray **z_out, int i_lower,
                                   int i_upper, int j_lower, int j_upper, int k_lower, int k_upper,
                                   Dimension mode) {
  // prepare output dimensions
  const int ndims = 3;
  int outdims[ndims] = {i_upper - i_lower + 1, j_upper - j_lower + 1, 1};
  if (mode == THREE) {
    outdims[2] = k_upper - k_lower + 1;
    if (outdims[1] < 1) {
      // full simulation (all components being computed) but in 2D - allow one cell in the y-direction to avoid NULL pointers
      outdims[1] = 1;
    }
  }

  // create the memory for the outputs
  *x_out = mxCreateNumericArray(ndims, (const mwSize *) outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *y_out = mxCreateNumericArray(ndims, (const mwSize *) outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *z_out = mxCreateNumericArray(ndims, (const mwSize *) outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  // allow memory to cast to our datatypes
  XYZTensor3D real_out, imag_out;
  real_out.x = cast_matlab_3D_array(mxGetPr(*x_out), outdims[0], outdims[1], outdims[2]);
  imag_out.x = cast_matlab_3D_array(mxGetPi(*x_out), outdims[0], outdims[1], outdims[2]);
  real_out.y = cast_matlab_3D_array(mxGetPr(*y_out), outdims[0], outdims[1], outdims[2]);
  imag_out.y = cast_matlab_3D_array(mxGetPi(*y_out), outdims[0], outdims[1], outdims[2]);
  real_out.z = cast_matlab_3D_array(mxGetPr(*z_out), outdims[0], outdims[1], outdims[2]);
  imag_out.z = cast_matlab_3D_array(mxGetPi(*z_out), outdims[0], outdims[1], outdims[2]);

  /* If we are not interpolating _all_ components, we need to use E- and H-field specific methods.
  Also (in full-field 3D simulations only) we need to check whether j_upper<j_lower as this indicates that the simulation is 2D and the y-axis doesn't exist!
  ALSO: note that there is a switch called within the for loop - this could be moved outside the loop, but this causes hideous code repetition.
  */
  if ((mode == THREE) && (j_upper < j_lower)) {
    // in a 2D simulation, interpolation can't occur in all 3 dimensions
    // beyond that however, interpolation is the same
    for (int i = i_lower; i <= i_upper; i++) {
      for (int k = k_lower; k <= k_upper; k++) {
        complex<double> x_at_centre = interpolate_to_centre_of(AxialDirection::X, i, 0, k),
                        z_at_centre = interpolate_to_centre_of(AxialDirection::Z, i, 0, k);
        real_out.x[k - k_lower][0][i - i_lower] = x_at_centre.real();
        imag_out.x[k - k_lower][0][i - i_lower] = x_at_centre.imag();
        // y interpolation doesn't take place, so use placeholder values
        real_out.y[k - k_lower][0][i - i_lower] = real.y[k][0][i];
        imag_out.y[k - k_lower][0][i - i_lower] = imag.y[k][0][i];
        real_out.z[k - k_lower][0][i - i_lower] = z_at_centre.real();
        imag_out.z[k - k_lower][0][i - i_lower] = z_at_centre.imag();
      }
    }
  } else {
    // 3D loop, and use the field-specific methods to ensure we interpolate components correctly
    for (int i = i_lower; i <= i_upper; i++) {
      for (int j = j_lower; j <= j_upper; j++) {
        for (int k = k_lower; k <= k_upper; k++) {
          complex<double> x_at_centre, y_at_centre, z_at_centre;
          switch (mode) {
            case THREE:
              // 3D interpolation is identical for both E and H fields
              x_at_centre = interpolate_to_centre_of(AxialDirection::X, i, j, k);
              y_at_centre = interpolate_to_centre_of(AxialDirection::Y, i, j, k);
              z_at_centre = interpolate_to_centre_of(AxialDirection::Z, i, j, k);
              break;
            case TE:
              interpolate_TE_components(i, j, k, &x_at_centre, &y_at_centre, &z_at_centre);
              break;
            case TM:
              interpolate_TM_components(i, j, k, &x_at_centre, &y_at_centre, &z_at_centre);
              break;
          }
          real_out.x[k - k_lower][j - j_lower][i - i_lower] = x_at_centre.real();
          imag_out.x[k - k_lower][j - j_lower][i - i_lower] = x_at_centre.imag();
          real_out.y[k - k_lower][j - j_lower][i - i_lower] = y_at_centre.real();
          imag_out.y[k - k_lower][j - j_lower][i - i_lower] = y_at_centre.imag();
          real_out.z[k - k_lower][j - j_lower][i - i_lower] = z_at_centre.real();
          imag_out.z[k - k_lower][j - j_lower][i - i_lower] = z_at_centre.imag();
        }
      }
    }
  }
}
void Field::interpolate_over_range(mxArray **x_out, mxArray **y_out, mxArray **z_out,
                            Dimension mode) {
  interpolate_over_range(x_out, y_out, z_out, 0, I_tot, 0, J_tot, 0, K_tot, mode);
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
