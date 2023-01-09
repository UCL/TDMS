#include "iterator_executor.h"

#include <complex>

#include "globals.h"

using tdms_math_constants::DCPI, tdms_math_constants::IMAGINARY_UNIT;
using namespace std;

bool Iterator_Executor::phasors_have_converged() {
  if ((dft_counter == Nsteps) && (params.run_mode == RunMode::complete) &&
      (params.source_mode == SourceMode::steadystate) && params.exphasorsvolume) {

    dft_counter = 0;

    double tol = E.normalised_difference(E_copy);
    if (tol < TOL) { return true; } //required accuracy obtained

    // if we get to here we need to setup the next iteration, as we have not converged
    spdlog::debug("Phasor convergence: {} (actual) > {} (required)", tol, TOL);
    E_copy.set_values_from(E);

    E.zero();
    H.zero();
    spdlog::debug("Zeroed the phasors");

    if (params.exphasorssurface) {
      surface_phasors.zero_surface_EH();
      spdlog::debug("Zeroed the surface components");
    }
  }
  // return false unless we exit early due to convergence
  return false;
}

void Iterator_Executor::extract_E_phasor_norm(int frequency_index, int Nt) {
  double omega = f_ex_vec[frequency_index] * 2 * DCPI;
  complex<double> phase_factor =
          exp(fmod(omega * ((double) (tind + 1)) * params.dt, 2 * DCPI) * IMAGINARY_UNIT) * 1. /
          ((double) Nt);
  E_norm[frequency_index] += E.ft * phase_factor;
}

void Iterator_Executor::extract_H_phasor_norm(int frequency_index, int Nt) {
  double omega = f_ex_vec[frequency_index] * 2 * DCPI;
  complex<double> phase_factor =
          exp(fmod(omega * ((double) tind + 0.5) * params.dt, 2 * DCPI) * IMAGINARY_UNIT) * 1. /
          ((double) Nt);
  H_norm[frequency_index] += H.ft * phase_factor;
}

void Iterator_Executor::extract_phasors_in_plane() {
  // short names for ease of writing
  int K1 = K0.index + 1;
  int Nt = params.Nt;

  complex<double> phaseTerm = 0., subResult = 0.;
  phaseTerm = fmod(params.omega_an * ((double) tind) * params.dt, 2 * DCPI);

  for (int j = 0; j < J_tot; j++)
    for (int i = 0; i < (I_tot + 1); i++) {
      //Eyz
      subResult = (E_s.yz[K1][j][i] + E_s.yx[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. /
                  ((double) Nt);

      iwave_lEy_Rbs[j][i] = iwave_lEy_Rbs[j][i] + real(subResult);
      iwave_lEy_Ibs[j][i] = iwave_lEy_Ibs[j][i] + imag(subResult);

      //Hxz
      subResult = (H_s.xz[K1 - 1][j][i] + H_s.xy[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. /
                  ((double) Nt);

      iwave_lHx_Rbs[j][i] = iwave_lHx_Rbs[j][i] + real(subResult);
      iwave_lHx_Ibs[j][i] = iwave_lHx_Ibs[j][i] + imag(subResult);
    }

  for (int j = 0; j < (J_tot + 1); j++)
    for (int i = 0; i < I_tot; i++) {
      //Exz
      subResult = (E_s.xz[K1][j][i] + E_s.xy[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. /
                  ((double) Nt);

      iwave_lEx_Rbs[j][i] = iwave_lEx_Rbs[j][i] + real(subResult);
      iwave_lEx_Ibs[j][i] = iwave_lEx_Ibs[j][i] + imag(subResult);

      //Hyz
      subResult = (H_s.yz[K1 - 1][j][i] + H_s.yx[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. /
                  ((double) Nt);

      iwave_lHy_Rbs[j][i] = iwave_lHy_Rbs[j][i] + real(subResult);
      iwave_lHy_Ibs[j][i] = iwave_lHy_Ibs[j][i] + imag(subResult);
    }
}
