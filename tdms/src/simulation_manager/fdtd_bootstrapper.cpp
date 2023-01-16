#include "fdtd_bootstrapper.h"

#include "globals.h"

using tdms_math_constants::DCPI, tdms_math_constants::IMAGINARY_UNIT;
using namespace std;

void FDTDBootstrapper::allocate_memory(IJKTotal IJK_tot) {
    int I_tot = IJK_tot.I_tot(), J_tot = IJK_tot.J_tot();
  // x electric field source phasor - boot strapping
  Ex.allocate(I_tot, J_tot + 1);
  // y electric field source phasor - boot strapping
  Ey.allocate(I_tot + 1, J_tot);
  // x magnetic field source phasor - boot strapping
  Hx.allocate(I_tot + 1, J_tot);
  // y magnetic field source phasor - boot strapping
  Hy.allocate(I_tot, J_tot + 1);
}

void FDTDBootstrapper::extract_phasors_in_plane(ElectricSplitField &E_s, MagneticSplitField &H_s,
                                                IJKTotal IJK_tot, int K1, int tind,
                                                SimulationParameters params) {
  int Nt = params.Nt, I_tot = IJK_tot.I_tot(), J_tot = IJK_tot.J_tot();

  complex<double> phaseTerm = fmod(params.omega_an * ((double) tind) * params.dt, 2 * DCPI);
  complex<double> subResult = 0.;

  for (int j = 0; j < J_tot; j++) {
    for (int i = 0; i < (I_tot + 1); i++) {
      //Eyz
      subResult = (E_s.yz[K1][j][i] + E_s.yx[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. /
                  ((double) Nt);
      Ey[i][j] += subResult;
      //Hxz
      subResult = (H_s.xz[K1 - 1][j][i] + H_s.xy[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. /
                  ((double) Nt);
      Hx[i][j] += subResult;
    }
  }
  for (int j = 0; j < (J_tot + 1); j++)
    for (int i = 0; i < I_tot; i++) {
      //Exz
      subResult = (E_s.xz[K1][j][i] + E_s.xy[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. /
                  ((double) Nt);
      Ex[i][j] += subResult;
      //Hyz
      subResult = (H_s.yz[K1 - 1][j][i] + H_s.yx[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. /
                  ((double) Nt);
      Hy[i][j] += subResult;
    }
}
