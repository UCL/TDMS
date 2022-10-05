#include <complex>

#include "field.h"
#include "globals.h"

using namespace std;
using namespace tdms_math_constants;

void MagneticField::add_to_angular_norm(double f, int n, int Nt, SimulationParameters &params) {
  angular_norm += phasor_norm(f, n, params.omega_an, params.dt, Nt);
}

complex<double> MagneticField::phasor_norm(double f, int n, double omega, double dt, int Nt){
  return f
         * exp( fmod(omega*((double) n + 0.5)*dt, 2*dcpi) * IMAGINARY_UNIT)
         * 1./((double) Nt);
}
