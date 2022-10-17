#include <complex>

#include "field.h"
#include "globals.h"

using namespace std;
using namespace tdms_math_constants;

void ElectricField::add_to_angular_norm(double f, int n, int Nt, SimulationParameters &params) {
  angular_norm += phasor_norm(f, n, params.omega_an, params.dt, Nt);
}

complex<double> ElectricField::phasor_norm(double f, int n, double omega, double dt, int Nt){
  return f
         * exp( fmod(omega*((double) (n+1))*dt, 2*DCPI) * IMAGINARY_UNIT)
         * 1./((double) Nt);
}
