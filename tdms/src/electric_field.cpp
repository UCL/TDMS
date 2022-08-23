#include <complex>
#include "globals.h"
#include "field.h"

using namespace std;


void ElectricField::add_to_angular_norm(double f, int n, int Nt, SimulationParameters &params) {
  angular_norm += phasor_norm(f, n, params.omega_an, params.dt, Nt);
}

complex<double> ElectricField::phasor_norm(double f, int n, double omega, double dt, int Nt){
  return f
         * exp( fmod(omega*((double) (n+1))*dt, 2*dcpi) * I)
         * 1./((double) Nt);
}
