#include <complex>
#include "globals.h"
#include "field.h"

using namespace std;


complex<double> MagneticField::phasor_norm(double f, int n, double omega, double dt, int Nt){
  return f
         * exp( fmod(omega*((double) n + 0.5)*dt, 2*dcpi) * I)
         * 1./((double) Nt);
}
