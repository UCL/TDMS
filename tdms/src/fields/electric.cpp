#include <complex>
#include "globals.h"
#include "field.h"

using namespace std;


complex<double> ElectricField::phasor_norm(double f, int n, double omega, double dt, int Nt){
  return f
         * exp( fmod(omega*((double) (n+1))*dt, 2*dcpi) * I)
         * 1./((double) Nt);
}
